"""统一导出前 N 条候选、中文报告、路线 SVG 和完整网络。"""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, is_dataclass
import json
from pathlib import Path

import networkx as nx

from .host_match import compare_pathway, match_network
from .kegg import get_compound_names
from .logging_config import get_logger
from .search import MinimumEcPathsResult, PathwayNetwork
from .schemas import Genome
from .svg_export import export_network_svg, route_colors_for


logger = get_logger("result_export")


class ResultExportError(RuntimeError):
    """部分图像导出失败；JSON 和报告保留，错误另列在 export_errors。"""


def json_value(value):
    """稳定转换集合、元组和数据对象，保留全部网络属性。"""
    if is_dataclass(value) and not isinstance(value, type):
        return json_value(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [json_value(item) for item in sorted(value, key=repr)]
    if isinstance(value, (tuple, list)):
        return [json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"无法序列化 {type(value).__name__}")


def check_output_directory(directory):
    """提前检查冲突，不删除或覆盖已有业务结果。"""
    directory = Path(directory).expanduser().resolve()
    if directory.exists():
        if not directory.is_dir():
            raise NotADirectoryError(f"输出路径不是目录：{directory}")
        conflicts = [directory / name for name in ("result.json", "report.md", "network.json", "network.svg")
                     if (directory / name).exists()]
        conflicts += sorted(directory.glob("pathway_*.svg"))
        if conflicts:
            raise FileExistsError(f"输出目录已有结果，请更换目录：{directory}（{conflicts[0].name}）")
    return directory


def _write_json(path, value, *, indent=2):
    path.write_text(json.dumps(json_value(value), ensure_ascii=False, indent=indent,
                               allow_nan=False) + "\n", encoding="utf-8")


def _raw_graph(graph):
    return {
        "metadata": graph.graph,
        "nodes": [{"id": cid, "entity": asdict(data["entity"]), "attributes": data}
                  for cid, data in sorted(graph.nodes(data=True))],
        "edges": [{"source": a, "target": b, "key": key,
                   "connection": asdict(data["connection"]),
                   "ko_ids": data.get("ko_ids", set()),
                   "diagram_entry_ids": data.get("diagram_entry_ids", set()), "attributes": data}
                  for a, b, key, data in sorted(graph.edges(keys=True, data=True),
                                                key=lambda edge: tuple(map(repr, edge[:3])))],
    }


def _candidate_record(network, genome, candidate, rank):
    pathway = candidate.pathway
    if pathway is None or len(pathway.connections) != len(candidate.trace):
        raise ValueError("候选必须包含有效路线，且连接数量与追溯步骤一致")
    matched = compare_pathway(pathway, genome)
    expected = matched["one_minimum_extra_ec_set_under_or_assumption"]
    if expected is None or set(expected) != candidate.extra_ecs:
        raise ValueError("候选需补 EC 与当前宿主的精确覆盖不一致，请重新搜索")
    if matched["route_review"]["circular_dependencies"]:
        raise ValueError("候选存在内部循环依赖，不能作为合格路线导出")
    start = pathway.connections[0].source.id if pathway.connections else network.target_id
    previous = start
    origins, edges = [], []
    for connection, event, step in zip(pathway.connections, candidate.trace, matched["steps"], strict=True):
        if (connection.source.id != previous or connection.source.id != event.source_id
                or connection.target.id != event.target_id):
            raise ValueError("候选路线不连续或与追溯信息不一致")
        data = network.graph[event.source_id][event.target_id][event.network_edge_key]
        if data["connection"] != connection:
            raise ValueError("候选连接与对应合并网络边不一致")
        previous = event.target_id
        origins.append(data.get("origins", []))
        edges.append((event.source_id, event.target_id, event.network_edge_key))
        step["selected_extra_ec"] = (None if step["is_endogenous"] else
                                      min(set(step["ec_candidates"]) & candidate.extra_ecs))
    if previous != network.target_id or start not in network.terminal_nodes:
        raise ValueError("候选必须从明确的 TCA 接点连接到目标")
    return {
        "rank": rank, "start_id": start, "target_id": network.target_id,
        "compound_ids": [start, *(c.target.id for c in pathway.connections)],
        "name": pathway.name, "reaction_steps": len(pathway.connections),
        "extra_ec_count": len(candidate.extra_ecs), "extra_ecs": sorted(candidate.extra_ecs),
        "connections": [asdict(connection) for connection in pathway.connections],
        "trace": [asdict(event) for event in candidate.trace], "origins": origins,
        "host_match": matched, "route_review": matched["route_review"],
        "svg_file": None,
    }, edges


def _route_graph(network, record, edges):
    # edge_subgraph 精确保留选中边；按节点取子图会错误带入未选择的平行反应。
    graph = deepcopy(network.graph.edge_subgraph(edges).copy())
    if not edges:
        cid = record["start_id"]
        graph.add_node(cid, **deepcopy(network.graph.nodes[cid]))
    graph.graph["terminal_nodes"] = frozenset({record["start_id"]})
    graph.graph["initial_nodes"] = (network.target_id,)
    graph.graph["candidate_rank"] = record["rank"]
    graph.graph["extra_ecs"] = record["extra_ecs"]
    graph.graph["terminal_origins"] = {
        record["start_id"]: deepcopy(network.graph.graph.get("terminal_origins", {}).get(record["start_id"], []))}
    involved = {origin["map_id"] for origins in record["origins"] for origin in origins}
    involved.update(item.get("map_id") for item in graph.graph["terminal_origins"][record["start_id"]])
    graph.graph["map_metadata"] = {
        mid: value for mid, value in graph.graph.get("map_metadata", {}).items() if mid in involved}
    return graph


def _cell(value):
    return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _compound_label(cid, names):
    name = names.get(cid)
    return f"{cid} ({name})" if name else cid


def _status(report):
    if not report["search_complete"]:
        return "前 N 条搜索尚未完成，当前候选未确认完整排名"
    if not report["data_is_complete"]:
        return "当前网络候选，数据不完整，不能确认搜索范围内全局最优"
    return "已确认当前搜索范围内的前 N 条候选" if report["candidates"] else "当前搜索范围内没有路线"


def _write_report(path, report):
    lines = [f"# {report['target_id']} 补酶候选报告", "", _status(report), "",
             f"宿主：{_cell(report['host']['name'] or '未命名宿主')}；EC 数：{report['host']['ec_count']}。", "",
             f"请求 {report['requested_count']} 条，实际返回 {report['returned_count']} 条；"
             f"前 N 条判定完成：{report['search_complete']}；全部路线穷尽：{report['exhausted']}。", "",
             "| 排名 | 起点 | 需补 EC 数 | 需补 EC | 反应步数 | 路线图 |",
             "| ---: | --- | ---: | --- | ---: | --- |"]
    for record in report["candidates"]:
        figure = f"[查看]({record['svg_file']})" if record["svg_file"] else "未生成"
        start = _cell(_compound_label(record["start_id"], record["compound_names"]))
        lines.append(f"| {record['rank']} | {start} | {record['extra_ec_count']} | "
                     f"{', '.join(record['extra_ecs']) or '无需补充'} | {record['reaction_steps']} | {figure} |")
    if report["network_svg_file"]:
        lines += ["", f"[查看完整网络 SVG]({report['network_svg_file']})。节点名称和反应 EC 在悬停时显示。"]
    if report["route_colors"]:
        lines += ["", "大图候选颜色：" + "；".join(
            f"候选 {rank} `{color}`" for rank, color in report["route_colors"].items()) + "。",
            "共享反应使用多色并列线；独立候选图仍以绿／红表示宿主匹配。"]
    if not report["candidates"]:
        lines += ["", "数据不完整，暂未找到路线。" if not report["is_complete"] else "没有符合条件的路线。"]
    for record in report["candidates"]:
        review = record["route_review"]
        background = sorted({cid for item in review["unverified_substrates"] for cid in item["compound_ids"]})
        names = record["compound_names"]
        lines += ["", f"## 候选 {record['rank']}", "",
                  " → ".join(_cell(_compound_label(cid, names)) for cid in record["compound_ids"]), "",
                  "内部循环检查：未发现内部循环依赖。",
                  f"其他原料供应未核查：{', '.join(background) or '图示记录未列出其他原料'}。", "",
                  "| 步骤 | 化合物 | R 编号 | EC 候选 | 宿主匹配 EC | 本步采用的补充 EC | 来源 map |",
                  "| ---: | --- | --- | --- | --- | --- | --- |"]
        for step, origins in zip(record["host_match"]["steps"], record["origins"], strict=True):
            cells = [step["step"], f"{_compound_label(step['source_id'], names)} → {_compound_label(step['target_id'], names)}",
                     ", ".join(step["reaction_ids"]), ", ".join(step["ec_candidates"]),
                     ", ".join(step["matched_ecs"]) or "未匹配", step["selected_extra_ec"] or "无需补充",
                     ", ".join(sorted({item["map_id"] for item in origins}))]
            lines.append("| " + " | ".join(_cell(cell) for cell in cells) + " |")
    lines += ["", "## 数据异常与范围", "", report["scope"], "", "失败分支：", ""]
    lines += [f"- {_cell(item['map_id'] or '初始查询')} / {_cell(item['compound_id'])} / "
              f"{_cell(item['stage'])}：{_cell(item['reason'])}" for item in report["failures"]] or ["无。"]
    lines += ["", "404 跳过（已排除在范围外）：", ""]
    lines += [f"- {mid}：{_cell(reason)}" for mid, reason in report["skipped_maps"].items()] or ["无。"]
    lines += ["", "导出错误：", ""]
    lines += [f"- {_cell(item['file'])}：{_cell(item['reason'])}" for item in report["export_errors"]] or ["无。"]
    lines += ["", "化合物名称查询提示（不影响搜索完整性）：", ""]
    lines += [f"- {item['compound_id']}：{_cell(item['reason'])}" for item in report["name_lookup_warnings"]] or ["无。"]
    lines += ["", *report["notes"], "", "完整网络与原图信息保存在 network.json。"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def export_results(network: PathwayNetwork, genome: Genome, result: MinimumEcPathsResult, output_dir,
                   *, cache_dir=None, refresh=False, include_network_svg=False):
    """导出前补齐候选名称；部分 SVG 失败时保留结果并抛 ResultExportError。

    名称查询复用 cache_dir 和 refresh，失败仅产生 name_lookup_warnings，
    不改变搜索完整性。成功名称写入已有 Entity.name。include_network_svg=True
    时查询全网名称，并在候选图之后导出按候选着色的完整 network.svg。
    """
    directory = check_output_directory(output_dir)
    match_network(network.graph, genome)
    records, paths = [], []
    for rank, candidate in enumerate(result.candidates, 1):
        record, edges = _candidate_record(network, genome, candidate, rank)
        records.append(record)
        paths.append(edges)
    logger.info("正在导出 %d 条候选的报告和路线图%s", len(records),
                "，并生成完整网络 SVG" if include_network_svg else "")
    compound_ids = sorted(network.graph) if include_network_svg else sorted(
        {cid for record in records for cid in record["compound_ids"]})
    names = get_compound_names(compound_ids, cache_dir=cache_dir, refresh=refresh) if compound_ids else {}
    name_warnings = []
    for cid in compound_ids:
        entity = network.graph.nodes[cid]["entity"]
        if cid in names:
            entity.name = names[cid]
        else:
            reason = "未获取到新名称；保留已有名称" if entity.name else "未获取到名称；仅显示 C 编号"
            name_warnings.append({"compound_id": cid, "reason": reason})
    for record, candidate in zip(records, result.candidates, strict=True):
        record["compound_names"] = {cid: network.graph.nodes[cid]["entity"].name for cid in record["compound_ids"]}
        record["connections"] = [asdict(connection) for connection in candidate.pathway.connections]
        for step in record["host_match"]["steps"]:
            step["source_name"] = record["compound_names"][step["source_id"]]
            step["target_name"] = record["compound_names"][step["target_id"]]
    report = {
        "schema": "pathfiller.minimum_ec_paths.v1", "schema_version": 1,
        "target_id": network.target_id,
        "host": {"name": genome.name, "file_path": genome.file_path, "ec_count": len(genome.ec_numbers)},
        "requested_count": result.requested_count, "returned_count": len(records),
        "data_is_complete": result.data_is_complete, "search_complete": result.search_complete,
        "exhausted": result.exhausted, "is_complete": result.is_complete,
        "objective": ["distinct_extra_ecs", "reaction_steps", "extra_ecs_lexicographic"],
        "map_count": len(network.graphs), "node_count": len(network.graph),
        "edge_count": network.graph.number_of_edges(), "candidates": records,
        "failures": [asdict(item) for item in result.failures], "skipped_maps": result.skipped_maps,
        "export_errors": [], "name_lookup_warnings": name_warnings, "network_file": "network.json",
        "network_svg_file": None,
        "route_colors": route_colors_for(range(1, len(records) + 1)) if include_network_svg else {},
        "scope": "按 C 编号合并的详细代谢图及明确 TCA 接点；排除 HTTP 404 图和路线内部循环依赖，仅搜索有 EC 注释的边。",
        "notes": ["同一步多个 EC 按 OR 选项处理，同一补充 EC 用于多步只计一次。",
                  "每条反应路线只输出一组最优 EC 覆盖，不表示所需基因数。",
                  "前 N 条已确定不表示全部路线穷尽；数据异常时只报告当前网络候选。",
                  "未验证其他共底物供应、计量、碳流、表达及酶的底物适用性；未发现内部循环不等于有效合成方案。"],
    }
    network_data = {
        "schema_version": 2, "representation": "compound", "target_id": network.target_id,
        "graphs": {mid: _raw_graph(graph) for mid, graph in sorted(network.graphs.items())},
        "graph": network.graph.graph,
        "nodes": [{"id": cid, "attributes": data} for cid, data in sorted(network.graph.nodes(data=True))],
        "edges": [{"source": a, "target": b, "key": key, "attributes": data}
                  for a, b, key, data in sorted(network.graph.edges(keys=True, data=True),
                                                key=lambda edge: tuple(map(repr, edge[:3])))],
        "initial_nodes": network.initial_nodes, "terminal_nodes": sorted(network.terminal_nodes),
    }
    directory.mkdir(parents=True, exist_ok=True)
    artifacts = {name: directory / name for name in ("result.json", "report.md", "network.json")}
    _write_json(artifacts["network.json"], network_data, indent=None)
    _write_json(artifacts["result.json"], report)
    _write_report(artifacts["report.md"], report)
    for record, edges in zip(records, paths, strict=True):
        name = f"pathway_{record['rank']:03d}.svg"
        try:
            graph = _route_graph(network, record, edges)
            export_network_svg(graph, directory / name,
                               title=f"Candidate {record['rank']}: {record['start_id']} to {network.target_id}")
            record["svg_file"] = name
            artifacts[name] = directory / name
        except (RuntimeError, OSError, ValueError) as error:
            report["export_errors"].append({"file": name, "reason": str(error)})
            logger.error("候选 %d SVG 导出失败：%s", record["rank"], error)
    if include_network_svg:
        logger.info("正在绘制完整网络 SVG：%d 个化合物、%d 条反应", len(network.graph),
                    network.graph.number_of_edges())
        try:
            graph = network.graph.copy()
            graph.graph["target_id"] = network.target_id
            export_network_svg(graph, directory / "network.svg", overview=True,
                               highlighted_paths={record["rank"]: edges for record, edges in
                                                  zip(records, paths, strict=True)},
                               title=f"Complete compound network: {network.target_id}")
            report["network_svg_file"] = "network.svg"
            artifacts["network.svg"] = directory / "network.svg"
        except (RuntimeError, OSError, ValueError) as error:
            report["export_errors"].append({"file": "network.svg", "reason": str(error)})
            logger.error("完整网络 SVG 导出失败：%s", error)
    _write_json(artifacts["result.json"], report)
    _write_report(artifacts["report.md"], report)
    logger.info("结果导出完成：%d 条候选，%d 个导出错误", len(records), len(report["export_errors"]))
    logger.debug("输出目录：%s", directory)
    if report["export_errors"]:
        raise ResultExportError(f"{len(report['export_errors'])} 张 SVG 导出失败；JSON 和报告已保存在 {directory}")
    return artifacts
