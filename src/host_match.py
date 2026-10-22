"""用 Genome 的 EC 集合匹配反应，并精确计算一条路线的最少 EC 覆盖。"""

import networkx as nx

from .logging_config import get_logger
from .pathway_validation import analyze_pathway_dependencies
from .schemas import Genome, PathWay


logger = get_logger("host_match")


def match_network(graph: nx.MultiDiGraph, genome: Genome) -> nx.MultiDiGraph:
    """原地重新匹配所有合并连接，返回同一网络；不修改原图及其来源记录。"""
    matched = unannotated = 0
    for _, _, data in graph.edges(data=True):
        connection = data["connection"]
        connection.is_endogenous = bool(connection.ec_numbers & genome.ec_numbers)
        matched += connection.is_endogenous
        unannotated += not connection.ec_numbers
    logger.info("宿主 EC 匹配完成：%s，%d 种 EC；%d 条反应匹配，%d 条未匹配，其中 %d 条无 EC 注释",
                genome.name or "未命名宿主", len(genome.ec_numbers), matched,
                graph.number_of_edges() - matched, unannotated)
    return graph


def minimum_ec_cover(choices) -> list[str] | None:
    """精确求每步 OR 选项的最少不同 EC；同大小解取编号字典序最小的一组。

    空路线返回 []；存在无 EC 的必需步骤返回 None。按步骤覆盖位图动态规划，
    保留每个覆盖状态的最少选择，不能用逐步选择 EC 的贪心算法替代。
    """
    choices = [set(options) for options in choices]
    if any(not options for options in choices):
        return None
    forced = {next(iter(options)) for options in choices if len(options) == 1}
    choices = [options for options in choices if not options & forced]
    coverage = {}
    for index, options in enumerate(choices):
        for ec in options:
            coverage[ec] = coverage.get(ec, 0) | (1 << index)
    best = {0: ()}
    for ec, mask in sorted(coverage.items()):
        for covered, selected in list(best.items()):
            new_mask = covered | mask
            candidate = (*selected, ec)
            previous = best.get(new_mask)
            if previous is None or (len(candidate), candidate) < (len(previous), previous):
                best[new_mask] = candidate
    return sorted(forced | set(best[(1 << len(choices)) - 1]))


def compare_pathway(pathway: PathWay, genome: Genome) -> dict:
    """只读比较一条 PathWay，不依赖 Connection 上可能属于其他宿主的匹配状态。

    多个 EC 按 OR 计算，部分 EC 只作字符串精确匹配，不作通配符。
    没有路线应传递 None 并报错；零步 PathWay 是有效的零补充路线。
    """
    if pathway is None:
        raise ValueError("没有候选路线，不能当作不需要补酶的零步路线")
    host_ecs, route_ecs, steps = genome.ec_numbers, set(), []
    for index, connection in enumerate(pathway.connections, 1):
        candidates = set(connection.ec_numbers)
        matched = candidates & host_ecs
        route_ecs.update(candidates)
        steps.append({
            "step": index, "source_id": connection.source.id, "target_id": connection.target.id,
            "reaction_ids": list(connection.reaction_ids),
            "substrate_ids": list(connection.substrate_ids), "product_ids": list(connection.product_ids),
            "ec_candidates": sorted(candidates), "matched_ecs": sorted(matched),
            "unmatched_ecs": sorted(candidates - host_ecs), "is_endogenous": bool(matched),
            "status": "matched" if matched else ("unmatched" if candidates else "unannotated"),
        })
    extra = minimum_ec_cover([step["ec_candidates"] for step in steps if not step["is_endogenous"]])
    return {
        "summary": {
            "host_non_pseudo_unique_ec_count": len(host_ecs),
            "reaction_steps": len(steps), "route_unique_ec_count": len(route_ecs),
            "matched_unique_ec_count": len(route_ecs & host_ecs),
            "unmatched_unique_ec_count": len(route_ecs - host_ecs),
            "matched_steps": sum(step["status"] == "matched" for step in steps),
            "unmatched_steps": sum(step["status"] == "unmatched" for step in steps),
            "unannotated_steps": sum(step["status"] == "unannotated" for step in steps),
            "minimum_extra_ec_count_under_or_assumption": len(extra) if extra is not None else None,
        },
        "matched_ecs": sorted(route_ecs & host_ecs), "unmatched_ecs": sorted(route_ecs - host_ecs),
        "one_minimum_extra_ec_set_under_or_assumption": extra, "steps": steps,
        "route_review": analyze_pathway_dependencies(pathway),
        "notes": [
            "按 Genome.ec_numbers 精确匹配非伪基因 EC；不推断未注释活性或扩展部分 EC 编号。",
            "每步多个 EC 按 OR 关系处理，同一种补充 EC 在整条路线中只计一次。",
            "无 EC 注释的步骤不能计算补充数，最少补酶搜索不使用此类边。",
            "结果是图示路线的 EC 编号覆盖候选，未验证共底物供应、表达或相同 EC 的底物适用性。",
        ],
    }
