"""按化合物合并详细代谢图，在统一网络中查找 TCA 路线并枚举无环候选。"""

from collections import deque
from collections.abc import Hashable, Iterator
from dataclasses import dataclass, field
from heapq import heappop, heappush
from itertools import count
import json
from pathlib import Path
import re
from time import perf_counter
from typing import Literal

import networkx as nx

from . import kegg
from .logging_config import get_logger
from .host_match import match_network
from .pathway_validation import analyze_pathway_dependencies
from .schemas import Connection, Entity, Genome, PathWay


logger = get_logger("search")
_STAGES = {"association": "查询关联图", "classification": "核查通路分类", "download": "下载数据",
           "parse": "解析数据", "compound": "核对链接化合物", "graph": "校验网络"}


@dataclass(frozen=True)
class SearchFailure:
    map_id: str | None
    compound_id: str
    stage: str
    reason: str


@dataclass(frozen=True)
class ReactionStep:
    """一步合并反应；原字段指代表性来源，network_edge_key 指向合并边。"""

    map_id: str
    source_id: str
    target_id: str
    edge_key: Hashable
    diagram_entry_ids: tuple[str, ...]
    network_edge_key: Hashable | None = None
    kind: Literal["reaction"] = field(default="reaction", init=False)
    cost: int = field(default=1, init=False)


@dataclass(frozen=True)
class MapSwitch:
    """旧报告的兼容类型；新的合并网络搜索不产生此事件。"""

    source_map: str
    target_map: str
    compound_id: str
    kind: Literal["map_switch"] = field(default="map_switch", init=False)
    cost: int = field(default=0, init=False)


@dataclass
class SearchResult:
    pathway: PathWay | None
    is_complete: bool
    failures: list[SearchFailure]
    trace: list[ReactionStep | MapSwitch]
    graphs: dict[str, nx.MultiDiGraph]
    skipped_maps: dict[str, str] = field(default_factory=dict)


@dataclass
class MinimumEcResult(SearchResult):
    """最少不同补充 EC 的候选；无路线看 pathway，extra_ecs 为空也可能是零补充路线。"""

    extra_ecs: set[str] = field(default_factory=set)


@dataclass
class MinimumEcPathsResult:
    """前 N 条不同反应路线；前 N 条已确定与全部路线已穷尽分别记录。"""

    candidates: list[MinimumEcResult]
    requested_count: int
    data_is_complete: bool
    search_complete: bool
    exhausted: bool
    is_complete: bool
    failures: list[SearchFailure]
    graphs: dict[str, nx.MultiDiGraph]
    skipped_maps: dict[str, str]


_MAP_ID = re.compile(r"map[0-9]{5}")
_TCA = "map00020"


def _special_map(map_id: str) -> bool:
    # KEGG 的结构图、总图和概览图不能提供本搜索所需的详细反应方向。
    # 01040、01056、01057、01059 是具有 KO 展开的普通代谢图，不能一并排除。
    structure = map_id.startswith("map010") and map_id not in {
        "map01040", "map01056", "map01057", "map01059",
    }
    return structure or map_id[3:6] in {"011", "012", "013"} or map_id.startswith("map07")


def _fields(text: str) -> dict[str, str]:
    """读取 KEGG 固定宽度条目，包括字段的续行。"""
    fields: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        label = line[:12].strip()
        if label:
            current = label
            fields.setdefault(current, [])
        if current is not None:
            fields[current].append(line[12:].strip())
    return {label: " ".join(parts) for label, parts in fields.items()}


def _metabolic_classification(map_id: str, cache: Path, refresh: bool) -> bool:
    path = cache / f"{map_id}.pathway.txt"
    cached = path.exists() and not refresh
    logger.debug("%s 分类信息：%s", map_id, "读取缓存" if cached else "从 KEGG 获取")
    text = path.read_text(encoding="utf-8") if cached else kegg._request_text(f"get/{map_id}")
    fields = _fields(text)
    if fields.get("ENTRY", "").split()[:1] != [map_id] or not fields.get("CLASS"):
        raise ValueError(f"{map_id} 的条目编号不匹配或缺少 CLASS 分类")
    if not cached:
        cache.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return fields["CLASS"].split(";", 1)[0].strip() == "Metabolism"


def _initial_maps(target: str, cache: Path, refresh: bool) -> list[str]:
    path = cache / f"{target}.maps.json"
    cached = path.exists() and not refresh
    if cached:
        logger.debug("读取 %s 的关联 map 缓存：%s", target, path.name)
        record = json.loads(path.read_text(encoding="utf-8"))
    else:
        logger.debug("获取 %s 的关联 map 列表", target)
        record = {
            "target_id": target,
            "maps": kegg.get_maps_by_compound(target, refresh=refresh),
        }
    if (
        not isinstance(record, dict) or record.get("target_id") != target
        or not isinstance(record.get("maps"), list)
        or any(not isinstance(mid, str) or _MAP_ID.fullmatch(mid) is None for mid in record["maps"])
    ):
        raise ValueError(f"{target} 的关联 map 列表格式错误或缓存编号不匹配")
    if not cached:
        cache.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    maps = sorted(set(record["maps"]))
    logger.debug("初始关联图共 %d 张：%s", len(maps), "、".join(maps) or "无")
    return maps


class _MapLoader:
    def __init__(self, target: str, cache: Path, refresh: bool):
        self.target, self.cache, self.refresh = target, cache, refresh
        self.entities: dict[str, Entity] = {}
        self.graphs: dict[str, nx.MultiDiGraph] = {}
        self.failures: list[SearchFailure] = []
        self.failed_maps: dict[str, tuple[str, str]] = {}
        self.skipped_maps: dict[str, str] = {}
        self.eligible: dict[str, bool] = {}

    def fail(self, map_id, compound, stage, reason):
        failure = SearchFailure(map_id, compound, stage, str(reason))
        if failure not in self.failures:
            self.failures.append(failure)
            logger.debug("%s / %s：%s失败：%s", map_id or "初始查询", compound,
                           _STAGES.get(stage, stage), reason)

    def map_error(self, map_id, compound, stage, error):
        if isinstance(error, kegg.KeggRequestError) and error.status_code == 404:
            self.skipped_maps[map_id] = f"{stage}: {error}"
            logger.debug("跳过 %s：%s返回 HTTP 404；该图排除在搜索范围外，不影响数据完整性",
                        map_id, _STAGES.get(stage, stage))
        else:
            self.failed_maps[map_id] = (stage, str(error))
            self.fail(map_id, compound, stage, error)

    def load(self, map_id, compound):
        if map_id in self.skipped_maps:
            logger.debug("%s 本轮已因 HTTP 404 跳过，不再请求", map_id)
            return None
        if map_id in self.failed_maps:
            logger.debug("%s 本轮已失败，不再请求", map_id)
            stage, reason = self.failed_maps[map_id]
            self.fail(map_id, compound, stage, reason)
            return None
        if map_id not in self.eligible:
            try:
                self.eligible[map_id] = not _special_map(map_id) and _metabolic_classification(
                    map_id, self.cache, self.refresh,
                )
            except Exception as error:
                self.map_error(map_id, compound, "classification", error)
                return None
        if not self.eligible[map_id]:
            logger.debug("筛除 %s：%s", map_id, "总图、概览图或结构图" if _special_map(map_id)
                        else "不属于普通详细代谢图")
            return None
        if map_id not in self.graphs:
            logger.debug("加载 %s，核对接入化合物 %s", map_id, compound)
            try:
                graph = kegg.parse_map(
                    map_id, entities=self.entities, cache_dir=self.cache, refresh=self.refresh,
                )
                if not isinstance(graph, nx.MultiDiGraph) or graph.graph.get("map_id") != map_id:
                    raise ValueError(f"{map_id} 返回了类型或编号不匹配的网络")
                self.graphs[map_id] = graph
            except Exception as error:
                stage = "download" if isinstance(error, RuntimeError) else "parse"
                self.map_error(map_id, compound, stage, error)
                return None
        graph = self.graphs[map_id]
        if compound not in graph:
            self.fail(map_id, compound, "compound", "链接指向的网络中没有对应化合物")
            return None
        return graph



def _edge_order(item):
    key, data = item
    return (tuple(sorted(data["connection"].reaction_ids)), repr(key))


def search_shortest_path(target_id: str, *, cache_dir: str | Path | None = None,
                         refresh: bool = False) -> SearchResult:
    """在与全部路线枚举相同的化合物网络中找最少反应步数的候选。

    先完成可达数据的组图闭包；所有反应成本 1，不计算宿主补充数。
    相同 C 编号自然连通，不产生 MapSwitch。失败/404/缓存规则与组图一致。
    """
    logger.info("开始搜索到 %s 的最少反应步数路线", target_id)
    network = build_pathway_network(target_id, cache_dir=cache_dir, refresh=refresh)
    target = network.target_id
    pathway, trace, _ = _find_acyclic_path(network, host_aware=False)
    if pathway is not None:
        start = pathway.connections[0].source.id if pathway.connections else target
        logger.info("%s：%s → %s，共 %d 步反应", "找到范围内最短候选" if network.is_complete else
                    "找到当前候选（数据不完整，尚不能确认最短）", start, target, len(pathway.connections))
    elif network.is_complete:
        logger.info("搜索结束：当前范围内没有连接 TCA 的路线")
    else:
        logger.warning("搜索结束：暂未找到连接 TCA 的路线，数据不完整，尚无法确定是否存在路线")
    return SearchResult(pathway, network.is_complete, network.failures, trace, network.graphs, network.skipped_maps)


@dataclass
class CandidatePath:
    """连接及反应追踪按合成方向排列，不包含换图事件。"""

    pathway: PathWay
    trace: list[ReactionStep]


@dataclass
class PathwayNetwork:
    """按 C 编号合并的代谢网络；maps 只用于下载与来源追溯。

    graph 的节点直接是 C 编号，所有边为 weight=1 的反应。
    graphs 保存未修改的原始图；合并边 origins 保留全部原图来源。
    is_complete 表示组图数据完整性，枚举完整性查看 iter_paths() 返回对象。
    """

    target_id: str
    graph: nx.MultiDiGraph
    graphs: dict[str, nx.MultiDiGraph]
    initial_nodes: tuple[str, ...]
    terminal_nodes: frozenset[str]
    is_complete: bool
    failures: list[SearchFailure]
    skipped_maps: dict[str, str]

    def iter_paths(self) -> "PathwayEnumeration":
        """逐条生成全部不重复化合物的边路径，不设隐式截断。"""
        return PathwayEnumeration(self)


def find_minimum_ec_path(network: PathwayNetwork, genome: Genome) -> MinimumEcResult:
    """兼容单路线调用；与批量接口使用同一个有序搜索内核。"""
    result = find_minimum_ec_paths(network, genome, top_n=1)
    if result.candidates:
        return result.candidates[0]
    return MinimumEcResult(None, result.is_complete, result.failures, [],
                           result.graphs, result.skipped_maps, set())


def find_minimum_ec_paths(network: PathwayNetwork, genome: Genome, *, top_n: int = 3) -> MinimumEcPathsResult:
    """直接搜索前 N 条反应路线，同一路线只保留最优 EC 覆盖方案。

    排名为不同补充 EC 数、反应步数、补充 EC 字典序，再按稳定遍历次序。
    top_n 是输出数量，不是枚举预览上限；无 EC 边不参加搜索。
    """
    if type(top_n) is not int or top_n < 1:
        raise ValueError("top_n 必须是正整数")
    started = perf_counter()
    match_network(network.graph, genome)
    logger.info("宿主匹配完成，开始搜索前 %d 条补酶候选：%d 个 TCA 接点 → %s",
                top_n, len(network.terminal_nodes), network.target_id)
    routes, exhausted = _find_acyclic_paths(network, host_aware=True, top_n=top_n)
    candidates = [MinimumEcResult(pathway, network.is_complete, list(network.failures), trace,
                                  network.graphs, dict(network.skipped_maps), ecs)
                  for pathway, trace, ecs in routes]
    logger.info("补酶搜索完成：请求 %d 条，返回 %d 条；全部穷尽：%s；数据%s；用时 %.2f 秒",
                top_n, len(candidates), "是" if exhausted else "否",
                "完整" if network.is_complete else "不完整，结果为当前网络候选", perf_counter() - started)
    return MinimumEcPathsResult(candidates, top_n, network.is_complete, True, exhausted,
                                network.is_complete, list(network.failures), network.graphs,
                                dict(network.skipped_maps))


def _find_acyclic_path(network: PathwayNetwork, *, host_aware: bool):
    routes, _ = _find_acyclic_paths(network, host_aware=host_aware, top_n=1)
    return routes[0] if routes else (None, [], set())


def _find_acyclic_paths(network: PathwayNetwork, *, host_aware: bool, top_n: int):
    """反向 A* 的有序前缀；状态保留边路径，避免等价状态吞掉其他路线。

    同一边序列不同 EC 选择中，第一次出队的终点就是该路线最优覆盖。
    只在输出时去掉重复路线；其状态仍可扩展到更上游的 TCA 接点。
    """
    graph, target = network.graph, network.target_id
    if target not in graph:
        return [], True
    eligible = nx.subgraph_view(graph, filter_edge=lambda a, b, key:
                               not host_aware or bool(graph[a][b][key]["connection"].ec_numbers))
    allowed = nx.ancestors(eligible, target) | {target}
    sources = network.terminal_nodes & allowed
    if not sources:
        return [], True
    incoming = {
        node: sorted((edge for edge in graph.in_edges(node, keys=True, data=True)
                      if edge[0] in allowed and (not host_aware or edge[3]["connection"].ec_numbers)),
                     key=lambda edge: (edge[0], _edge_order((edge[2], edge[3]))))
        for node in allowed
    }
    distances = {node: 0 for node in sources}
    todo = deque(sorted(sources))
    while todo:
        node = todo.popleft()
        for following in eligible.successors(node):
            if following in allowed and following not in distances:
                distances[following] = distances[node] + 1
                todo.append(following)
    if target not in distances:
        return [], True

    def ec_lower_bound(node, selected, downstream):
        if not host_aware:
            return 0
        # 放宽为只要求宿主/已补 EC 的可达性，不检查未决依赖。若仍不可达，至少需一个新 EC。
        # 不在此剪掉可能由更上游副产物消解的依赖，以免高估。
        pending, reached = [node], {node}
        while pending:
            current = pending.pop()
            if current in sources:
                return 0
            for source, _, _, data in incoming[current]:
                connection = data["connection"]
                if source in downstream or source in reached:
                    continue
                if connection.is_endogenous or connection.ec_numbers & selected:
                    reached.add(source)
                    pending.append(source)
        return 1

    serial = count()
    queue, seen = [], set()
    routes, emitted = [], set()

    def push(state, steps):
        node, selected, downstream, _, _ = state
        extra = ec_lower_bound(node, selected, downstream)
        # 未知新 EC 使用空字符串作字典序下界；全局步数下界不考虑循环约束。
        heappush(queue, (len(selected) + extra, steps + distances[node],
                        ("",) * extra + tuple(sorted(selected)), next(serial), steps, state))

    initial = (target, frozenset(), frozenset({target}), frozenset(), ())
    seen.add(initial)
    push(initial, 0)
    settled = rejected = 0
    last_progress = perf_counter()
    while queue:
        _, _, _, _, steps, state = heappop(queue)
        node, selected, downstream, obligations, path = state
        settled += 1
        if settled == 1 or settled % 10000 == 0:
            logger.debug("%s已处理 %d 个状态；当前 %d 种补充 EC、%d 步，待处理 %d 个状态；已排除 %d 个循环候选",
                        "最少补酶搜索" if host_aware else "最少步数搜索", settled, len(selected), steps, len(queue), rejected)
            now = perf_counter()
            if now - last_progress >= 10:
                logger.info("搜索进行中：已找到 %d/%d 条候选", len(routes), top_n)
                last_progress = now
        if node in sources and not obligations:
            if path not in emitted:
                emitted.add(path)
                connections = [graph[a][b][key]["connection"] for a, b, key in path]
                trace = [graph[a][b][key]["event"] for a, b, key in path]
                name = f"TCA → {target}" + ("（最少补酶）" if host_aware else "")
                routes.append((PathWay(name, connections), trace, set(selected)))
                logger.debug("已找到候选 %d/%d：%s → %s；需补 %d 种 EC（%s），共 %d 步",
                            len(routes), top_n, node, target, len(selected),
                            "、".join(sorted(selected)) or "无需补充", steps)
                if len(routes) == top_n:
                    # 尚未展开本终点的上游，不能因队列暂空就声称全部穷尽。
                    return routes, False
        if node in sources and obligations:
            rejected += 1
        for source, _, edge_key, data in incoming[node]:
            if source not in distances or source in downstream:
                continue
            connection = data["connection"]
            substrates, products = set(connection.substrate_ids), set(connection.product_ids)
            unresolved = (obligations - (products - substrates)) | ((substrates - {source}) & downstream)
            unresolved = unresolved - {source}
            options = (None,) if not host_aware or connection.is_endogenous else sorted(connection.ec_numbers)
            for ec in options:
                updated = selected if ec is None else selected | {ec}
                successor = (source, updated, downstream | {source}, frozenset(unresolved),
                             ((source, node, edge_key),) + path)
                if successor in seen:
                    continue
                seen.add(successor)
                push(successor, steps + 1)
    logger.debug("路线依赖检查结束：处理 %d 个状态，排除 %d 个循环候选，穷尽后共 %d 条候选", settled, rejected, len(routes))
    return routes, True


class PathwayEnumeration(Iterator[CandidatePath]):
    def __init__(self, network: PathwayNetwork):
        self.network = network
        self.exhausted = False
        self.count = 0
        self._iterator = self._walk()

    @property
    def is_complete(self) -> bool:
        return self.exhausted and self.network.is_complete

    def __next__(self) -> CandidatePath:
        try:
            candidate = next(self._iterator)
        except StopIteration:
            if not self.exhausted:
                logger.info("候选枚举结束：共 %d 条，已穷尽当前网络的化合物简单路径；数据%s",
                            self.count, "完整" if self.network.is_complete else "不完整，不能确认全局最优")
            self.exhausted = True
            raise
        self.count += 1
        if self.count == 1 or self.count % 1000 == 0:
            logger.info("已生成 %d 条候选路线，仍在枚举中", self.count)
        return candidate

    def _walk(self):
        graph = self.network.graph
        logger.info("开始枚举 %s 的候选路线：%d 个 TCA 接点，默认不设条数上限",
                    self.network.target_id, len(self.network.terminal_nodes))
        allowed = set(self.network.terminal_nodes)
        queue = deque(sorted(allowed))
        while queue:
            for successor in graph.successors(queue.popleft()):
                if successor not in allowed:
                    allowed.add(successor)
                    queue.append(successor)

        def incoming(node):
            return iter(sorted((edge for edge in graph.in_edges(node, keys=True, data=True)
                                if edge[0] in allowed), key=lambda edge: (edge[0], repr(edge[2]))))

        for start in self.network.initial_nodes:
            if start not in allowed:
                continue
            visited, edges = {start}, []
            stack = [(start, incoming(start))]
            if start in self.network.terminal_nodes:
                yield CandidatePath(PathWay(f"TCA → {self.network.target_id}", []), [])
            while stack:
                node, iterator = stack[-1]
                edge = next(iterator, None)
                if edge is None:
                    stack.pop()
                    visited.remove(node)
                    if edges:
                        edges.pop()
                    continue
                source = edge[0]
                if source in visited:
                    continue
                visited.add(source)
                edges.append(edge)
                stack.append((source, incoming(source)))
                # 一个接点可能位于另一个接点的下游；两者起始的简单路径都属于候选。
                if source in self.network.terminal_nodes:
                    ordered = list(reversed(edges))
                    pathway = PathWay(f"TCA → {self.network.target_id}",
                                      [item[3]["connection"] for item in ordered])
                    if not analyze_pathway_dependencies(pathway)["circular_dependencies"]:
                        yield CandidatePath(pathway, [item[3]["event"] for item in ordered])


def _reaction_key(connection):
    return (tuple(sorted(set(connection.reaction_ids))), connection.direction,
            tuple(sorted(connection.substrate_ids)), tuple(sorted(connection.product_ids)))


def _origin_order(origin):
    return (origin["map_id"], repr(origin["edge_key"]), tuple(sorted(origin["diagram_entry_ids"])))


class _NetworkBuilder(_MapLoader):
    def __init__(self, target, cache, refresh):
        super().__init__(target, cache, refresh)
        self.network = nx.MultiDiGraph(target_id=target, representation="compound", schema_version=2,
                                       map_metadata={}, terminal_origins={})
        self.assembled = set()
        self.terminals = set()
        self.processed_links = set()
        self.invalid_links = set()

    def assemble(self, map_id, compound):
        graph = self.load(map_id, compound)
        if graph is None or map_id in self.assembled:
            return graph
        try:
            records = []
            for cid, data in graph.nodes(data=True):
                entity = data.get("entity")
                if not isinstance(entity, Entity) or entity.id != cid:
                    raise ValueError(f"{map_id} 的实体与节点 {cid} 不一致")
            for source, target, raw_key, data in graph.edges(keys=True, data=True):
                connection = data["connection"]
                if not isinstance(connection, Connection) or (connection.source.id, connection.target.id) != (source, target):
                    raise ValueError(f"{map_id} 的连接对象与网络边 {source} → {target} 不一致")
                origin = {"map_id": map_id, "edge_key": raw_key,
                          "diagram_entry_ids": set(data["diagram_entry_ids"]),
                          "ec_numbers": set(connection.ec_numbers), "ko_ids": set(data["ko_ids"])}
                records.append((source, target, _reaction_key(connection), connection, origin))
        except Exception as error:
            self.failed_maps[map_id] = ("graph", str(error))
            self.fail(map_id, compound, "graph", error)
            self.graphs.pop(map_id, None)
            return None
        # 所有原图记录先校验；仅对新容器求并集，不改原图 Connection/注释集合。
        before_nodes, before_edges = len(self.network), self.network.number_of_edges()
        for cid, data in graph.nodes(data=True):
            entity = self.entities.setdefault(cid, data["entity"])
            if cid not in self.network:
                self.network.add_node(cid, entity=entity, compound_id=cid, map_ids=set())
            self.network.nodes[cid]["map_ids"].add(map_id)
        for source, target, key, connection, origin in records:
            if not self.network.has_edge(source, target, key):
                merged = Connection(self.entities[source], self.entities[target], reaction_ids=list(key[0]),
                                    ec_numbers=set(), direction=connection.direction,
                                    substrate_ids=list(key[2]), product_ids=list(key[3]))
                self.network.add_edge(source, target, key=key, connection=merged, ko_ids=set(),
                                      origins=[], map_ids=set(), kind="reaction", weight=1)
            data = self.network[source][target][key]
            data["connection"].ec_numbers.update(connection.ec_numbers)
            data["ko_ids"].update(origin["ko_ids"])
            data["map_ids"].add(map_id)
            data["origins"].append(origin)
            data["origins"].sort(key=_origin_order)
            first = data["origins"][0]
            data["map_id"] = first["map_id"]  # 仅兼容代表性来源，不参与寻路。
            data["diagram_entry_ids"] = set(first["diagram_entry_ids"])
            data["event"] = ReactionStep(first["map_id"], source, target, first["edge_key"],
                                          tuple(sorted(first["diagram_entry_ids"])), network_edge_key=key)
        self.network.graph["map_metadata"][map_id] = graph.graph
        self.assembled.add(map_id)
        new_nodes = len(self.network) - before_nodes
        new_edges = self.network.number_of_edges() - before_edges
        logger.debug("已接入 %s：新增 %d 个化合物、%d 条反应，合并 %d 个共享化合物、%d 条重复反应；"
                    "当前共 %d 个化合物、%d 条反应", map_id, new_nodes, new_edges,
                    len(graph) - new_nodes, graph.number_of_edges() - new_edges,
                    len(self.network), self.network.number_of_edges())
        return graph

    def run(self) -> PathwayNetwork:
        try:
            maps = _initial_maps(self.target, self.cache, self.refresh)
        except Exception as error:
            self.fail(None, self.target, "association", error)
            maps = []
        for map_id in maps:
            if self.assemble(map_id, self.target) is not None and map_id == _TCA:
                logger.debug("目标 %s 本身位于 TCA 图 map00020，允许零反应路线", self.target)
                self.terminals.add(self.target)
                self.network.graph["terminal_origins"][self.target] = [{"map_id": _TCA, "kind": "initial_map"}]
        round_number = 0
        while self.target in self.network:
            # 同一化合物在其他已加载图中补出的反应也计入全局上游，循环直到无新图。
            reachable = nx.ancestors(self.network, self.target) | {self.target}
            before = len(self.assembled)
            round_number += 1
            logger.info("组网进度：第 %d 轮，已接入 %d 张 map，上游化合物 %d 个",
                        round_number, before, len(reachable))
            for map_id in sorted(self.assembled):
                if map_id == _TCA or (map_id, None) in self.invalid_links:
                    continue
                try:
                    links = self.graphs[map_id].graph["map_links"]
                    contacts = sorted(reachable & links.keys())
                except Exception as error:
                    self.fail(map_id, self.target, "graph", error)
                    self.invalid_links.add((map_id, None))
                    continue
                for contact in contacts:
                    if (map_id, contact) in self.invalid_links:
                        continue
                    try:
                        destinations = links[contact]
                        if not isinstance(destinations, (list, tuple, set)) or any(
                            not isinstance(mid, str) or _MAP_ID.fullmatch(mid) is None for mid in destinations
                        ):
                            raise ValueError(f"{contact} 的跨图链接格式错误")
                    except Exception as error:
                        self.fail(map_id, contact, "graph", error)
                        self.invalid_links.add((map_id, contact))
                        continue
                    for destination in sorted(set(destinations)):
                        link = (map_id, contact, destination)
                        if link in self.processed_links:
                            continue
                        self.processed_links.add(link)
                        if destination == _TCA:
                            if contact not in self.terminals:
                                logger.debug("找到 TCA 接点：%s 中的 %s 明确链接 map00020", map_id, contact)
                            self.terminals.add(contact)
                            self.network.graph["terminal_origins"].setdefault(contact, []).append(
                                {"map_id": map_id, "linked_map": _TCA, "kind": "map_link"})
                        else:
                            if destination not in self.assembled and destination not in self.failed_maps and destination not in self.skipped_maps:
                                logger.debug("发现上游图：%s 的 %s 链接到 %s", map_id, contact, destination)
                            self.assemble(destination, contact)
            if len(self.assembled) == before:
                logger.debug("上游扩展结束：没有新的可达 map 需要接入")
                break
        initial = (self.target,) if self.target in self.network else ()
        terminal = frozenset(self.terminals)
        self.network.graph.update(initial_nodes=initial, terminal_nodes=terminal)
        originals = dict(sorted(self.graphs.items()))
        return PathwayNetwork(self.target, self.network, originals, initial, terminal,
                              not self.failures, self.failures, self.skipped_maps)


def build_pathway_network(target_id: str, *, cache_dir: str | Path | None = None,
                          refresh: bool = False) -> PathwayNetwork:
    """按 C 编号合并整张 map，反复发现全局上游链接直到组图闭包。

    map_links 只决定继续获取哪些图以及明确的 TCA 接点，不限制已加载反应的衔接。
    同一化合物只有一个节点，等价反应合并 EC/KO 并保留所有 origins；不添加换图边。
    graphs 保留原图，原有数据获取/404/重试规则不变；不分析 GBFF、不设隐式截断。
    """
    if not isinstance(target_id, str) or re.fullmatch(r"C[0-9]{5}", target_id.strip().upper()) is None:
        raise ValueError("化合物编号应为 C 加五位数字，例如 C05432")
    cache = Path(cache_dir) if cache_dir is not None else kegg.DEFAULT_CACHE_DIR
    target = target_id.strip().upper()
    started = perf_counter()
    logger.info("开始构建 %s 的合并代谢网络", target)
    logger.debug("缓存目录：%s；%s", cache,
                "强制刷新数据" if refresh else "优先复用缓存")
    result = _NetworkBuilder(target, cache, refresh).run()
    logger.info("组网完成：接入 %d 张 map、%d 个化合物、%d 条反应，找到 %d 个 TCA 接点；"
                "404 跳过 %d 张，异常记录 %d 条；数据%s；用时 %.2f 秒",
                len(result.graph.graph["map_metadata"]), len(result.graph), result.graph.number_of_edges(),
                len(result.terminal_nodes), len(result.skipped_maps), len(result.failures),
                "完整" if result.is_complete else "不完整", perf_counter() - started)
    if result.failures:
        logger.warning("有 %d 条异常记录，结果仅适用于当前网络；详情见报告或使用 --verbose 查看",
                       len(result.failures))
    return result
