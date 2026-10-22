# PathFiller

基于 KEGG 通路图和宿主 GBFF 注释查找候选合成路线、筛选需要补充的酶。

安装项目依赖：

```powershell
uv sync
```

程序默认实时输出中文控制台日志，包含基因组解析、关联 map 查询、缓存读取、下载与重试等待、筛图、
上游扩展、图合并数量、TCA 接点、宿主匹配、补酶搜索状态数、候选枚举进度和搜索结果。HTTP 404 跳过与异常失败分别记录。
时间格式为 `16:30:00 [信息] 已接入 map00900：……`，默认写到标准错误流。

可在调用前调整日志级别，配置只作用于 PathFiller，重复调用不会重复输出：

```python
from src.logging_config import configure_logging

configure_logging("DEBUG")    # 显示反应方程缓存等详细信息
# configure_logging("WARNING")  # 只查看警告和错误
```

在项目的 Python 环境中使用：

先解析宿主 GBFF，生成保存宿主名称、来源路径和 EC 集合的 `Genome`：

```python
from src.genome_parser import parse_genome

genome = parse_genome("data/genomic.gbff")
print(genome.name)
print(genome.file_path)
print(len(genome.ec_numbers))
```

解析器使用 Biopython 读取文件中的全部染色体/质粒记录，提取并去重非伪基因 CDS 的
`EC_number`。CDS 自身或对应 gene 标有 `pseudo`／`pseudogene` 时排除；
名称来自文件的 organism 注释，来源路径保存为绝对路径，解析过程输出中文日志。
有效 GBFF 没有 EC 注释时返回空集合；文件没有可解析记录或格式错误时明确报错。
`Genome` 只保存 `file_path`、`name`、`ec_numbers`；解析一次即可用于多个目标的网络匹配。

用宿主寻找需要补充的不同 EC 数量最少的路线：

```python
from src.genome_parser import parse_genome
from src.host_match import match_network, compare_pathway
from src.search import build_pathway_network, find_minimum_ec_path

genome = parse_genome("data/genomic.gbff")
network = build_pathway_network("C05432")
result = find_minimum_ec_path(network, genome)  # 内部自动完成网络匹配

if result.pathway is not None:
    print("需要补充的 EC：", sorted(result.extra_ecs))
    print("反应步数：", len(result.pathway.connections))
    print("数据是否完整：", result.is_complete)
    comparison = compare_pathway(result.pathway, genome)
    print(comparison["summary"])
```

`match_network(network.graph, genome)` 原地重新设置所有合并 `Connection.is_endogenous`，返回同一网络。
反应任一 EC 命中 `genome.ec_numbers` 即为 `True`，其余为 `False`；不修改保留的原图及来源记录。
更换宿主时重新匹配。`compare_pathway(pathway, genome)` 只读比较，不依赖连接上次匹配的宿主，
返回逐步结果和一组精确的最少 EC 覆盖 `one_minimum_extra_ec_set_under_or_assumption`。
没有候选路线时不能调用比较；零步 `PathWay` 的补充集合为空。

`find_minimum_ec_path()` 从全部明确的 TCA 接点沿箭头搜索目标，优先最少不同补充 EC 数，其次最少反应步数。
同一步多个 EC 为 OR 选项，同一 EC 用于多步只计一次。搜索状态保存化合物、已补 EC 集合、
下游路线节点及待上游消解的依赖，不能只按当前化合物合并状态。
反向 A* 使用乐观的步数和补酶数下界直接寻找候选，不预先枚举全部路线，不设隐式数量上限；
同成本按 EC 编号、节点及边稳定选择。
普通 NetworkX 边权累加不能直接表达“同一 EC 只算一次”，因此该接口在 NetworkX 网络上管理上述搜索状态。
状态数可能随不同 EC 组合快速增长，中文日志会显示处理进度。

三个寻路入口（最少补酶、最少反应步数、全部路线枚举）都排除内部循环依赖：
如果某步需要的化合物只能在本路线后面得到，这个起点的候选被排除，继续比较其他路线。
起点和前面反应的全部图示产物视为已产生；由早期副产物满足的需求不会误判为循环。
其他未知原料仅列入 `route_review.unverified_substrates`，不推断其供应。
可用 `src.pathway_validation.analyze_pathway_dependencies(pathway)` 单独检查保存的路线；
返回 `no_internal_cycle_detected` 也不表示计量、碳流和全部原料供应已验证。

`MinimumEcResult` 保留 `SearchResult` 的 `pathway`、`trace`、`graphs`、`is_complete`、`failures`、`skipped_maps`，
新增 `extra_ecs: set[str]`。无路线以 `pathway is None` 判断；空 EC 集合也可能表示有效的零补充路线。
无 EC 注释的边保留在网络中用于绘图，但不参与最少补酶搜索；因此“没有路线”针对有注释的搜索范围。
`is_complete=False` 时只能称当前最佳候选，不能确认范围内全局最优。尚未验证共底物、表达和酶的底物适用性。

直接搜索的可运行示例：

```powershell
.venv/Scripts/python.exe tests/inspect_minimum_ec_path.py --fixture-snapshot
# 使用在线数据和 cache 缓存，不加 --fixture-snapshot
.venv/Scripts/python.exe tests/inspect_minimum_ec_path.py C05432
```

报告写入独立的 `tests/output/minimum_ec_C05432/`；`result.json` 保留路线、EC、来源和完整性，
`result.md` 便于阅读。三图快照与真实宿主 GBFF 的当前结果为 16 步、3 种 EC：
`1.3.99.26`、`2.5.1.29`、`2.5.1.32`。选择来自实际搜索，未硬编码；快照缺少其他分支，数据不完整。

查询与解析 KEGG 通路：

```python
from src.kegg import get_maps_by_compound, parse_map

maps = get_maps_by_compound("C05432")
graph = parse_map("map00906")

print(graph.number_of_nodes())
print(graph.number_of_edges())
print(list(graph.predecessors("C05432")))

entity = graph.nodes["C05432"]["entity"]
for source, target, key, data in graph.in_edges("C05432", keys=True, data=True):
    connection = data["connection"]
    print(source, target, connection.reaction_ids, connection.ec_numbers)
```

`parse_map` 联合读取同一通路的 RN、EC、KO 三种 KGML：
RN 提供化合物连接与图示方向，EC 提供酶编号，KO 提供功能组编号。
返回 NetworkX `MultiDiGraph`：节点以 C 编号标识并保存 `Entity`，
边保存 `Connection`。图示单向反应只建立一条边，双向反应建立两个方向的边。
方向始终来自 RN KGML，不被其他视图或 R 条目的方程方向覆盖。
注释按图内节点编号和 R 编号核对，等价的重复反应记录合并关联 EC、KO 和原图节点编号。
缺少 EC 的连接仍保留，其 `ec_numbers` 为空，`is_endogenous=False`。
多底物反应的边表示图上的连接，完整底物、产物列表保存在对应 `Connection` 中。

部分 KGML 的反应节点列有 R 编号，却没有完整的 `<reaction>` 记录。
`map00906` 中缺失的 6 条记录已对照原图补充到
`src/map_supplements/map00906_reactions.xml`，补充箭头方向来自原图，R 方程只用于区分
同一圆点代表的不同异构体。补充数据是本图的人工核对记录，不能替代所有通路的自动解析。
补充数据按 map 编号保存，解析器检查其图内节点、R 编号、参与物是否与 RN 一致。
同一节点仅部分反应被上游补齐时，只补仍缺失的 R，不覆盖已有反应。
其他缺失记录未获补齐时明确报错，
避免将漏解析造成的断点当作通路起点。

多个酶框共用同一条连接时，补充文件可用 `reaction_alias` 引用已有的 RN 反应，
不要求每个酶框都有独立的 `<reaction>`。例如已对照原图核对的
`src/map_supplements/map00650_reactions.xml` 保存 `205 → 71`：

```xml
<reaction_alias id="205" source_id="71" name="rn:R03026" type="reversible">
  <substrate id="125" name="cpd:C00877" />
  <product id="137" name="cpd:C01144" />
</reaction_alias>
```

记录中的参与物和方向是核对条件。解析器检查通路编号、两端节点声明的 R 编号、
完整原始反应的参与物和方向，然后复用原始连接，合并 EC、KO 和 `diagram_entry_ids`。
`map00650` 因此解析为 47 个化合物、111 条有向边；`R03026` 的两个方向均保存
EC `4.2.1.17`、`4.2.1.150` 以及图内节点 `71`、`205`。
复用关系必须明确记录，不按相同 R 编号或坐标接近自动推断，也不依赖复用记录之间的链式引用。
记录不符合当前数据时明确报错；KEGG 已提供完整反应的节点／R 组合始终使用原始数据。
共用关系通过中文日志输出，并保存在 `graph.graph["reaction_aliases"]` 中。

重新测试丁酸并保留之前的失败报告：

```powershell
.venv/Scripts/python.exe tests/inspect_minimum_ec_path.py C00246 --output-dir tests/output/minimum_ec_C00246_after_alias
```

本次丁酸实测已接入 54 张 map，得到 2188 个化合物、4242 条反应。
当前图搜索候选为 `C00024 → C00010 → C00136 → C00246`，3 步、需补 EC `2.8.3.8`，
但前两步分别使用同一反应 `R01177` 的正向和反向，第一步还需要 `C00136` 作共底物，存在循环依赖。
该结果暴露了当前多底物反应拆成两两连接的限制，不能作为有效丁酸合成方案。
共用连接解析已修复；本次未调整反应建边和搜索规则。其他 map 数据失败使 `is_complete=False`。
本次结果的独立核查及说明保存在新报告目录的 `verification.json`，旧失败报告保留。

加入内部循环排除后，旧的 `R01177 → R01177 → R01179` 候选不再接受。
重新搜索得到 `R03778 → R01177 → R01179`，化合物链仍为
`C00024 → C00010 → C00136 → C00246`，3 步、1 种补充 EC `2.8.3.8`。
该候选仍需图外原料 `C01944`、`C05269`，本轮没有核查它们的供应，不能直接据此认定只需补一个酶。
结果保存在独立目录 `tests/output/minimum_ec_C00246_cycle_checked`，此前的失败及循环候选报告保留。

图中的跨图链接也会保留，例如：

```python
graph.graph["map_links"].get("C00353", [])  # ["map00900"]
```

这里的 `map00900` 对应原图的“Terpenoid backbone biosynthesis”圆角框。
链接表示可前往另一张图继续查找，不作为化合物之间的反应边。

先组装一个大的 NetworkX 网络，再枚举所有无环候选路线：

```python
from src.search import build_pathway_network

result = build_pathway_network("C05432")
graph = result.graph  # nx.MultiDiGraph，反应边按合成方向保存
routes = result.iter_paths()
for candidate in routes:
    print(candidate.pathway.name, len(candidate.pathway.connections))
    for connection in candidate.pathway.connections:
        print(connection.source.id, "->", connection.target.id, sorted(connection.ec_numbers))
    # trace 保留代表性原图来源和指向合并边的 network_edge_key
print("是否穷尽：", routes.exhausted, "整体是否完整：", routes.is_complete)
```

所有已接入图按 C 编号合并为一个 `MultiDiGraph`。例如 `graph.nodes["C00022"]` 只有一个节点，
糖酵解与 MEP 的反应可在这里直接衔接，引用同一个共享 `Entity`。
网络仅有 `kind="reaction", weight=1` 的反应边，没有换图边或虚拟 TCA 节点。
`initial_nodes` 为目标编号的元组，`terminal_nodes` 为明确 TCA 接点的 C 编号集合。

只有起终点、R 编号集合、图示方向、完整底物和产物均一致时才合并反应。
合并边使用新的 `Connection`，汇总 EC 和 `ko_ids`；原图对象不变。
每条边的 `origins` 保存所有来源 map、原图边编号及图内反应节点编号，避免不同图的节点编号混淆。
`event` 保留固定排序选出的代表性来源；`event.network_edge_key` 可直接定位合并边。
原图保存在 `result.graphs`，坐标和图形元数据保存在 `graph.graph["map_metadata"]`。

map 只作为下载和来源追溯单位。开始时查询目标关联图，后续通过已加载图的 `map_links` 发现新图。
每张图只解析一次；每次接入新图后重新计算整个网络中目标的上游节点，
并重新检查所有已加载图的链接，直到没有新的可达图。不在第一次找到 TCA 时停止组图。
明确指向 `map00020` 的接点直接成为终点，不为该链接额外下载 TCA 网络。

`result.iter_paths()` 在统一网络中反向枚举，保留 MVA、MEP 和平行反应。
每条候选不重复化合物且没有内部循环依赖，不同合并边算不同候选；
即使路径经过另一个 TCA 接点，也保留各起点满足依赖条件的候选。
输出逐条生成，默认没有路线数或反应步数上限。
`result.is_complete` 表示组图数据完整性；`routes.exhausted` 表示枚举正常结束，
`routes.is_complete` 要求二者同时满足。提前退出循环不能称为已穷尽全部候选。
404 跳过、重试、筛图、共享实体及缓存规则与下方的单路线接口一致。

正式匹配逻辑位于 `src.host_match`。如需查看全部候选及排名，可运行枚举示例；GBFF 只解析一次，后续传递 `Genome`：

```powershell
uv run --no-project --python .venv/Scripts/python.exe tests/inspect_pathway_network.py
# 固定三图快照，离线复现；未提供的分支记为失败，整体不会误报为完整
uv run --no-project --python .venv/Scripts/python.exe tests/inspect_pathway_network.py --fixture-snapshot
# 仅显式预览前 1000 条，报告标记未穷尽；默认无此上限
uv run --no-project --python .venv/Scripts/python.exe tests/inspect_pathway_network.py --fixture-snapshot --max-paths 1000
```

报告默认输出到 `tests/output/genome_pathway_network_C05432`，旧 `merged_pathway_network_C05432` 结果保留：
`network.json` 使用 `schema_version=2, representation="compound"`，直接保存完整的
`graph`、`nodes[].attributes`、`edges[].attributes` 和原图信息，绘图无需再从单图恢复边属性；
`pathways.jsonl` 每行保存一条路线；`ranking.tsv` 和 `ranking.sqlite` 按需补 EC 数、反应步数排序；
`summary.json` 单独记录数据完整性、枚举完整性、失败和 404 跳过；`best_candidate.json` 保存当前首位候选。
GBFF 提取沿用非伪基因 CDS 的精确 EC 匹配，同一步候选 EC 暂按 OR 关系计算编号覆盖，
重复 EC 只算一次；无 EC 注释的步骤标为未知。结果不是已验证的基因引入方案，仍未检查共底物供应。
包含无 EC 注释步骤的路线无法计算补充数，排在可计算路线之后。
新报告仅保留宿主名称、来源文件、EC 匹配与排名，不再导出 `ec_index`、`host_evidence` 或基因、蛋白、CDS 明细。
单独比较已保存路线可运行 `tests/inspect_host_pathway_match.py --search-report <报告路径>`，
输出到 `tests/output/genome_host_pathway_<目标>/`；旧 GBFF 检查脚本已删除，历史明细文件仍保留。
穷尽枚举的候选数量可能很大，完整性只对上述简单路径定义和已获取数据范围成立。

当前三图固定快照合并为 **201 个化合物、311 条反应边**。
原始简单边路径共 **14,736 条**，排除内部循环依赖后保留 **7,056 条候选**；历史报告保留原数量。
按宿主需要补充的不同 EC 数、反应步数排序，当前首位为 16 步、3 种 EC 的编号覆盖候选。
MEP 与 CrtI 分支保留，参考集合 `2.5.1.29`、`2.5.1.32`、`1.3.99.31` 也能覆盖对应路线；
排序未硬编码这组答案。同长同覆盖数的集合按编号排序，因此首位集合可能选择其他 CrtI 候选 EC。
快照未提供的其他分支记为失败，`enumeration_exhausted=True`、`data_is_complete=False`，
不能将这个排名称为所有 KEGG 通路中的已确认最优。

联合网络可导出为保留数据的 SVG（需要系统 PATH 中可用的 Graphviz `dot`）：

```python
from src.svg_export import export_network_svg

export_network_svg(result.graph, "figures/merged_pathway_network_C05432/merged_pathway_network_C05432.svg")
```

同一 C 编号只画一个圆，中心标 `Cxxxxx`，不按 map 分区。
反应边带箭头并标出全部 `EC x.x.x.x`，无注释标 `EC unknown`；标签始终绘制在路径上方。
布局由 Graphviz 计算，SVG 保留真实文字与独立图形对象，不复刻 KEGG 原图布局。
节点、边有稳定独立 ID、悬停详情和 `data-*` 属性；`metadata#pathfiller-network` 包含 JSON，
保留图元数据、全部节点和边属性，以及 Entity、Connection、EC、KO、全部 origins 和原图坐标信息。
导出不修改输入网络；原有单 map 网络也可使用这个接口。
新图和实际 SVG 渲染预览保存在 `figures/merged_pathway_network_C05432/`，可用其中的
`render_network.py` 从新版 `network.json` 重现。
`tests/audit_merged_network.py` 可离线核对固定快照报告：将保存的全部候选与 NetworkX
`all_simple_edge_paths()` 逐条比较，并检查排名、三 EC 参考覆盖和 SVG 中的完整属性。

原来的单条最短候选接口保持兼容，继续按反应步数选择路线，不考虑宿主 EC：

```python
from src.search import search_shortest_path

result = search_shortest_path("C05432")
print("是否完成最短判定：", result.is_complete)
if result.pathway is not None:
    print("反应步数：", len(result.pathway.connections))
    for connection in result.pathway.connections:
        print(connection.source.id, "->", connection.target.id,
              connection.reaction_ids, sorted(connection.ec_numbers))
else:
    print("范围内无路线" if result.is_complete else "数据不完整，无法确定")

for failure in result.failures:
    print(failure.map_id, failure.compound_id, failure.stage, failure.reason)
for map_id, reason in result.skipped_maps.items():
    print("404 跳过：", map_id, reason)
```

搜索先查询目标的关联 map，排除总图、概览图和结构图，再按 KEGG `CLASS` 选择
普通详细代谢图。地图类型参照 [KEGG PATHWAY 分类](https://www.kegg.jp/kegg/pathway.html)，
`map01040`、`map01056`、`map01057`、`map01059` 等可作 KO 展开的图保留，
不能把所有 `010` 编号一律视为结构图。多张合格初始图都会参加搜索。
后续只沿 `map_links` 获取新图，不再逐一查询中间化合物的关联图。

`search_shortest_path()` 与枚举接口共同使用上述化合物合并网络。
完成组图闭包后，在 NetworkX 网络上用反向 A* 比较候选，并记录下游节点及依赖。
广度优先距离只作乐观步数下界；存在未消解依赖的 TCA 接点不能提前结束搜索。
每条反应成本为 1，相同化合物自然连接，无需换图成本。
单向、双向及平行反应保留原网络方向；同长路线按固定顺序选一条。
不要求接点的入度为 0；图链接循环不会造成重复加载，反应环不影响最短距离计算。

只有目标本身位于成功加载的初始 `map00020`，或当前路线能沿箭头到达一个明确链接 `map00020` 的接点，
才算到达 TCA。例如在 `map00900` 找到 C00024 时仍需进入 `map00010`，
确认 C00024 在该图的链接指向 `map00020`。比较所有已发现 TCA 接点的距离后才选择候选。
`Entity.is_tca` 不参与停止判定。

`SearchResult` 的数据用于展示路线和后续绘图：

| 属性 | 含义 |
| --- | --- |
| `pathway` | 一条候选 `PathWay`，连接按“上游 → 目标”排列；未找到为 `None` |
| `is_complete` | 是否完成当前范围内最短判定；404 图排除在范围外，其他已探索分支失败会使其为 `False` |
| `failures` | 失败的 map、化合物、阶段和原因；阶段包括 association、classification、download、parse、compound、graph |
| `trace` | 按合成方向保存 `ReactionStep`，含代表性原图来源和 `network_edge_key`；不再产生 `MapSwitch` |
| `graphs` | 成功加载的各图，完整保留 EC、KO、坐标和原图节点信息 |
| `skipped_maps` | 因加载资源返回 HTTP 404 而排除的 map 及原因，默认空字典 |

加载某张图所需的分类、RN／EC／KO 视图或 R 条目返回 **HTTP 404** 时，停止扩展整图，
原因放入 `skipped_maps`，不放入 `failures`，也不会仅因此使 `is_complete=False`。
最短判定范围排除这些跳过图。缺少本地缓存本身不代表 404。
初始目标关联查询失败时没有具体 map 可以跳过，仍记入 `failures`。

网络连接失败、超时、HTTP 408／429／5xx 会在统一请求层自动重试，最多 **5 次尝试，
包含首次请求**，每次请求超时仍为 30 秒。重试采用 Full Jitter 指数退避：

```python
# r = 1、2、3、4，等待上限依次为 1、2、4、8 秒
delay = random.uniform(0, min(30, 1 * 2 ** (r - 1)))
```

每次尝试都遵守现有至少 0.4 秒的请求间隔；退避等待在限速锁之外执行。
404 和其他不可重试 HTTP 错误直接返回；解析错误及错误视图不进行网络重试。
请求最终失败抛出 `src.kegg.KeggRequestError`，它兼容 `RuntimeError`，携带
`endpoint`、`status_code`（网络错误为 `None`）和 `attempts`，同时保留原始异常原因。
直接调用 `parse_map()` 时仍会收到这个异常；搜索层负责将对应图跳过或记录失败。

重试耗尽、其他下载错误、分类或解析异常会进入 `failures`，继续搜索其他分支。
同一次搜索内，已跳过或最终失败的图不再启动新的请求轮次；下一次独立搜索可以重新请求。
有候选但 `is_complete=False` 时，只能称为“当前最佳候选”，不能确认最短。
`pathway=None, is_complete=True` 表示范围内没有路线；
`pathway=None, is_complete=False` 表示数据不完整，尚无法确定。
零反应路线用非空 `PathWay` 对象及空 `connections` 表示。

反应追踪记录可以定位代表性原图边，也可以在 `PathwayNetwork.graph` 中定位合并边的全部来源：

```python
network = build_pathway_network("C05432")
candidate = next(network.iter_paths(), None)
if candidate is not None:
    for step in candidate.trace:
        edge = network.graph[step.source_id][step.target_id][step.network_edge_key]
        print(edge["origins"], edge["ko_ids"])
```

一次搜索内每张图最多解析一次，跨图共享同一个 `Entity` 字典。
`cache_dir` 默认为项目根目录的 `cache/`，不受启动工作目录影响。
所有 API 和检查脚本共同使用 `src.kegg.DEFAULT_CACHE_DIR`；显式传入 `cache_dir` 或
`--cache-dir` 时使用指定目录。`Cxxxxx.maps.json` 缓存初始关联列表、
`mapxxxxx.pathway.txt` 缓存分类条目；`refresh=True` 会重新获取这些数据及加载图的
RN、EC、KO 和所需 R 方程。所需数据均已缓存时支持离线复用。
不设置隐式步数或图数上限。最短判定只覆盖排除 404 图后获取的详细代谢图的合并网络，
结果是图示候选，不表示共底物、宿主酶或 GBFF 分析已满足。

运行搜索示例并导出核查报告：

```powershell
uv run python tests/inspect_search.py C05432
# 可选：--cache-dir cache --refresh --output-dir tests/output/merged_search
```

报告保存为 `tests/output/merged_search/search_C05432.json` 和 Markdown，包含完整追踪记录、失败分支、
404 跳过图及成功加载图的注释和布局。固定三图快照仍得到从 C00024 出发的 12 步最短候选，
每步原图来源保留，不产生换图事件。未提供快照的其他分支记为失败，不能误报已确认最短。
404 与网络重试通过离线回归单独验证。本次仅提供 Python API 和示例，主程序入口不变。

三种原始 KGML 和消除多 C 编号歧义所需的 R 方程缓存到 `cache/`；
添加 `refresh=True` 可以重新获取三种视图及本次需要的 R 方程。
三种视图和所需方程均已缓存时可以离线使用。RN 获取失败或视图不匹配时明确报错。
离线快照缓存单独保存在 `cache/snapshots/<报告目录名>/`，不与在线缓存混用；
指定 `--cache-dir` 时，快照写入该目录的 `snapshots/<报告目录名>/`。
宿主 GBFF 仍放在 `data/`，报告与 SVG 保持原输出位置。历史输出中的旧路径仅作来源记录。
化合物名称默认留空，内源状态默认 `False`；TCA 标记使用
`src/map_parser.py` 中约定的名单，包含入口乙酰辅酶 A。

跨图解析时，可传入共享的实体字典，让同一 C 编号复用同一个对象：

```python
entities = {}
graph = parse_map("map00906", entities=entities)
```

TCA 循环（`map00020`）可通过专用入口获取，返回相同类型的网络：

```python
from src.kegg import get_tca_graph

tca_graph = get_tca_graph()
print(tca_graph.number_of_nodes())
print(tca_graph.number_of_edges())

# 与其他通路复用实体时，传入同一个字典。
entities = {}
tca_graph = get_tca_graph(entities=entities)
target_graph = parse_map("map00906", entities=entities)
```

图保留 KGML 中显示的全部化合物，包括部分辅因子。
`Entity.is_tca` 是预览标记；跨图搜索通过明确的 `map00020` 链接判断是否到达 TCA。
双向反应在网络中有两个方向的边，边数不等于反应数。

显示 TCA 网络并保存 PNG、SVG：

```powershell
uv run python -m src.visualization
```

窗口支持缩放、平移和保存。添加 `--no-show` 只导出图片，默认保存在
`figures/tca_network/`；使用 `--map map00906` 可查看其他通路。
这是自动布局的网络预览，节点位置与 KEGG 原图不同。

在 Python 中显示已有网络：

```python
import matplotlib.pyplot as plt
from src.visualization import draw_network

figure = draw_network(tca_graph)
plt.show()
```

蓝色节点表示项目约定的 TCA 核心中间体及入口，供预览使用。连接颜色读取
`Connection.is_endogenous`：绿色表示内源，红色表示未匹配到宿主 EC。
新连接默认 `False`，调用 `match_network` 后颜色反映指定宿主的匹配结果；
绘图导入旧报告时仍兼容历史 `None` 状态的深灰色显示。

已有 KGML 数据可以直接解析，避免联网。第一个参数必须是 RN 视图，
`ec_xml`、`ko_xml` 为可选的注释视图；省略时只构图，不填相应注释。
图节点含多个 C 编号时，需提供返回
对应 R 方程字符串的 `reaction_lookup` 回调来确认参与物，不能任意推断连接：

```python
from src.map_parser import parse_kgml

graph = parse_kgml(
    rn_xml_text,
    ec_xml=ec_xml_text,
    ko_xml=ko_xml_text,
    reaction_lookup=lookup_equation,
)
```

KO 和原图布局信息保存在网络属性中，三个数据对象的字段不变：

```python
for source, target, data in graph.edges(data=True):
    print(data["connection"].ec_numbers)
    print(data["ko_ids"])             # 关联的 Kxxxxx 集合
    print(data["diagram_entry_ids"])  # 对应的原图反应节点编号集合

layout = graph.graph["diagram_entries"]
print(layout["327"]["compound_ids"])  # ["C05421", "C05413"]
print(layout["327"]["graphics"])      # KGML 中的坐标、形状、尺寸等原始属性
```

`diagram_entries` 按图内节点编号保存，允许多个图内节点对应同一个 `Entity`，
也允许一个圆点关联多个异构体。坐标只用于后续布局，不代表全部原图曲线都已提取。
`source_views` 记录参与构图的视图编号，`supplemented_entry_ids` 记录使用补充数据的节点。
`reaction_aliases` 单独记录实际启用的共用连接，每项保存 `entry_id`、`source_entry_id` 和 `reaction_ids`；
这些酶框也保留在 `diagram_entries` 及相应边的 `diagram_entry_ids` 中。

使用生产入口检查 map00906：

```powershell
uv run python tests/inspect_map00906.py
```

结果保存到 `tests/output/map00906_joint_network.json`、CSV 和 Markdown，原始快照和旧报告保留。
固定快照下，map00906 有 128 个化合物、149 个不同 R 编号、150 条有向连接；
其中 4 种反应依赖已人工核对的补充记录。

运行测试：

```powershell
uv run python -m unittest discover -s tests -p "test_*.py" -v
```
