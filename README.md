
# PathFiller

基于 **KEGG 通路图与宿主基因组注释** 的代谢通路补全工具。

输入目标化合物的 KEGG 编号（`Cxxxxx`）和宿主的 GBFF 文件，PathFiller 会自动获取关联通路、合并代谢网络，并搜索从 TCA 接点到目标的候选路线，给出需要补充的酶的 **EC 编号**。

默认返回 **3 套不同的补酶方案**，按需补 EC 集合区分，每套展示一条最短合法代表路线。优先比较补酶数量，再比较反应步数。结果包含中文报告、结构化 JSON 和 SVG 路线图。

## 功能

- **自动组网**：通过目标编号查询初始 map，沿可达的跨图链接扩展，联合 RN、EC、KO 三种 KGML 构建网络。
- **化合物合并**：相同 `Cxxxxx` 在 NetworkX `MultiDiGraph` 中对应同一个节点，保留平行反应和各 map 的来源信息。
- **宿主匹配**：从 GBFF 的非伪基因 CDS 中提取、去重 EC 注释，与网络中的反应匹配。
- **前 N 套补酶方案**：默认 3 套，可自定义数量；相同需补 EC 集合只计一套，同一种 EC 用于多步只计一次。
- **结果可视化**：每条候选独立生成 SVG，显示化合物编号、完整名称、反应 EC 和宿主匹配颜色。
- **可选完整大图**：额外导出全部已加载网络，用不同颜色突出各候选，悬停查看名称和反应注释。
- **缓存与容错**：缓存 KEGG 数据，支持刷新、失败重试、404 跳过和数据完整性记录。

## 快速开始

### 1. 准备环境

| 组件 | 要求 |
| --- | --- |
| Python | 3.12 或更高；项目默认使用 3.12 |
| uv | 用于安装依赖和运行项目，见 [安装说明](https://docs.astral.sh/uv/getting-started/installation/) |
| Graphviz | SVG 绘图需要独立安装，并将 `dot` 加入 `PATH`；完整大图还需要 `sfdp`，见 [下载页面](https://graphviz.org/download/) |
| 网络 | 首次获取尚未缓存的 KEGG 数据时需要访问 `rest.kegg.jp` |

克隆仓库并进入项目根目录后，安装锁定的 Python 依赖：

```bash
uv sync --locked
```

主要 Python 依赖为 Biopython、NetworkX 和 Matplotlib，完整配置见 [pyproject.toml](pyproject.toml)。依赖锁定与同步方式见 [uv 文档](https://docs.astral.sh/uv/concepts/projects/sync/)。

确认 Graphviz 可用：

```bash
dot -V
```

### 2. 准备宿主基因组

在项目根目录创建 `data/`，将宿主的 GenBank 格式基因组文件放入其中，例如：

```text
data/
└── genomic.gbff
```

文件需要包含 CDS 的 `EC_number` 注释。程序会读取全部序列记录，排除标记为伪基因的 CDS，并合并为宿主的 EC 集合；仅有基因组序列的 FASTA 文件不能用于当前的 EC 匹配。

宿主基因组不随仓库提供，也可通过 `--gbff` 指定自己的文件路径。

### 3. 运行

以番茄红素 `C05432` 为目标，默认输出 3 套不同补酶方案：

```bash
uv run --locked python main.py C05432 --gbff data/genomic.gbff
```

以丁酸 `C00246` 为目标，输出最多 5 套不同补酶方案：

```bash
uv run --locked python main.py C00246 --gbff data/genomic.gbff --top-n 5
```

省略目标与宿主参数时，使用 `C05432` 和项目根目录的 `data/genomic.gbff`：

```bash
uv run --locked python main.py
```

默认日志显示阶段进度、组网汇总、搜索结果和重要异常。排查问题时开启详细日志：

```bash
uv run --locked python main.py C05432 --verbose
```

同时导出完整网络 SVG，并用不同颜色突出返回的全部候选：

```bash
uv run --locked python main.py C05432 --export-network-svg
```

## 命令行参数

```text
main.py [target] [--gbff FILE] [--top-n N] [--cache-dir DIR]
        [--refresh] [--verbose] [--export-network-svg] [--output-dir DIR]
```

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `target` | `C05432` | 目标化合物编号，格式为 `C` 加五位数字 |
| `--gbff` | `data/genomic.gbff` | 宿主 GBFF 文件；默认路径基于项目根目录 |
| `--top-n` | `3` | 请求的不同补酶方案数量，必须为正整数；不足时返回实际数量，相同 EC 集合不重复计数 |
| `--cache-dir` | 项目根目录 `cache/` | KEGG 缓存目录，可显式指定其他位置 |
| `--refresh` | 关闭 | 忽略已有缓存，重新请求本次所需的数据 |
| `--verbose` | 关闭 | 显示请求、重试、缓存和逐图处理等调试日志 |
| `--export-network-svg` | 关闭 | 额外生成完整的 `network.svg`，按候选分别着色，名称和 EC 仅在悬停时显示 |
| `--output-dir` | `outputs/<目标>/<唯一运行编号>/` | 指定输出目录，必须位于项目根目录的 `outputs/` 内 |

例如，将结果写入指定目录：

```bash
uv run --locked python main.py C05432 --output-dir outputs/lycopene_run
```

相对的 `--output-dir` 按项目根目录解释。指定目录已有同名业务结果时会报错，避免覆盖历史文件。

## 输出结果

一次运行的输出结构：

```text
outputs/
└── C05432/
    └── <运行编号>/
        ├── report.md
        ├── result.json
        ├── network.json
        ├── network.svg          # 仅启用 --export-network-svg 时生成
        ├── pathway_001.svg
        ├── pathway_002.svg
        └── pathway_003.svg
```

| 文件 | 内容 |
| --- | --- |
| `report.md` | 中文方案对比表、代表路线逐步反应、宿主匹配、需补 EC 和未核查的原料供应 |
| `result.json` | 批量结果版本 1，`candidate_unit="extra_ec_set"`；包含方案排名、代表路线、选定 EC、来源、完整性及异常 |
| `network.json` | 完整网络版本 2，包含合并节点、反应边及原图数据与注释 |
| `network.svg` | 可选的完整网络大图，包含全部合并节点与反应，按方案代表路线着色 |
| `pathway_*.svg` | 每套方案的独立代表路线图，可用浏览器打开 |

每套方案的独立 SVG 显示代表路线及选定的需补 EC；每个化合物对应一个圆，圆内显示 C 编号，圆下显示完整 KEGG 名称，长名称自动换行；箭头仍标注全部候选 EC。**绿色**表示宿主已有 EC 注释匹配，**红色**表示未匹配，深蓝色圆表示该路线使用的 TCA 接点。图中保留可检索文字和节点、边的追溯元数据，布局由 Graphviz 自动计算。

### 完整网络大图

`--export-network-svg` 绘制所有已加载 map 合并后的完整网络；`--top-n 5` 返回五套时，大图会突出全部五套方案的代表路线。前三套依次为蓝、橙、紫，其他方案自动分配颜色，图例标明方案编号，悬停图例查看需补 EC。同一路径可代表不同 EC 方案，共享反应使用多色并列线，非候选反应使用深灰细线；节点、业务反应边的数量和编号保持不变。

大图只常驻显示圆内 C 编号。用浏览器直接打开 `network.svg`，悬停节点查看化合物名称，悬停边查看 EC、R 编号、候选归属和宿主匹配状态。图中使用颜色区分候选，独立候选图继续保留绿／红的匹配含义。

大图使用固定种子的 `sfdp` 布局，高亮线会绕开无关化合物圆点，独立候选图使用 `dot`。完整名称、EC、KO、原图来源和候选归属保存在 SVG 元数据中。启用大图会为整个网络查询名称，首次运行通常比只输出候选图更慢，后续复用缓存；默认关闭时仍只查询候选路线的化合物名称。

`result.json` 的 `network_svg_file` 保存成功生成的大图文件名，`route_colors` 保存方案编号到颜色的映射。默认关闭时分别为 `null` 和空对象。大图导出失败时保留独立方案图及 JSON、报告，错误记入 `export_errors`，命令返回非零退出码；搜索完整性不受绘图错误影响。

没有候选时仍输出报告和网络；零反应路线会绘制为单节点。名称查询失败会保留已有名称或仅显示编号，记录在 `name_lookup_warnings` 中。SVG 导出失败时保留其他结果，记录 `export_errors`，命令返回非零退出码。

### 如何判断结果完整性

| 字段 | 含义 |
| --- | --- |
| `data_is_complete` | 当前搜索范围内的数据获取和解析是否完整，已跳过的 404 图排除在范围外 |
| `search_complete` | 当前网络内前 N 套方案的判定是否完成，或搜索队列是否已经耗尽 |
| `is_complete` | 数据和前 N 套搜索是否同时完整 |
| `exhausted` | 是否已经穷尽当前网络内符合规则的不同补酶方案 |
| `failures` | 数据获取、解析或链接校验的失败记录 |
| `skipped_maps` | 因 HTTP 404 跳过的 map 及原因 |

取得前 N 套不代表全部方案已经穷尽。存在数据异常时，返回的是**当前网络内的候选方案**，不能据此确认缺失分支中不存在更优方案。

## 工作原理与适用边界

```text
宿主 GBFF → 提取 EC 集合
                         ↘
目标 C 编号 → KEGG 关联 map → 合并网络 → 宿主匹配 → 前 N 套方案搜索 → 名称补齐 → 结果导出
```

**组网与起点。** 仅在开始查询目标关联 map，排除总图、概览图、结构图和非普通代谢图；后续通过图内 `map_links` 发现新图。RN 提供反应参与物和图示方向，EC、KO 提供注释。搜索起点由图中明确指向 TCA 图 `map00020` 的接点确定；目标本身位于成功加载的 TCA 图时，也支持零反应路线。固定化合物名单和 `Entity.is_tca` 不作为搜索终止条件。

**宿主匹配。** 反应的任一候选 EC 命中 `Genome.ec_numbers`，即标记 `Connection.is_endogenous=True`；否则为 `False`。多个 EC 视为可选关系，部分 EC 编号按字符串精确匹配。

**方案排名。** 使用反向 A* 搜索，依次比较补充 EC 数、代表路线的反应步数和 EC 编号字典序。相同 EC 集合只输出一次，保留最短合法代表路线；完全同分时按固定遍历顺序选择。同一种 EC 用于多步只计一次。同一路线保留全部并列最少的 EC 组合，例如 A＋B、A＋C 都最少时可分别输出；如果该路线补 A 已经足够，补 B＋C 的更多酶选择不会输出。不同路线可以有不同的最少补酶数，统一参与排名。无 EC 注释的反应保留在网络中，但不参与补酶搜索。

每套方案的代表路径可能相同，区别在于选用的 EC。空 EC 集合也只计一套；存在零步路线时优先以它作为空集合的代表。去重只影响方案输出，原网络的方向、平行反应和 `network.iter_paths()` 的全量路径枚举保持原样。寻找不同方案可能比取得同数量的反应路线需要探索更多状态。

**当前模型。** 每套方案以一条按箭头连接的化合物主链为代表。程序保留反应中列出的多个底物和产物，并排除路线内部的循环依赖；其他前体的供应分支会列为未核查，尚未联合搜索。补充 EC 数不等同于所需基因数，当前结果也未验证计量、碳流、表达和酶的实际底物适用性。

搜索范围由可达的详细 map 及其显式链接确定，不覆盖 KEGG 中的所有反应；默认没有隐式图数或步数上限，较大的网络和较多的候选可能需要更长时间。

## Python API

在项目根目录使用：

```python
from src.genome_parser import parse_genome
from src.search import build_pathway_network, find_minimum_ec_paths
from src.result_export import export_results

genome = parse_genome("data/genomic.gbff")
network = build_pathway_network("C05432")
result = find_minimum_ec_paths(network, genome, top_n=3)

for rank, candidate in enumerate(result.candidates, start=1):
    print(rank, sorted(candidate.extra_ecs), len(candidate.pathway.connections))

export_results(network, genome, result, "outputs/my_api_run")
# 可选：额外输出完整大图
# export_results(network, genome, result, "outputs/my_network_run", include_network_svg=True)
```

`network.graph` 为按化合物合并的 `MultiDiGraph`，`network.graphs` 保存各来源 map。组网与导出接口均支持 `cache_dir`、`refresh`；自定义缓存时应将相同参数传递给两者。

其他常用接口：

| 接口 | 用途 |
| --- | --- |
| `src.kegg.parse_map()` | 联合解析单张 map 的 RN、EC、KO 数据 |
| `src.kegg.get_tca_graph()` | 获取 TCA 图 `map00020` |
| `src.kegg.get_compound_names()` | 批量获取并缓存化合物名称 |
| `src.host_match.match_network()` | 使用 Genome 更新网络的宿主匹配状态 |
| `src.host_match.compare_pathway()` | 比较一条路线与宿主，返回一组按编号选定的最少 EC 覆盖；方案搜索会保留其他并列最优组合 |
| `src.search.find_minimum_ec_path()` | 返回单条最少补酶候选；无路线时 `pathway=None` |
| `network.iter_paths()` | 逐条枚举候选，用于独立分析；默认不截断 |

Python API 的详细日志可通过 `src.logging_config.configure_logging("DEBUG")` 开启。

## 数据与目录

```text
PathFiller/
├── main.py                  # 命令行入口
├── pyproject.toml           # 项目与依赖配置
├── uv.lock                  # 依赖锁文件
├── src/
│   ├── schemas/             # Entity、Connection、PathWay、Genome
│   ├── genome_parser.py     # GBFF 与宿主 EC 提取
│   ├── kegg.py              # KEGG 请求、缓存和重试
│   ├── map_parser.py        # KGML 联合解析
│   ├── map_supplements/     # 已核对的图示反应补充记录
│   ├── search.py            # 合并网络与路线搜索
│   ├── host_match.py        # 宿主匹配和精确 EC 覆盖
│   ├── result_export.py     # 统一结果导出
│   └── svg_export.py        # SVG 绘图
├── data/                    # 自行准备的宿主文件
├── cache/                   # 运行时生成的 KEGG 缓存
└── outputs/                 # 运行结果
```

`data/`、`cache/`、`outputs/` 已被 Git 忽略，不随仓库分发。KEGG 缓存包括 KGML、反应方程、关联 map、分类信息和 `Cxxxxx.name.txt` 名称；已有有效缓存可复用，缺失的数据仍需要联网获取。

请求失败的处理规则：HTTP 404 跳过对应 map；连接错误、超时及 HTTP 408、429、5xx 最多尝试 5 次（包含首次），每次超时 30 秒，使用带随机扰动的指数退避。重试耗尽或数据校验失败会记录异常并继续其他分支。

KEGG 数据及 API 的使用范围与访问限制请遵循 [KEGG 官方说明](https://www.kegg.jp/kegg/rest/)。图数据格式参考 [KGML 文档](https://www.kegg.jp/kegg/xml/docs/)。
