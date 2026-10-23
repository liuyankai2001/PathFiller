
# PathFiller

基于 **KEGG 通路图与宿主基因组注释** 的代谢通路补全工具。

输入目标化合物的 KEGG 编号（`Cxxxxx`）和宿主的 GBFF 文件，PathFiller 会自动获取关联通路、合并代谢网络，并搜索从 TCA 接点到目标的候选路线，给出需要补充的酶的 **EC 编号**。

默认返回 **3 条候选路线**，优先选择需要补充的不同 EC 数量最少的路线，再比较反应步数。结果包含中文报告、结构化 JSON 和 SVG 路线图。

## 功能

- **自动组网**：通过目标编号查询初始 map，沿可达的跨图链接扩展，联合 RN、EC、KO 三种 KGML 构建网络。
- **化合物合并**：相同 `Cxxxxx` 在 NetworkX `MultiDiGraph` 中对应同一个节点，保留平行反应和各 map 的来源信息。
- **宿主匹配**：从 GBFF 的非伪基因 CDS 中提取、去重 EC 注释，与网络中的反应匹配。
- **前 N 条补酶路线**：默认 3 条，可自定义数量；同一种补充 EC 在整条路线中只计一次。
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

以番茄红素 `C05432` 为目标，默认输出 3 条候选：

```bash
uv run --locked python main.py C05432 --gbff data/genomic.gbff
```

以丁酸 `C00246` 为目标，输出最多 5 条候选：

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
| `--top-n` | `3` | 请求的候选数量，必须为正整数；不足时返回实际找到的数量 |
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
| `report.md` | 中文候选对比表、逐步反应、宿主匹配、需补 EC 和未核查的原料供应 |
| `result.json` | 批量结果版本 1，包含候选排名、连接、EC 覆盖、来源、完整性及异常信息 |
| `network.json` | 完整网络版本 2，包含合并节点、反应边及原图数据与注释 |
| `network.svg` | 可选的完整网络大图，包含全部合并节点与反应，按候选路线着色 |
| `pathway_*.svg` | 每条候选的独立路线图，可用浏览器打开 |

独立候选 SVG 中，每个化合物对应一个圆，圆内显示 C 编号，圆下显示完整 KEGG 名称，长名称自动换行；箭头标注 EC。**绿色**表示宿主已有 EC 注释匹配，**红色**表示未匹配，深蓝色圆表示该路线使用的 TCA 接点。图中保留可检索文字和节点、边的追溯元数据，布局由 Graphviz 自动计算。

### 完整网络大图

`--export-network-svg` 绘制所有已加载 map 合并后的完整网络；`--top-n 5` 返回五条时，大图会突出全部五条。前三条依次为蓝、橙、紫，其他候选自动分配颜色，图例标明候选编号。共享反应使用多色并列线，非候选反应使用深灰细线；节点、业务反应边的数量和编号保持不变。

大图只常驻显示圆内 C 编号。用浏览器直接打开 `network.svg`，悬停节点查看化合物名称，悬停边查看 EC、R 编号、候选归属和宿主匹配状态。图中使用颜色区分候选，独立候选图继续保留绿／红的匹配含义。

大图使用固定种子的 `sfdp` 布局，高亮线会绕开无关化合物圆点，独立候选图使用 `dot`。完整名称、EC、KO、原图来源和候选归属保存在 SVG 元数据中。启用大图会为整个网络查询名称，首次运行通常比只输出候选图更慢，后续复用缓存；默认关闭时仍只查询候选路线的化合物名称。

`result.json` 的 `network_svg_file` 保存成功生成的大图文件名，`route_colors` 保存候选编号到颜色的映射。默认关闭时分别为 `null` 和空对象。大图导出失败时保留独立候选图及 JSON、报告，错误记入 `export_errors`，命令返回非零退出码；搜索完整性不受绘图错误影响。

没有候选时仍输出报告和网络；零反应路线会绘制为单节点。名称查询失败会保留已有名称或仅显示编号，记录在 `name_lookup_warnings` 中。SVG 导出失败时保留其他结果，记录 `export_errors`，命令返回非零退出码。

### 如何判断结果完整性

| 字段 | 含义 |
| --- | --- |
| `data_is_complete` | 当前搜索范围内的数据获取和解析是否完整，已跳过的 404 图排除在范围外 |
| `search_complete` | 当前网络内前 N 条候选的判定是否完成，或搜索队列是否已经耗尽 |
| `is_complete` | 数据和前 N 条搜索是否同时完整 |
| `exhausted` | 是否已经穷尽当前网络内符合规则的候选路线 |
| `failures` | 数据获取、解析或链接校验的失败记录 |
| `skipped_maps` | 因 HTTP 404 跳过的 map 及原因 |

取得前 N 条不代表所有路线已经穷尽。存在数据异常时，返回的是**当前网络内的候选**，不能据此确认缺失分支中不存在更优路线。

## 工作原理与适用边界

```text
宿主 GBFF → 提取 EC 集合
                         ↘
目标 C 编号 → KEGG 关联 map → 合并网络 → 宿主匹配 → 前 N 条搜索 → 名称补齐 → 结果导出
```

**组网与起点。** 仅在开始查询目标关联 map，排除总图、概览图、结构图和非普通代谢图；后续通过图内 `map_links` 发现新图。RN 提供反应参与物和图示方向，EC、KO 提供注释。搜索起点由图中明确指向 TCA 图 `map00020` 的接点确定；目标本身位于成功加载的 TCA 图时，也支持零反应路线。固定化合物名单和 `Entity.is_tca` 不作为搜索终止条件。

**宿主匹配。** 反应的任一候选 EC 命中 `Genome.ec_numbers`，即标记 `Connection.is_endogenous=True`；否则为 `False`。多个 EC 视为可选关系，部分 EC 编号按字符串精确匹配。

**路线排名。** 使用反向 A* 搜索，依次比较不同补充 EC 数、反应步数和补充 EC 编号字典序。同一种 EC 用于多步时只计一次，同一路线的多个 EC 选择只保留一组最优覆盖；不同反应边组成的路线分别计数。无 EC 注释的反应保留在网络中，但不参与最少补酶搜索。

**当前模型。** 每条候选是一条按箭头连接的化合物主链。程序保留反应中列出的多个底物和产物，并排除路线内部的循环依赖；其他前体的供应分支会列为未核查，尚未联合搜索。补充 EC 数不等同于所需基因数，当前结果也未验证计量、碳流、表达和酶的实际底物适用性。

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
| `src.host_match.compare_pathway()` | 比较一条路线与宿主，精确计算最少补充 EC 覆盖 |
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
