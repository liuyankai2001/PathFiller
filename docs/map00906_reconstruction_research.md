# map00906 信息重构调研

调研日期：2026-10-22。目标是为 PathFiller 重建可搜索、可关联 EC、方向符合 KEGG 图示的网络，并为后续绘图保留布局信息。

## 结论

推荐以 **RN KGML 为反应骨架、原图 CONF 为节点与标注清单、R 条目为参与物核对来源**，对缺失连接单独核对原图。当前 map00906 的数据足以构建覆盖原图全部 R 编号的网络；其中 4 种反应仍需使用人工核对的补充记录。

此前选择 EC KGML 是为了直接取得 EC，但它不适合作为完整原图的唯一骨架。此次对照发现，EC 版未声明原图中的另外 45 种反应。这与“已声明反应中缺失 4 种完整 reaction 记录”是两个不同问题。

“重构网络”与“逐像素复刻原图”需要分别验收：前者在此次离线试验中已实现 R 编号全覆盖，后者的全图连线形状和全部箭头尚未完成验证。

## 1. 对实际文件的核对

输入是上次核查时直接从 KEGG REST 下载的原始文件。下载 URL、UTC 时间和 SHA-256 保存在 `tests/output/map00906_audit/report.json`，此次没有把人工补充记录算入原始 KGML 的统计。

原图 [map00906 CONF](https://rest.kegg.jp/get/map00906/conf) 中有 **128 个不同 C 编号、149 个不同 R 编号**。这些是编号数，不是圆点或箭头数量：例如一个 Phytoene 圆点同时关联两个 C 编号。

| 数据源 | XML reaction 元素数 | 完整 reaction 元素覆盖的不同 R 编号 | entry 中声明的不同 R 编号 | 相对原图 R 清单缺少完整记录 |
| --- | ---: | ---: | ---: | ---: |
| [ec00906 KGML](https://rest.kegg.jp/get/ec00906/kgml) | 115 | 100 | 104 | 49 |
| [ko00906 KGML](https://rest.kegg.jp/get/ko00906/kgml) | 138 | 117 | 121 | 32 |
| [rn00906 KGML](https://rest.kegg.jp/get/rn00906/kgml) | 166 | 145 | 149 | 4 |

统计方法是从 CONF 提取去重后的 `Rxxxxx` 集合，与各 KGML 的 `reaction.name`、`entry.reaction` 集合分别比较；三种 KGML 都包含原图清单中的 128 个 C 编号。完整差集见 `tests/output/map00906_research/source_coverage.json`。

这一差异与 KEGG 对各视图的定义一致：EC、KO、Reaction 视图分别突出原图中对应类型的条目；不能假定三种视图包含相同的完整反应清单。[KEGG 查看器说明](https://www.kegg.jp/kegg/document/help_pathway.html)

RN 缺失完整记录的 4 种反应是：

| R 编号 | 图上主连接（经原图核对的方向） | R 条目用于核对的内容 |
| --- | --- | --- |
| [R04218](https://www.kegg.jp/entry/R04218) | C03427 → C05421 | Prephytoene-PP 与 15-cis-Phytoene |
| [R07270](https://www.kegg.jp/entry/R07270) | C03427 → C05413 | Prephytoene-PP 与 all-trans-Phytoene |
| [R09653](https://www.kegg.jp/entry/R09653) | C05421 → C19765 | 15-cis-Phytoene 与 15,9'-dicis-Phytofluene |
| [R09692](https://www.kegg.jp/entry/R09692) | C05421 → C05414 | 15-cis-Phytoene 与 all-trans-Phytofluene |

这 4 种反应对应 6 个图形节点，记录已保存在 `src/map_supplements/map00906_reactions.xml`。单向方向来自先前提供的原图截图，R 方程只用于区分参与物；方程的 `<=>` 不能代替图示方向。

## 2. 已做的离线可行性验证

在内存中将 6 条已人工核对的补充记录加入 RN 原始 XML，再使用现有解析器解析。没有修改原始文件或当前生产入口。

结果见 `tests/output/map00906_research/rn_topology_probe.json`：

- 128 个化合物 Entity。
- 149 个不同 R 编号，与原图 CONF 的 R 清单数量一致。
- 150 条有向 Connection；一个图示双向反应生成两个方向的边。
- 没有孤立化合物节点。
- C00353 能到达 C05432；番茄红素分支的无前驱边界为 C00353。
- 保留 C00353、C00448 与 map00900 的跨图链接。

这个试验依赖现有人工补充记录，不是“算法自动从图片恢复了缺失反应”的结果。试验暂未给 RN 节点补 EC，因此不能直接把这个试验网络用于宿主 EC 匹配。

此外，对 CONF 中的 27 条显式反应折线，按端点坐标匹配 RN 的底物、产物图形节点，并核对箭头多边形落在哪个端点。27 条均与 RN 图示方向一致，包括双向情况。结果见 `tests/output/map00906_research/conf_direction_crosscheck.json`。这只验证了 CONF 明确保存的这部分连线，不代表整张 PNG 的所有箭头都已验证。

## 3. 现有工具与研究能解决什么

### KGML 解析和修补工具

[KEGGgraph](https://bioconductor.org/packages/release/bioc/html/KEGGgraph.html) 将 KGML 转为图模型。[KNeXT](https://github.com/everest-castaneda/knext) 使用 Python 和 NetworkX，可输出基因网络、混合网络及节点坐标，但其输入仍然是 KGML。这些工具可供实现参考，无法仅凭更换解析器取得原文件未记录的箭头。

[KEGGParser 原始论文](https://academic.oup.com/bioinformatics/article/29/4/518/200267) 明确讨论 KGML 与静态图片不一致的问题，支持部分自动修复及人工增删节点、边和反转方向。它提供了可借鉴的半自动工作流，但论文也明确指出，KGML 内部信息不总足以支持全部自动修复。该工具主要基于 MATLAB/Scilab，不需要为本项目专门迁移运行环境。

### 直接从图片提取图结构

[Amdex：Automated meta-data extraction from KEGG pathways](https://www.cs.uit.no/hdl/papers/bachelor-kenneth.pdf)，Kenneth Knudsen，2014，直接研究用 OpenCV 从 KEGG PNG 中提取元素。细读其摘要、评估及结论后，不能将其视作完整网络提取器：论文主要验证基因框识别，报告平均约 80% 的识别结果；结论明确说明尚不能识别完整的基因交互连接。附录说明代码随 CD-ROM 提供，本次没有找到可核验的公开 Amdex 代码仓库。它证明问题可研究，但不提供已验证的 map00906 完整还原方案。

[pathCLIP 论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11363067/) 与[作者代码](https://github.com/yangyang-69/pathCLIP) 展示了从通路图识别基因及有向激活、抑制关系的方法。其发表版本在检测到的基因关系图像片段上，有方向关系提取的平均 F1 为 68.97%。这是基因关系任务指标，不能当作代谢通路或 map00906 的准确率。代码还依赖检测/OCR 模块与模型权重。对本项目，适合作为箭头检测方法参考，需要额外验证化合物、酶框和反应分支的适配效果。

直接图像处理可以采用二值化、轮廓检测、线条骨架化与端点追踪。对应基础组件是 [OpenCV 轮廓检测](https://docs.opencv.org/4.x/d4/d73/tutorial_py_contours_begin.html) 和 [scikit-image skeletonize](https://scikit-image.org/docs/stable/auto_examples/edges/plot_skeleton.html)。这些函数不直接解决图形语义：交叉线是否相连、多个箭头是否共用一段线、一个圆点对应哪个异构体，都需要额外约束。

因此，本次未找到并验证“输入任意 KEGG PNG 即可完整、准确输出 map00906 有向网络”的现成工具。这里的结论限于本次检索和验证范围。

## 4. 适合 PathFiller 的实现路线

1. **原图清单**：解析 CONF/网页热点，读取 C、R、EC、通路链接及图形坐标，用它检查重构是否漏项。CONF 中未列出的信息不能默认不存在于整个代谢数据库。
2. **反应骨架**：主要读取 `rn00906/kgml`，按其中 substrate、product 和 reaction.type 创建连接，沿用已有 MultiDiGraph。
3. **编号与异构体**：用 R 条目核对一个圆点含多个 C 编号时的参与物，逐个 R 建边。R 方程不决定图示单、双向。
4. **EC 标注**：从 CONF 同一酶框的 EC/R 共现信息及 R 条目关联补 EC。无可靠 EC 的反应保留未确定状态，不制造 EC 编号。
5. **缺失处理**：对本图的 4 种缺失反应使用已核对的补充记录；其他图先检测差集，再做局部图像识别或人工核对。文件缺记录不能作为生物学上的搜索停止条件。
6. **布局保存**：在独立布局数据中保存“图形节点 ID → 坐标及 C/R 编号列表”。Entity 仍按 C 编号唯一，绘图层保留一个圆点可能对应多个 Entity 的关系，无需给 Entity 重新添加 appearances。

对于本图，仅需在缺失连接附近研究图片箭头识别，已有 CONF 坐标可用于定位，节点编号不必依靠 OCR。优先验证这 4 种连接，再决定是否扩大图像识别范围。

若绘图只要求 KEGG 风格，保留节点坐标并重新布线即可；若要求所有曲线、分叉位置与原图完全一致，还需要提取未在 CONF 中保存的线条路径，或用原图作背景提供交互覆盖层。后一种展示方式保留原图外观，但本身不证明背后的网络已正确重构。

## 5. 验收边界

“覆盖 149 个 R 编号”只代表目录覆盖，不代表 149 种反应的每一条连接和方向均已逐一与 PNG 核验。接入前应检查：

- 128 个 C 编号、149 个 R 编号与原图清单逐项对应。
- 图示反应与 Connection 的对应关系明确，重复酶框不重复计数。
- 对一个圆点的多个化合物，不把异构体随意混连。
- 缺失反应、未知 EC、未知方向分别保留可检查的状态。
- C00353 到 C05432 的两条主要分支，以及到 map00900 的跨图链接保持连通。
- 共底物可供应性另行评估；图示主连接连通不等于宿主可以完成全部反应。

实现更新：生产 `parse_map` 已改用 RN 提供连接、EC/KO 提供注释的联合流程，并保留原图节点布局和跨图链接。补充记录按 map 校验，只填补缺失的图内节点/R 对。固定快照回归覆盖 map00906 的 128 个化合物、149 个不同 R 编号、150 条有向连接，以及 TCA 的 20 个化合物、64 条有向连接。上文离线试验文件保留其历史含义，不代表自动恢复了缺失箭头或逐像素复刻原图。
