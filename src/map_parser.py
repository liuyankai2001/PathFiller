"""以 RN KGML 构图，并合并 EC、KO 视图的注释。"""

from collections.abc import Callable
from itertools import product
from pathlib import Path
import re
from typing import Literal
import xml.etree.ElementTree as ET

import networkx as nx

from .schemas import Connection, Entity
from .logging_config import get_logger


logger = get_logger("map_parser")


# 本项目约定的 TCA 核心中间体，以及入口乙酰辅酶 A。
TCA_COMPOUND_IDS = frozenset({
    "C00024", "C00158", "C00417", "C00311", "C00026",
    "C00091", "C00042", "C00122", "C00149", "C00036",
})


def _identifiers(text: str, prefix: str) -> list[str]:
    return list(dict.fromkeys(
        token.removeprefix(prefix) for token in text.split() if token.startswith(prefix)
    ))


def _read_view(xml: str | bytes, view: str, number: str | None = None) -> tuple[ET.Element, dict]:
    try:
        root = ET.fromstring(xml)
    except ET.ParseError as error:
        raise ValueError(f"无法解析 {view.upper()} KGML XML") from error
    actual_number = root.get("number", "")
    if (
        root.tag != "pathway" or re.fullmatch(r"[0-9]{5}", actual_number) is None
        or root.get("name") != f"path:{view}{actual_number}" or root.get("org") != view
        or (number is not None and actual_number != number)
    ):
        raise ValueError(f"KGML 不是预期通路的 {view.upper()} 视图")
    entries = {}
    for entry in root.findall("entry"):
        entry_id = entry.get("id")
        if not entry_id or entry_id in entries:
            raise ValueError(f"图内节点编号缺失或重复：{entry_id}")
        entries[entry_id] = entry
    return root, entries


def _declared_reactions(entry: ET.Element) -> set[str]:
    ids = set(_identifiers(entry.get("reaction", ""), "rn:"))
    if entry.get("type") == "reaction":
        ids.update(_identifiers(entry.get("name", ""), "rn:"))
    if any(re.fullmatch(r"R[0-9]{5}", rid) is None for rid in ids):
        raise ValueError(f"图内节点的 R 编号无效：{entry.get('id')}")
    return ids


def _annotations(xml: str | bytes | None, view: str, root: ET.Element, entries: dict) -> dict:
    """注释按图内节点及 R 编号关联，不将同一 R 的其他图形节点混入。"""
    if xml is None:
        return {}
    _, annotation_entries = _read_view(xml, view, root.get("number"))
    result = {}
    pattern = r"[0-9]+(?:\.(?:[0-9]+|-)){3}" if view == "ec" else r"K[0-9]{5}"

    def position(element):
        return [
            tuple(graphic.get(key) for key in ("x", "y", "type"))
            for graphic in element.findall("graphics")
        ]

    for entry_id, entry in annotation_entries.items():
        identifiers = set(_identifiers(entry.get("name", ""), f"{view}:"))
        if not identifiers:
            continue
        if any(re.fullmatch(pattern, identifier) is None for identifier in identifiers):
            raise ValueError(f"无效的 {view.upper()} 注释：{entry_id}")
        reaction_ids = _declared_reactions(entry)
        original = entries.get(entry_id)
        if original is None or not reaction_ids or not reaction_ids <= _declared_reactions(original):
            raise ValueError(f"{view.upper()} 注释的图内节点或 R 编号与 RN 不一致：{entry_id}")
        # 有坐标时也检查位置，避免把更新后重用的图内编号静默关联到旧快照。
        annotated_position, original_position = position(entry), position(original)
        if annotated_position and original_position and annotated_position != original_position:
            raise ValueError(f"{view.upper()} 注释的图内节点位置与 RN 不一致：{entry_id}")
        result[entry_id] = {rid: set(identifiers) for rid in reaction_ids}
    return result


def _participants(reaction: ET.Element, role: str, entries: dict) -> list[list[str]]:
    groups = []
    for participant in reaction.findall(role):
        entry_id = participant.get("id")
        if entry_id not in entries or entries[entry_id].get("type") != "compound":
            raise ValueError(f"反应引用了不存在的图内节点：{entry_id}")
        name = participant.get("name") or entries[entry_id].get("name", "")
        compounds = _identifiers(name, "cpd:")
        if not compounds or any(re.fullmatch(r"C[0-9]{5}", cid) is None for cid in compounds):
            raise ValueError(f"当前解析器只支持 C 编号的反应参与物：{name}")
        if not set(compounds) <= set(_identifiers(entries[entry_id].get("name", ""), "cpd:")):
            raise ValueError(f"反应参与物与图内节点的 C 编号不一致：{entry_id}")
        groups.append(compounds)
    if not groups:
        raise ValueError(f"反应缺少 {role}：{reaction.get('name')}")
    return groups


def _resolve_participants(
    substrates: list[list[str]], products: list[list[str]], equation: str,
) -> tuple[list[str], list[str]]:
    """用 R 方程消除多 C 编号歧义，仅核对参与物，不采用方程方向。"""
    sides = re.split(r"\s*(?:<=>|=>|<=|->|=)\s*", equation, maxsplit=1)
    if len(sides) != 2:
        raise ValueError("无法从 R 条目读取反应两侧的化合物")
    left, right = (set(re.findall(r"\bC[0-9]{5}\b", side)) for side in sides)
    options = set()
    # 两种对应方式都检查；最终 source/target 方向仍由 KGML 决定。
    for substrate_side, product_side in ((left, right), (right, left)):
        selected_substrates = [[cid for cid in group if cid in substrate_side] for group in substrates]
        selected_products = [[cid for cid in group if cid in product_side] for group in products]
        if all(selected_substrates) and all(selected_products):
            options.add((
                tuple(dict.fromkeys(cid for group in selected_substrates for cid in group)),
                tuple(dict.fromkeys(cid for group in selected_products for cid in group)),
            ))
    if len(options) != 1:
        raise ValueError("多 C 编号节点与 R 条目的对应关系仍不明确，未创建猜测连接")
    resolved_substrates, resolved_products = options.pop()
    return list(resolved_substrates), list(resolved_products)


def _add_connection(
    graph: nx.MultiDiGraph,
    source_id: str,
    target_id: str,
    reaction_ids: list[str],
    ec_numbers: set[str],
    direction: Literal["irreversible", "reversible"],
    substrate_ids: list[str],
    product_ids: list[str],
    ko_ids: set[str],
    diagram_entry_id: str,
) -> None:
    key = (
        tuple(sorted(reaction_ids)), direction,
        tuple(sorted(substrate_ids)), tuple(sorted(product_ids)),
    )
    if graph.has_edge(source_id, target_id, key):
        data = graph.edges[source_id, target_id, key]
        data["connection"].ec_numbers.update(ec_numbers)
        data["ko_ids"].update(ko_ids)
        data["diagram_entry_ids"].add(diagram_entry_id)
        return
    connection = Connection(
        source=graph.nodes[source_id]["entity"],
        target=graph.nodes[target_id]["entity"],
        reaction_ids=list(reaction_ids),
        ec_numbers=set(ec_numbers),
        direction=direction,
        substrate_ids=list(substrate_ids),
        product_ids=list(product_ids),
    )
    graph.add_edge(
        source_id, target_id, key=key, connection=connection,
        ko_ids=set(ko_ids), diagram_entry_ids={diagram_entry_id},
    )


def _reaction_layout(reaction: ET.Element, entries: dict) -> tuple:
    """核对图内参与物编号及 C 编号，忽略 XML 子元素顺序。"""
    direction = reaction.get("type")
    if direction not in {"irreversible", "reversible"}:
        raise ValueError("缺少有效的图示方向")
    roles = []
    for role in ("substrate", "product"):
        groups = _participants(reaction, role, entries)
        roles.append(tuple(sorted(
            (item.get("id"), tuple(sorted(group)))
            for item, group in zip(reaction.findall(role), groups, strict=True)
        )))
    return direction, *roles


def _reuse_reactions(supplement: ET.Element, originals: list[ET.Element], entries: dict,
                    missing: dict) -> tuple[list[ET.Element], list[dict]]:
    """只复用经核对的原始 RN 连接；不从 R 编号、坐标或其他复用记录推断连接。"""
    reactions, reused, claimed = [], [], set()
    required = {entry_id: set(ids) for entry_id, ids in missing.items()}
    for alias in supplement.findall("reaction_alias"):
        entry_id, source_id = alias.get("id"), alias.get("source_id")
        if not required.get(entry_id):
            continue  # KEGG 已补齐该节点，采用原始记录，不覆盖或验证旧的备用映射。
        ids = set(_identifiers(alias.get("name", ""), "rn:"))
        if not ids or not ids <= _declared_reactions(entries[entry_id]):
            raise ValueError(f"反应复用记录的 R 编号与目标节点不一致：{entry_id}")
        active = ids & required[entry_id]
        if not active:
            continue
        if source_id == entry_id or source_id not in entries or not active <= _declared_reactions(entries[source_id]):
            raise ValueError(f"反应复用记录的来源节点或 R 编号不一致：{entry_id} → {source_id}")
        if any((entry_id, rid) in claimed for rid in active):
            raise ValueError(f"反应复用记录重复指定同一节点/R：{entry_id}")
        try:
            expected = _reaction_layout(alias, entries)
            candidates = [reaction for reaction in originals
                          if reaction.get("id") == source_id
                          and active & set(_identifiers(reaction.get("name", ""), "rn:"))]
            available = set().union(*(set(_identifiers(item.get("name", ""), "rn:")) for item in candidates))
            if not active <= available:
                raise ValueError("来源没有所需 R 的完整原始 reaction 记录")
            if any(_reaction_layout(item, entries) != expected for item in candidates):
                raise ValueError("来源底物、产物或方向与核对记录不一致")
        except ValueError as error:
            raise ValueError(f"反应复用记录失效：{entry_id} → {source_id}：{error}") from error
        copied = set()
        for original in candidates:
            selected_ids = tuple(sorted(active & set(_identifiers(original.get("name", ""), "rn:"))))
            if selected_ids in copied:
                continue
            copied.add(selected_ids)
            selected = ET.fromstring(ET.tostring(original))
            selected.set("id", entry_id)
            selected.set("name", " ".join(f"rn:{rid}" for rid in selected_ids))
            reactions.append(selected)
        claimed.update((entry_id, rid) for rid in active)
        missing[entry_id].difference_update(active)
        reused.append({"entry_id": entry_id, "source_entry_id": source_id, "reaction_ids": sorted(active)})
    return reactions, reused


def _complete_reactions(root: ET.Element, entries: dict) -> tuple[list[ET.Element], list[str], list[dict]]:
    """补齐缺失记录或复用已核对的原始连接，不从 R 方程猜箭头方向。"""
    reactions = list(root.findall("reaction"))
    originals = list(reactions)
    recorded = {}
    for reaction in reactions:
        entry_id = reaction.get("id")
        ids = set(_identifiers(reaction.get("name", ""), "rn:"))
        if entry_id not in entries or not ids or not ids <= _declared_reactions(entries[entry_id]):
            raise ValueError(f"反应的 R 编号与图内节点不一致：{entry_id}")
        recorded.setdefault(entry_id, set()).update(ids)
    missing = {
        entry_id: absent for entry_id, entry in entries.items()
        if (absent := _declared_reactions(entry) - recorded.get(entry_id, set()))
    }
    if not missing:
        return reactions, [], []
    supplement_file = Path(__file__).with_name("map_supplements") / f"map{root.get('number')}_reactions.xml"
    if not supplement_file.exists():
        raise ValueError(
            f"KGML 缺少完整 reaction 记录：{list(missing)}；"
            "需要与原图核对的补充数据，不能将断点当作通路起点"
        )
    supplement = ET.fromstring(supplement_file.read_bytes())
    if (
        supplement.tag != "pathway" or supplement.get("number") != root.get("number")
        or supplement.get("name") != f"path:map{root.get('number')}"
    ):
        raise ValueError("反应补充数据与当前通路不一致")
    completed = []
    for reaction in supplement.findall("reaction"):
        entry_id = reaction.get("id")
        if entry_id not in missing:
            continue
        expected = _declared_reactions(entries[entry_id])
        actual = set(_identifiers(reaction.get("name", ""), "rn:"))
        if not actual or actual != expected:
            raise ValueError(f"反应补充数据的 R 编号与图内节点不一致：{entry_id}")
        for participant in list(reaction.findall("substrate")) + list(reaction.findall("product")):
            entry = entries.get(participant.get("id"))
            compounds = set(_identifiers(participant.get("name", ""), "cpd:"))
            if entry is None or entry.get("type") != "compound" or not compounds or not compounds <= set(
                _identifiers(entry.get("name", ""), "cpd:")
            ):
                raise ValueError(f"反应补充数据的参与物与图内节点不一致：{entry_id}")
        active = actual & missing[entry_id]
        if not active:
            continue
        # 只补尚未记录的 R；即使同一图形节点的另一个 R 已被 KEGG 修复也不覆盖。
        selected = ET.fromstring(ET.tostring(reaction))
        selected.set("name", " ".join(f"rn:{rid}" for rid in sorted(active)))
        reactions.append(selected)
        missing[entry_id].difference_update(active)
        if entry_id not in completed:
            completed.append(entry_id)
    aliases, reused = _reuse_reactions(supplement, originals, entries, missing)
    reactions.extend(aliases)
    unresolved = {entry_id for entry_id, absent in missing.items() if absent}
    if unresolved:
        raise ValueError(f"KGML 仍缺少完整 reaction 记录：{sorted(unresolved)}；不能判断通路起点")
    for item in reused:
        logger.info("map%s 共用反应连接：酶框 %s 复用酶框 %s 的 %s；保留已有图示方向并合并 EC、KO 注释",
                    root.get("number"), item["entry_id"], item["source_entry_id"], "、".join(item["reaction_ids"]))
    return reactions, completed, reused


def _map_links(root: ET.Element, entries: dict) -> dict[str, list[str]]:
    """保留通过共同化合物连接的 map，供后续跨图搜索使用。"""
    links = {}
    for relation in root.findall("relation"):
        if relation.get("type") != "maplink":
            continue
        for entry_id in (relation.get("entry1"), relation.get("entry2")):
            entry = entries.get(entry_id)
            if entry is None or entry.get("type") != "map":
                continue
            maps = [
                "map" + name[-5:] for name in _identifiers(entry.get("name", ""), "path:")
                if re.fullmatch(r"(?:map|ec|ko|rn)[0-9]{5}", name)
            ]
            for subtype in relation.findall("subtype"):
                if subtype.get("name") != "compound":
                    continue
                compound = entries.get(subtype.get("value"))
                if compound is None or compound.get("type") != "compound":
                    raise ValueError("跨图链接引用了不存在的化合物节点")
                for cid in _identifiers(compound.get("name", ""), "cpd:"):
                    related = links.setdefault(cid, [])
                    related.extend(map_id for map_id in maps if map_id not in related)
    return links


def parse_kgml(
    xml: str | bytes,
    *,
    ec_xml: str | bytes | None = None,
    ko_xml: str | bytes | None = None,
    entities: dict[str, Entity] | None = None,
    reaction_lookup: Callable[[str], str] | None = None,
) -> nx.MultiDiGraph:
    """读取 RN KGML，返回以 C 编号为节点标识的化合物网络。

    节点的 entity 属性保存 Entity，边的 connection 属性保存 Connection。
    可传入已有实体字典，在多张图之间复用同一 Entity；返回图只含本图节点。
    一个反应涉及多个底物/产物时建立图示参与物之间的边，并在 Connection
    中保留完整列表；边表示图上的连接，不单独证明全部共底物可供应。
    一个节点包含多个 C 编号时，通过 reaction_lookup 获取 R 方程来核对。
    该回调只用于消除参与物歧义，单/双向按 KGML reaction.type 创建。
    ec_xml、ko_xml 可提供同图的注释视图，不参与决定连接与方向。
    缺少 reaction 的节点使用经原图核对的 map_supplements 补充或共用连接记录；未补齐时报错。
    graph.graph['map_links'] 保存化合物到其他 map 的链接，不作为化合物反应边。
    graph.graph['diagram_entries'] 保存原图节点及 graphics，不将坐标合并到 Entity。
    边的 ko_ids、diagram_entry_ids 分别保存 KO 和对应的图内反应节点编号。
    """
    root, entries = _read_view(xml, "rn")
    number = root.get("number")
    ecs = _annotations(ec_xml, "ec", root, entries)
    kos = _annotations(ko_xml, "ko", root, entries)
    reactions, supplemented, reused = _complete_reactions(root, entries)

    def annotations(index: dict, entry_id: str, reaction_ids) -> set[str]:
        return set().union(*(index.get(entry_id, {}).get(rid, set()) for rid in reaction_ids))

    # 先解析所有连接，确保源数据和参与物歧义检查通过后再使用共享实体池。
    resolved_reactions = []
    for reaction in reactions:
        entry_id = reaction.get("id")
        reaction_ids = _identifiers(reaction.get("name", ""), "rn:")
        raw_direction = reaction.get("type")
        if raw_direction not in ("irreversible", "reversible"):
            raise ValueError(f"反应缺少有效的图示方向：{reaction_ids}")
        substrate_groups = _participants(reaction, "substrate", entries)
        product_groups = _participants(reaction, "product", entries)
        if any(len(group) > 1 for group in substrate_groups + product_groups):
            if reaction_lookup is None:
                raise ValueError(f"多 C 编号节点需要 R 方程核对：{reaction_ids}")
            resolved = [
                ([rid], *_resolve_participants(substrate_groups, product_groups, reaction_lookup(rid)))
                for rid in reaction_ids
            ]
        else:
            resolved = [(
                reaction_ids, list(dict.fromkeys(group[0] for group in substrate_groups)),
                list(dict.fromkeys(group[0] for group in product_groups)),
            )]
        for ids, substrates, products in resolved:
            resolved_reactions.append((
                entry_id, ids, substrates, products, raw_direction,
                annotations(ecs, entry_id, ids), annotations(kos, entry_id, ids),
            ))

    graph = nx.MultiDiGraph(
        map_id=f"map{number}", name=root.get("title"),
        kgml_id=root.get("name"), image=root.get("image"), link=root.get("link"),
    )
    graph.graph["source_views"] = {
        view: f"{view}{number}" for view, document in (("rn", xml), ("ec", ec_xml), ("ko", ko_xml))
        if document is not None
    }
    graph.graph["supplemented_entry_ids"] = supplemented
    graph.graph["reaction_aliases"] = reused
    graph.graph["map_links"] = _map_links(root, entries)
    graph.graph["diagram_entries"] = {
        entry_id: {
            "type": entry.get("type"), "name": entry.get("name"), "link": entry.get("link"),
            "compound_ids": _identifiers(entry.get("name", ""), "cpd:"),
            "reaction_ids": sorted(_declared_reactions(entry)),
            "pathway_ids": _identifiers(entry.get("name", ""), "path:"),
            "ec_numbers": annotations(ecs, entry_id, _declared_reactions(entry)),
            "ko_ids": annotations(kos, entry_id, _declared_reactions(entry)),
            "graphics": [graphic.attrib.copy() for graphic in entry.findall("graphics")],
        }
        for entry_id, entry in entries.items()
    }
    entity_pool = {} if entities is None else entities

    def add_entity(compound_id: str) -> None:
        if re.fullmatch(r"C[0-9]{5}", compound_id) is None:
            raise ValueError(f"无效的化合物编号：{compound_id}")
        if compound_id not in entity_pool:
            entity_pool[compound_id] = Entity(
                id=compound_id, is_tca=compound_id in TCA_COMPOUND_IDS,
            )
        entity = entity_pool[compound_id]
        if entity.id != compound_id:
            raise ValueError(f"实体字典的键与 Entity.id 不一致：{compound_id}")
        if compound_id in TCA_COMPOUND_IDS:
            entity.is_tca = True
        graph.add_node(compound_id, entity=entity)

    for entry in entries.values():
        if entry.get("type") == "compound":
            for compound_id in _identifiers(entry.get("name", ""), "cpd:"):
                add_entity(compound_id)

    for entry_id, ids, substrates, products, direction, ec_numbers, ko_ids in resolved_reactions:
        for source_id, target_id in product(substrates, products):
            _add_connection(
                graph, source_id, target_id, ids, ec_numbers,
                direction, substrates, products, ko_ids, entry_id,
            )
            if direction == "reversible" and source_id != target_id:
                _add_connection(
                    graph, target_id, source_id, ids, ec_numbers,
                    direction, products, substrates, ko_ids, entry_id,
                )
    return graph
