"""Export a KEGG MultiDiGraph as labelled, machine-readable SVG using Graphviz."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import asdict, is_dataclass
from html import escape
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import networkx as nx

from .logging_config import get_logger


logger = get_logger("svg")

SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
_INK = "#303030"
_BLUE = "#1F497D"
_EDGE_COLORS = {True: "#26804A", False: "#B83D3D", None: _INK}


def _json_value(value):
    """Convert dataclasses and NetworkX containers without discarding attributes."""
    if is_dataclass(value) and not isinstance(value, type):
        return _json_value(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return [_json_value(item) for item in sorted(value, key=repr)]
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"SVG metadata cannot serialize {type(value).__name__}")


def _json_text(value):
    return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"), allow_nan=False)


def _dot_quote(value):
    # JSON string escaping is compatible with quoted DOT labels and identifiers.
    return json.dumps(str(value), ensure_ascii=False)


def _field(value, name, default=None):
    if isinstance(value, Mapping):
        return value.get(name, default)
    return getattr(value, name, default)


def _prepare_export(graph, title=None):
    if not isinstance(graph, nx.MultiDiGraph):
        raise TypeError("graph must be a networkx.MultiDiGraph")
    nodes = sorted(graph.nodes(data=True), key=lambda item: repr(item[0]))
    edges = sorted(graph.edges(keys=True, data=True), key=lambda item: tuple(map(repr, item[:3])))
    terminal = set(graph.graph.get("terminal_nodes", ()))
    map_metadata = graph.graph.get("map_metadata", {})
    compound_graph = graph.graph.get("representation") == "compound"
    records = {}
    clusters = defaultdict(list)
    used_ids = set()
    for key, attributes in nodes:
        entity = attributes.get("entity")
        cid = attributes.get("compound_id") or _field(entity, "id")
        if cid is None:
            raise ValueError(f"Node {key!r} has no entity.id or compound_id")
        map_id = "" if compound_graph else attributes.get("map_id") or graph.graph.get("map_id", "")
        if not compound_graph and not map_id and isinstance(key, tuple) and len(key) == 2:
            map_id = str(key[0])
        if compound_graph and key != cid:
            raise ValueError("Compound network node keys must match their Entity C identifier")
        map_ids = sorted(set(attributes.get("map_ids", ())) | ({str(map_id)} if map_id else set()))
        identifier = "node-" + re.sub(r"[^A-Za-z0-9_.-]", "-", f"{map_id}-{cid}".strip("-"))
        base = identifier
        suffix = 2
        while identifier in used_ids:
            identifier = f"{base}-{suffix}"
            suffix += 1
        used_ids.add(identifier)
        record = {"id": _json_value(key), "svg_id": identifier, "attributes": _json_value(attributes),
                  "compound_id": str(cid), "map_id": str(map_id), "map_ids": map_ids, "terminal": key in terminal}
        records[key] = record
        clusters[str(map_id)].append(key)
    edge_records = []
    for index, (source, target, key, attributes) in enumerate(edges, 1):
        kind = attributes.get("kind", "reaction")
        if kind not in {"reaction", "map_switch"}:
            raise ValueError(f"Unsupported network edge kind: {kind!r}")
        if compound_graph and kind != "reaction":
            raise ValueError("Compound networks must contain reaction edges only")
        connection = attributes.get("connection")
        if kind == "reaction" and connection is None:
            raise ValueError(f"Reaction edge {source!r} -> {target!r} has no connection")
        ecs = sorted(_field(connection, "ec_numbers", ())) if connection is not None else []
        rids = sorted(_field(connection, "reaction_ids", ())) if connection is not None else []
        origins = attributes.get("origins", ())
        map_ids = sorted({str(origin["map_id"]) for origin in origins if origin.get("map_id")})
        if not map_ids and attributes.get("map_id"):
            map_ids = [str(attributes["map_id"])]
        edge_records.append({
            "source": _json_value(source), "target": _json_value(target), "key": _json_value(key),
            "svg_id": f"edge-{index:06d}", "attributes": _json_value(attributes),
            "kind": kind, "ec_numbers": ecs, "reaction_ids": rids,
            "map_id": attributes.get("map_id") or records[source]["map_id"],
            "map_ids": map_ids,
        })
    heading = title or f"{'Merged KEGG compound network' if compound_graph else 'Combined KEGG network'}: {graph.graph.get('target_id', graph.graph.get('map_id', 'pathway'))}"
    subtitle = (f"{len(nodes)} unique compounds | {len(edges)} reaction edges | filled circles: TCA-linked compounds"
                if compound_graph else
                f"{len(nodes)} compound nodes | {len(edges)} directed edges | solid: reactions; dashed: map links")
    lines = ["digraph PathFiller {", "  graph [rankdir=LR, compound=true, newrank=true, splines=spline,",
             "         nodesep=0.5, ranksep=1.0, pad=0.22, bgcolor=white,",
             f"         fontname=Arial, fontsize=20, fontcolor={_dot_quote(_INK)}, labelloc=t,",
             f"         label={_dot_quote(heading + chr(10) + subtitle)}];",
             f"  node [shape=circle, fixedsize=true, width=1.08, height=1.08, fontname=Arial, fontsize=13,",
             f"        color={_dot_quote(_BLUE)}, fontcolor={_dot_quote(_INK)}, penwidth=1.4, style=filled, fillcolor=white];",
             f"  edge [fontname=Arial, fontsize=13, color={_dot_quote(_INK)}, fontcolor={_dot_quote(_INK)},",
             "        penwidth=1.2, arrowsize=0.75];"]
    drawing_groups = [("", list(records))] if compound_graph else sorted(clusters.items())
    for cluster_index, (map_id, keys) in enumerate(drawing_groups):
        info = map_metadata.get(map_id, {})
        virtual_tca = map_id == "map00020" and map_id not in map_metadata and graph.graph.get("map_id") != "map00020"
        if virtual_tca:
            label = "map00020 / TCA link endpoints"
        else:
            label = map_id or "Compound network"
            name = info.get("name") or (graph.graph.get("name") if not map_metadata else None)
            if name:
                label += " / " + str(name)
        if not compound_graph:
            lines += [f"  subgraph cluster_{cluster_index} {{",
                      f"    graph [id={_dot_quote('cluster-' + (map_id or 'network'))}, label={_dot_quote(label)},",
                      f"           color={_dot_quote(_BLUE)}, fontcolor={_dot_quote(_BLUE)}, fontsize=18, penwidth=1.3, margin=22, "
                      f"labelloc={'b' if virtual_tca else 't'}];"]
        for key in keys:
            record = records[key]
            fill = _BLUE if record["terminal"] else "white"
            font = "white" if record["terminal"] else _INK
            width = 2.5 if record["compound_id"] == graph.graph.get("target_id") else 1.4
            lines.append(f"    {_dot_quote(record['svg_id'])} [id={_dot_quote(record['svg_id'])}, label={_dot_quote(record['compound_id'])}, "
                         f"fillcolor={_dot_quote(fill)}, fontcolor={_dot_quote(font)}, penwidth={width}];")
        if not compound_graph:
            lines.append("  }")
    for (source, target, _, attributes), record in zip(edges, edge_records, strict=True):
        switch = record["kind"] == "map_switch"
        assessed = _field(attributes.get("connection"), "is_endogenous")
        if assessed not in _EDGE_COLORS:
            raise ValueError("Connection.is_endogenous must be True, False or None")
        color = _BLUE if switch else _EDGE_COLORS[assessed]
        labels = ["map link"] if switch else [f"EC {ec}" for ec in record["ec_numbers"]] or ["EC unknown"]
        label = ('<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2" BGCOLOR="white">' +
                 ''.join(f'<TR><TD><FONT FACE="Arial" POINT-SIZE="13" COLOR="{color}">{escape(text)}</FONT></TD></TR>'
                         for text in labels) + '</TABLE>>')
        style = "dashed" if switch else "solid"
        lines.append(f"  {_dot_quote(records[source]['svg_id'])} -> {_dot_quote(records[target]['svg_id'])} "
                     f"[id={_dot_quote(record['svg_id'])}, label={label}, style={style}, "
                     f"color={_dot_quote(color)}, fontcolor={_dot_quote(color)}];")
    lines.append("}")
    metadata = {"schema": "pathfiller.network.svg.v2" if compound_graph else "pathfiller.network.svg.v1", "directed": True, "multigraph": True,
                "graph": _json_value(graph.graph), "nodes": list(records.values()), "edges": edge_records}
    if compound_graph:
        metadata.update(representation="compound", schema_version=2)
    return "\n".join(lines) + "\n", metadata, records, edge_records, heading


def _annotate_svg(svg_text, metadata, nodes, edges, heading):
    root = ET.fromstring(svg_text)
    root.set("data-node-count", str(len(nodes)))
    root.set("data-edge-count", str(len(edges)))
    if metadata.get("representation") == "compound":
        root.set("data-representation", "compound")
        root.set("data-schema-version", "2")
    root.set("role", "img")
    root.set("aria-labelledby", "pathfiller-title pathfiller-description")
    root.insert(0, ET.Element(f"{{{SVG}}}title", {"id": "pathfiller-title"}))
    root[0].text = heading
    root.insert(1, ET.Element(f"{{{SVG}}}desc", {"id": "pathfiller-description"}))
    root[1].text = ("Each KEGG compound is drawn once as a circle. Arrows show reactions and EC identifiers; "
                    "map memberships and every original map edge are retained as provenance in SVG metadata."
                    if metadata.get("representation") == "compound" else
                    "Circular compound nodes labelled with KEGG C identifiers. Reaction arrows show EC identifiers; "
                    "dashed arrows are explicit map links. Full node, edge and diagram attributes are embedded as JSON metadata.")
    element = ET.Element(f"{{{SVG}}}metadata", {"id": "pathfiller-network"})
    element.text = json.dumps(metadata, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    root.insert(2, element)
    groups = {element.get("id"): element for element in root.iter(f"{{{SVG}}}g")}
    drawing = next((group for group in root.findall(f"{{{SVG}}}g")
                    if group.get("class") == "graph"), None)
    if drawing is None:
        raise RuntimeError("Graphviz SVG contains no graph drawing group")
    labels_layer = ET.Element(f"{{{SVG}}}g", {"id": "pathfiller-edge-labels", "class": "connection-labels"})
    for record in nodes.values():
        group = groups.get(record["svg_id"])
        if group is None:
            raise RuntimeError(f"Graphviz omitted SVG node {record['svg_id']}")
        group.set("class", "node compound")
        group.set("data-node-key", _json_text(record["id"]))
        group.set("data-compound-id", record["compound_id"])
        group.set("data-map-id", record["map_id"])
        group.set("data-map-ids", _json_text(record["map_ids"]))
        group.set("data-is-terminal", str(record["terminal"]).lower())
        _set_title(group, f"{record['compound_id']} / {', '.join(record['map_ids'])}\n" + _json_text(record["attributes"]))
        for ellipse in list(group.findall(f"{{{SVG}}}ellipse")):
            if abs(float(ellipse.get("rx", "0")) - float(ellipse.get("ry", "0"))) < 0.01:
                ellipse.tag = f"{{{SVG}}}circle"
                ellipse.set("r", ellipse.attrib.pop("rx"))
                ellipse.attrib.pop("ry", None)
    for record in edges:
        group = groups.get(record["svg_id"])
        if group is None:
            raise RuntimeError(f"Graphviz omitted SVG edge {record['svg_id']}")
        group.set("class", "edge " + record["kind"])
        for attr, key in (("data-source", "source"), ("data-target", "target"), ("data-edge-key", "key"),
                          ("data-ec-numbers", "ec_numbers"), ("data-reaction-ids", "reaction_ids")):
            group.set(attr, _json_text(record[key]))
        group.set("data-kind", record["kind"])
        group.set("data-map-id", str(record["map_id"]))
        group.set("data-map-ids", _json_text(record["map_ids"]))
        group.set("data-origins", _json_text(record["attributes"].get("origins", [])))
        group.set("data-network-edge-key", _json_text(record["key"]))
        for attr, key in (("data-diagram-entry-ids", "diagram_entry_ids"), ("data-ko-ids", "ko_ids")):
            values = record["attributes"].get(key)
            if values is None:
                values = sorted({value for origin in record["attributes"].get("origins", [])
                                 for value in origin.get(key, [])})
            group.set(attr, _json_text(values))
        group.set("data-is-endogenous", _json_text(_field(record["attributes"].get("connection"), "is_endogenous")))
        _set_title(group, _json_text(record))
        # Graphviz emits each edge's path followed by its label. A later edge
        # can therefore cross an earlier label. Keep editable native text in
        # the original edge group, but render every label together in the last
        # SVG layer so its white background masks all intersecting paths.
        label_items = [child for child in group if child.tag == f"{{{SVG}}}text" or (
            child.tag == f"{{{SVG}}}polygon" and child.get("fill") == "white" and child.get("stroke") == "none"
        )]
        if label_items:
            definitions = ET.SubElement(group, f"{{{SVG}}}defs")
            label_id = record["svg_id"] + "-label"
            definition = ET.SubElement(definitions, f"{{{SVG}}}g", {"id": label_id, "data-edge-id": record["svg_id"]})
            for item in label_items:
                group.remove(item)
                definition.append(item)
            instance = ET.SubElement(labels_layer, f"{{{SVG}}}use", {
                "href": "#" + label_id,
                "{http://www.w3.org/1999/xlink}href": "#" + label_id,
                "data-edge-id": record["svg_id"], "class": "connection-label",
            })
            ET.SubElement(instance, f"{{{SVG}}}title").text = _json_text(record)
    drawing.append(labels_layer)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def _set_title(group, text):
    title = group.find(f"{{{SVG}}}title")
    if title is None:
        title = ET.Element(f"{{{SVG}}}title")
        group.insert(0, title)
    title.text = text


def export_network_svg(graph: nx.MultiDiGraph, output_path, *, title: str | None = None) -> Path:
    """Export all network nodes and edges as SVG without modifying the input graph.

    Graphviz's ``dot`` executable must be installed and available on PATH.
    Nodes use circular C labels, reactions show every EC annotation, and map
    links are dashed for legacy map-state networks. A graph marked with
    ``representation='compound'`` is flat: each C identifier is drawn once,
    with original map memberships and merged edge origins retained in metadata.
    The SVG includes ``metadata#pathfiller-network`` JSON and
    stable node/edge groups with ``data-*`` attributes for later applications.
    """
    source, metadata, nodes, edges, heading = _prepare_export(graph, title)
    logger.info("开始导出 SVG：%d 个节点、%d 条边，正在计算 Graphviz 布局", len(nodes), len(edges))
    executable = shutil.which("dot")
    if executable is None:
        raise RuntimeError("SVG export requires Graphviz 'dot' on PATH; install Graphviz first")
    try:
        result = subprocess.run([executable, "-Tsvg"], input=source, encoding="utf-8",
                                capture_output=True, check=False, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise RuntimeError(f"Graphviz SVG export failed: {error}") from error
    if result.returncode != 0:
        raise RuntimeError(f"Graphviz SVG export failed: {result.stderr.strip()}")
    try:
        svg = _annotate_svg(result.stdout, metadata, nodes, edges, heading)
    except ET.ParseError as error:
        raise RuntimeError("Graphviz returned invalid SVG XML") from error
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(svg)
    logger.info("SVG 导出完成：%s", output.resolve())
    return output
