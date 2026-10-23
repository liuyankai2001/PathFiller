"""Export a KEGG MultiDiGraph as labelled, machine-readable SVG using Graphviz."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
import colorsys
from dataclasses import asdict, is_dataclass
from html import escape
from heapq import heappop, heappush
from itertools import count
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import textwrap
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
_ROUTE_COLORS = {1: "#2563EB", 2: "#D97706", 3: "#7C3AED"}
_BACKGROUND_EDGE = "#606060"


def route_colors_for(ranks) -> dict[int, str]:
    """Return stable route colours without imposing a candidate-count limit."""
    ranks = set(ranks)
    if any(type(rank) is not int or rank < 1 for rank in ranks):
        raise ValueError("Candidate ranks must be positive integers")
    colours = {}
    used = set(_ROUTE_COLORS.values())
    for rank in sorted(ranks):
        if rank in _ROUTE_COLORS:
            colour = _ROUTE_COLORS[rank]
        else:
            hue = (0.46 + ((rank - 4) % 1000003) * 0.618033988749895) % 1
            # A fixed-saturation HSV ring has only finitely many rounded RGB
            # colours. Vary saturation/value on collisions rather than looping
            # forever once a large candidate set fills that ring.
            for trial in range(64):
                saturation = 0.70 if trial == 0 else 0.55 + 0.25 * (
                    ((rank % 1000003) * 0.754877666 + trial * 0.618033989) % 1)
                value = 0.66 if trial == 0 else 0.52 + 0.22 * (
                    ((rank % 1000003) * 0.569840296 + trial * 0.754877666) % 1)
                rgb = colorsys.hsv_to_rgb((hue + trial * 0.037) % 1, saturation, value)
                colour = "#" + "".join(f"{round(value * 255):02X}" for value in rgb)
                if colour not in used:
                    break
            else:
                # The odd stride visits every RGB value before repeating.
                # A bounded scan keeps colour assignment finite even for
                # pathological inputs; no candidate is dropped for its colour.
                for offset in range(min(len(used) + 1, 1 << 24)):
                    colour = f"#{((rank + offset) * 0x9E3779) & 0xFFFFFF:06X}"
                    if colour not in used:
                        break
        colours[rank] = colour
        used.add(colour)
    return colours


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


def _name_lines(name):
    """Wrap visible labels without dropping any part of a compound name."""
    return textwrap.wrap(name, width=28, break_long_words=True, break_on_hyphens=True)


def _supplement_ecs(value):
    """Normalize the selected EC set without altering reaction annotations."""
    if isinstance(value, (str, bytes)):
        raise ValueError("Supplement ECs must be a collection of EC identifiers")
    result = sorted(set(value))
    if any(not isinstance(ec, str) or not ec for ec in result):
        raise ValueError("Supplement ECs must contain nonempty EC identifiers")
    return result


def _supplement_label(ecs):
    return "Supplement ECs: " + (", ".join(ecs) or "none")


def _prepare_export(graph, title=None, *, overview=False, highlighted_paths=None):
    if not isinstance(graph, nx.MultiDiGraph):
        raise TypeError("graph must be a networkx.MultiDiGraph")
    if type(overview) is not bool:
        raise TypeError("overview must be a boolean")
    if not overview and highlighted_paths is not None:
        raise ValueError("highlighted_paths requires overview=True")
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
        name = _field(entity, "name")
        if name is not None and not isinstance(name, str):
            raise ValueError(f"Node {key!r} has a compound name that is not text")
        name = name.strip() if name else None
        record = {"id": _json_value(key), "svg_id": identifier, "attributes": _json_value(attributes),
                  "compound_id": str(cid), "compound_name": name or None,
                  "map_id": str(map_id), "map_ids": map_ids, "terminal": key in terminal}
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
    if overview:
        return _prepare_overview(graph, nodes, edges, records, edge_records, title, highlighted_paths)
    heading = title or f"{'Merged KEGG compound network' if compound_graph else 'Combined KEGG network'}: {graph.graph.get('target_id', graph.graph.get('map_id', 'pathway'))}"
    subtitle = (f"{len(nodes)} unique compounds | {len(edges)} reaction edges | filled circles: TCA-linked compounds"
                if compound_graph else
                f"{len(nodes)} compound nodes | {len(edges)} directed edges | solid: reactions; dashed: map links")
    scheme_mode = graph.graph.get("candidate_unit") == "extra_ec_set"
    if scheme_mode:
        supplement = _supplement_label(_supplement_ecs(graph.graph.get("extra_ecs", ())))
        subtitle += "\n" + "\n".join(textwrap.wrap(supplement, width=88, break_long_words=False, break_on_hyphens=False))
    node_spacing = "0.16" if any(record["compound_name"] for record in records.values()) else "0.5"
    lines = ["digraph PathFiller {", "  graph [rankdir=LR, compound=true, newrank=true, splines=spline,",
             f"         nodesep={node_spacing}, ranksep=1.0, pad=0.22, bgcolor=white,",
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
            if record["compound_name"]:
                # A layout-only plain node reserves the complete name's bounds.
                # In LR layouts a same-rank invisible edge orders it below the
                # circle while preserving the real reaction's attachment point.
                name_id = "compound-name-" + record["svg_id"]
                label = _dot_quote("\n".join(_name_lines(record["compound_name"])))
                lines.append(f"    {_dot_quote(name_id)} [id={_dot_quote(name_id)}, shape=plain, fixedsize=false, "
                             f"width=0, height=0, margin=0, style=\"\", fontcolor={_dot_quote(_INK)}, label={label}];")
                lines.append(f"    {{ rank=same; {_dot_quote(record['svg_id'])}; {_dot_quote(name_id)}; "
                             f"{_dot_quote(record['svg_id'])} -> {_dot_quote(name_id)} [style=invis, weight=100]; }}")
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
    if scheme_mode:
        metadata["candidate_unit"] = "extra_ec_set"
        metadata["extra_ecs"] = _supplement_ecs(graph.graph.get("extra_ecs", ()))
    return "\n".join(lines) + "\n", metadata, records, edge_records, heading


def _prepare_overview(graph, nodes, edges, records, edge_records, title, highlighted_paths):
    """Create a flat sfdp drawing while keeping all business nodes and edges."""
    if any(key != record["compound_id"] for key, record in records.items()):
        raise ValueError("Overview nodes must use their compound C identifiers as keys")
    if highlighted_paths is None:
        highlighted_paths = {}
    if not isinstance(highlighted_paths, Mapping):
        raise TypeError("highlighted_paths must map candidate ranks to ordered edge sequences")
    colours = route_colors_for(highlighted_paths)
    members = defaultdict(set)
    node_members = defaultdict(set)
    route_records = []
    scheme_mode = graph.graph.get("candidate_unit") == "extra_ec_set"
    scheme_ecs = graph.graph.get("candidate_extra_ecs", {})
    if scheme_mode and not isinstance(scheme_ecs, Mapping):
        raise ValueError("candidate_extra_ecs must map scheme ranks to selected EC sets")
    for rank in sorted(highlighted_paths):
        selected, previous = [], None
        for edge in highlighted_paths[rank]:
            if not isinstance(edge, (tuple, list)) or len(edge) != 3:
                raise ValueError("Highlighted edges must contain source, target and merged edge key")
            source, target, key = edge
            try:
                exists = graph.has_edge(source, target, key)
            except TypeError as error:
                raise ValueError("Highlighted edge identifiers must be hashable") from error
            if not exists or (previous is not None and source != previous):
                raise ValueError(f"Candidate {rank} contains a missing or discontinuous network edge")
            selected.append((source, target, key))
            members[source, target, key].add(rank)
            node_members[source].add(rank)
            node_members[target].add(rank)
            previous = target
        if not selected:
            target = graph.graph.get("target_id")
            if target is None and len(graph) == 1:
                target = next(iter(graph))
            if target in graph:
                node_members[target].add(rank)
        route = {"rank": rank, "color": colours[rank], "reaction_steps": len(selected),
                 "edges": _json_value(selected)}
        if scheme_mode:
            if rank not in scheme_ecs and str(rank) not in scheme_ecs:
                raise ValueError(f"Missing supplementary EC set for scheme {rank}")
            route["extra_ecs"] = _supplement_ecs(scheme_ecs.get(rank, scheme_ecs.get(str(rank))))
        route_records.append(route)
    for key, record in records.items():
        record["candidate_ranks"] = sorted(node_members[key])
    for (source, target, key, attributes), record in zip(edges, edge_records, strict=True):
        if record["kind"] != "reaction":
            raise ValueError("Overview networks must contain reaction edges only")
        assessed = _field(attributes.get("connection"), "is_endogenous")
        if assessed not in _EDGE_COLORS:
            raise ValueError("Connection.is_endogenous must be True, False or None")
        record["candidate_ranks"] = sorted(members[source, target, key])
        record["route_colors"] = [colours[rank] for rank in record["candidate_ranks"]]
    heading = title or f"Merged KEGG compound network: {graph.graph.get('target_id', 'pathway')}"
    subtitle = f"{len(nodes)} unique compounds | {len(edges)} reaction edges | hover: names and EC"
    lines = ["digraph PathFiller {", "  graph [overlap=prism, sep=\"+14\", splines=line, start=22,",
             "         pad=0.22, bgcolor=white, fontname=Arial, fontsize=20, labelloc=t,",
             f"         fontcolor={_dot_quote(_INK)}, label={_dot_quote(heading + chr(10) + subtitle)}];",
             f"  node [shape=circle, fixedsize=true, width=1.08, height=1.08, fontname=Arial, fontsize=13,",
             f"        color={_dot_quote(_BLUE)}, fontcolor={_dot_quote(_INK)}, penwidth=1.4, style=filled, fillcolor=white];",
             "  edge [arrowsize=0.75];"]
    for record in records.values():
        fill = _BLUE if record["terminal"] else "white"
        font = "white" if record["terminal"] else _INK
        width = 2.5 if record["compound_id"] == graph.graph.get("target_id") else 1.4
        lines.append(f"  {_dot_quote(record['svg_id'])} [id={_dot_quote(record['svg_id'])}, "
                     f"label={_dot_quote(record['compound_id'])}, fillcolor={_dot_quote(fill)}, "
                     f"fontcolor={_dot_quote(font)}, penwidth={width}];")
    for (source, target, _, _), record in zip(edges, edge_records, strict=True):
        colour = ":".join(record["route_colors"]) or _BACKGROUND_EDGE
        width = 2.2 if record["candidate_ranks"] else 0.8
        lines.append(f"  {_dot_quote(records[source]['svg_id'])} -> {_dot_quote(records[target]['svg_id'])} "
                     f"[id={_dot_quote(record['svg_id'])}, color={_dot_quote(colour)}, penwidth={width}];")
    lines.append("}")
    metadata = {"schema": "pathfiller.network.svg.v2", "schema_version": 2, "representation": "compound",
                "directed": True, "multigraph": True, "overview": True, "layout_engine": "sfdp",
                "graph": _json_value(graph.graph), "nodes": list(records.values()), "edges": edge_records,
                "route_colors": {str(rank): colour for rank, colour in colours.items()}, "routes": route_records}
    if scheme_mode:
        metadata["candidate_unit"] = "extra_ec_set"
    return "\n".join(lines) + "\n", metadata, records, edge_records, heading


def _annotate_svg(svg_text, metadata, nodes, edges, heading):
    root = ET.fromstring(svg_text)
    overview = metadata.get("overview", False)
    root.set("data-node-count", str(len(nodes)))
    root.set("data-edge-count", str(len(edges)))
    if metadata.get("representation") == "compound":
        root.set("data-representation", "compound")
        root.set("data-schema-version", "2")
    root.set("role", "img")
    if overview:
        root.set("data-overview", "true")
        root.set("data-route-colors", _json_text(metadata["route_colors"]))
    if metadata.get("candidate_unit") == "extra_ec_set":
        root.set("data-candidate-unit", "extra_ec_set")
        if not overview:
            root.set("data-extra-ecs", _json_text(metadata["extra_ecs"]))
    root.set("aria-labelledby", "pathfiller-title pathfiller-description")
    root.insert(0, ET.Element(f"{{{SVG}}}title", {"id": "pathfiller-title"}))
    root[0].text = heading
    root.insert(1, ET.Element(f"{{{SVG}}}desc", {"id": "pathfiller-description"}))
    root[1].text = ("Each KEGG compound is drawn once as a circle. Arrows show reactions and EC identifiers; "
                    "map memberships and every original map edge are retained as provenance in SVG metadata."
                    if metadata.get("representation") == "compound" else
                    "Circular compound nodes labelled with KEGG C identifiers. Reaction arrows show EC identifiers; "
                    "dashed arrows are explicit map links. Full node, edge and diagram attributes are embedded as JSON metadata.")
    if overview:
        root[1].text = ("Whole merged metabolic network. Hover compounds for names and reactions for EC, R, candidate ranks "
                        "and host matching. Coloured parallel strokes mark shared candidate reactions.")
    element = ET.Element(f"{{{SVG}}}metadata", {"id": "pathfiller-network"})
    element.text = json.dumps(metadata, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    root.insert(2, element)
    groups = {element.get("id"): element for element in root.iter(f"{{{SVG}}}g")}
    parents = {child: parent for parent in root.iter() for child in parent}
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
        if record["compound_name"]:
            group.set("data-compound-name", record["compound_name"])
        group.set("data-map-id", record["map_id"])
        group.set("data-map-ids", _json_text(record["map_ids"]))
        group.set("data-is-terminal", str(record["terminal"]).lower())
        display = record['compound_id'] + (" / " + record["compound_name"] if record["compound_name"] else "")
        if overview:
            group.set("data-candidate-ranks", _json_text(record["candidate_ranks"]))
            _set_title(group, record["compound_id"] + "\n" + (record["compound_name"] or "名称未获取"))
        else:
            _set_title(group, f"{display} / {', '.join(record['map_ids'])}\n" + _json_text(record["attributes"]))
        for ellipse in list(group.findall(f"{{{SVG}}}ellipse")):
            if abs(float(ellipse.get("rx", "0")) - float(ellipse.get("ry", "0"))) < 0.01:
                ellipse.tag = f"{{{SVG}}}circle"
                ellipse.set("r", ellipse.attrib.pop("rx"))
                ellipse.attrib.pop("ry", None)
        if record["compound_name"] and not overview:
            name_group = groups.get("compound-name-" + record["svg_id"])
            if name_group is None:
                raise RuntimeError(f"Graphviz omitted compound name for SVG node {record['svg_id']}")
            # Keep text in the same editable/inspectable chemical node, without
            # exposing the auxiliary layout node as a second business entity.
            parents[name_group].remove(name_group)
            name_group.set("class", "compound-name")
            name_group.set("data-compound-id", record["compound_id"])
            name_group.set("data-compound-name", record["compound_name"])
            _set_title(name_group, record["compound_name"])
            group.append(name_group)
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
        if overview:
            group.set("data-candidate-ranks", _json_text(record["candidate_ranks"]))
            group.set("data-route-colors", _json_text(record["route_colors"]))
            matched = _field(record["attributes"].get("connection"), "is_endogenous")
            unit = "方案" if metadata.get("candidate_unit") == "extra_ec_set" else "候选"
            _set_title(group, "EC: " + (", ".join(record["ec_numbers"]) or "未知") +
                       "\nR: " + (", ".join(record["reaction_ids"]) or "未知") +
                       f"\n{unit}: " + (", ".join(map(str, record["candidate_ranks"])) or "无") +
                       "\n宿主匹配: " + ("是" if matched else ("否" if matched is False else "未评估")))
            continue
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
    if overview:
        _route_overview_highlights(root, drawing, metadata, nodes, edges, groups)
        # Base reactions, highlighted reactions, and compounds are explicit
        # layers, so route colours cannot be hidden by a later background edge.
        for record in edges:
            group = groups[record["svg_id"]]
            parents[group].remove(group)
        for record in nodes.values():
            group = groups[record["svg_id"]]
            parents[group].remove(group)
        for highlighted in (False, True):
            for record in edges:
                if bool(record["candidate_ranks"]) == highlighted:
                    drawing.append(groups[record["svg_id"]])
        for record in nodes.values():
            drawing.append(groups[record["svg_id"]])
        _overview_legend(root, metadata)
        element.text = json.dumps(metadata, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
    else:
        drawing.append(labels_layer)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


_GEOMETRY_ERROR = 0.03
_ROUTE_CLEARANCE = 2.0


def _segment_distance_squared(point, start, end):
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = dx * dx + dy * dy
    t = 0 if length == 0 else max(0.0, min(1.0, ((point[0] - start[0]) * dx + (point[1] - start[1]) * dy) / length))
    return (point[0] - start[0] - t * dx) ** 2 + (point[1] - start[1] - t * dy) ** 2


def _svg_path_lines(value):
    """Flatten Graphviz cubic curves with a conservative distance error bound.

    Every cubic lies in its control-point convex hull. Once both inner
    controls are within GEOMETRY_ERROR of the chord, the entire curve is
    inside that chord's error tube; circle tests expand by that error.
    """
    tokens = re.findall(r"[MLCZmlcz]|[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?", value)
    result, position, command, index = [], None, None, 0

    def flatten(a, b, c, d):
        pending = [(a, b, c, d)]
        while pending:
            a, b, c, d = pending.pop()
            if max(_segment_distance_squared(b, a, d), _segment_distance_squared(c, a, d)) <= _GEOMETRY_ERROR ** 2:
                result.append((a, d))
                continue
            ab, bc, cd = tuple((x + y) / 2 for x, y in zip(a, b)), tuple((x + y) / 2 for x, y in zip(b, c)), tuple((x + y) / 2 for x, y in zip(c, d))
            abc, bcd = tuple((x + y) / 2 for x, y in zip(ab, bc)), tuple((x + y) / 2 for x, y in zip(bc, cd))
            middle = tuple((x + y) / 2 for x, y in zip(abc, bcd))
            pending.append((middle, bcd, cd, d))
            pending.append((a, ab, abc, middle))

    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
        if command == "M":
            position = (float(tokens[index]), float(tokens[index + 1]))
            index += 2
            command = "L"
        elif command == "L":
            end = (float(tokens[index]), float(tokens[index + 1]))
            result.append((position, end))
            position, index = end, index + 2
        elif command == "C":
            points = [(float(tokens[i]), float(tokens[i + 1])) for i in range(index, index + 6, 2)]
            flatten(position, *points)
            position, index = points[-1], index + 6
        else:
            raise RuntimeError(f"Unsupported Graphviz SVG path command: {command}")
    return result


class _CircleIndex:
    """A spatial bucket index with exact line-segment/circle collision tests."""

    def __init__(self, circles):
        self.circles = list(circles)
        self.buckets = defaultdict(list)
        self.radius = max((circle[2] for circle in self.circles), default=0)
        for circle in self.circles:
            self.buckets[math.floor(circle[0] / 128), math.floor(circle[1] / 128)].append(circle)

    def nearby(self, start, end, extra=0.0):
        radius = self.radius + extra
        left, right = math.floor((min(start[0], end[0]) - radius) / 128), math.floor((max(start[0], end[0]) + radius) / 128)
        top, bottom = math.floor((min(start[1], end[1]) - radius) / 128), math.floor((max(start[1], end[1]) + radius) / 128)
        if (right - left + 1) * (bottom - top + 1) > len(self.circles):
            return self.circles
        return [circle for x in range(left, right + 1) for y in range(top, bottom + 1)
                for circle in self.buckets.get((x, y), ())]

    def clear(self, start, end, extra=0.0):
        return all(_segment_distance_squared(circle[:2], start, end) > (circle[2] + extra) ** 2
                   for circle in self.nearby(start, end, extra))


def _polygon_points(element):
    values = [float(value) for value in re.findall(r"[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?", element.get("points", ""))]
    return list(zip(values[::2], values[1::2]))


def _point_in_polygon(point, polygon):
    inside = False
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        if ((start[1] > point[1]) != (end[1] > point[1]) and
                point[0] < (end[0] - start[0]) * (point[1] - start[1]) / (end[1] - start[1]) + start[0]):
            inside = not inside
    return inside


def _highlight_clear(group, index):
    for path in group.findall(f"{{{SVG}}}path"):
        extra = float(path.get("stroke-width", "1")) / 2 + _ROUTE_CLEARANCE + _GEOMETRY_ERROR
        if any(not index.clear(a, b, extra) for a, b in _svg_path_lines(path.get("d", ""))):
            return False
    for polygon in group.findall(f"{{{SVG}}}polygon"):
        points = _polygon_points(polygon)
        extra = float(polygon.get("stroke-width", "1")) / 2 + _ROUTE_CLEARANCE
        if any(not index.clear(a, b, extra) for a, b in zip(points, points[1:] + points[:1])):
            return False
        if any(_point_in_polygon(circle[:2], points) for circle in index.circles):
            return False
    return True


def _route_grid(start, goal, index, bounds):
    """Deterministic eight-neighbour A*, followed by exact visibility shortcuts."""
    if index.clear(start, goal):
        return [start, goal]
    for step in (12.0, 6.0, 3.0):
        xmin = math.floor((bounds[0] - start[0]) / step)
        xmax = math.ceil((bounds[2] - start[0]) / step)
        ymin = math.floor((bounds[1] - start[1]) / step)
        ymax = math.ceil((bounds[3] - start[1]) / step)
        serial, queue, best, previous, clear_cache = count(), [], {(0, 0): 0.0}, {}, {}

        def position(key):
            return start[0] + key[0] * step, start[1] + key[1] * step

        def distance(point):
            return math.hypot(point[0] - goal[0], point[1] - goal[1])

        heappush(queue, (distance(start), next(serial), 0.0, (0, 0)))
        end = None
        while queue:
            _, _, cost, key = heappop(queue)
            if cost != best.get(key):
                continue
            point = position(key)
            if distance(point) <= step * 2 and index.clear(point, goal):
                end = key
                break
            for dx, dy in ((-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)):
                following = key[0] + dx, key[1] + dy
                if not (xmin <= following[0] <= xmax and ymin <= following[1] <= ymax):
                    continue
                next_cost = cost + step * (math.sqrt(2) if dx and dy else 1)
                if next_cost >= best.get(following, float("inf")):
                    continue
                target = position(following)
                if following not in clear_cache:
                    clear_cache[following] = index.clear(target, target)
                if not clear_cache[following] or not index.clear(point, target):
                    continue
                best[following], previous[following] = next_cost, key
                heappush(queue, (next_cost + distance(target), next(serial), next_cost, following))
        if end is None:
            continue
        keys = [end]
        while keys[-1] != (0, 0):
            keys.append(previous[keys[-1]])
        raw = [position(key) for key in reversed(keys)] + [goal]
        simplified, cursor = [raw[0]], 0
        while cursor < len(raw) - 1:
            following = len(raw) - 1
            while following > cursor + 1 and not index.clear(raw[cursor], raw[following]):
                following -= 1
            simplified.append(raw[following])
            cursor = following
        return simplified
    raise RuntimeError("Unable to route a highlighted reaction around compound nodes")


def _offset_polyline(points, offset):
    """Bevel joins keep the whole offset line inside the bounded route tube."""
    normals = []
    for start, end in zip(points, points[1:]):
        length = math.hypot(end[0] - start[0], end[1] - start[1])
        if length < 1e-8:
            raise RuntimeError("Degenerate highlighted reaction geometry")
        normals.append((-(end[1] - start[1]) / length, (end[0] - start[0]) / length))
    result = [(points[0][0] + offset * normals[0][0], points[0][1] + offset * normals[0][1])]
    for i in range(1, len(points) - 1):
        for normal in (normals[i - 1], normals[i]):
            point = points[i][0] + offset * normal[0], points[i][1] + offset * normal[1]
            if math.dist(point, result[-1]) > 1e-8:
                result.append(point)
    result.append((points[-1][0] + offset * normals[-1][0], points[-1][1] + offset * normals[-1][1]))
    return result


def _free_port(circle, initial, index, *, arrow=False):
    angle = math.atan2(initial[1] - circle[1], initial[0] - circle[0])
    for offset in [0] + [sign * i * math.pi / 48 for i in range(1, 49) for sign in (1, -1)]:
        direction = math.cos(angle + offset), math.sin(angle + offset)
        point = (circle[0] + direction[0] * (circle[2] + (12 if arrow else 1)),
                 circle[1] + direction[1] * (circle[2] + (12 if arrow else 1)))
        tip = circle[0] + direction[0] * (circle[2] + 0.75), circle[1] + direction[1] * (circle[2] + 0.75)
        approach = point[0] + direction[0] * 8, point[1] + direction[1] * 8
        if index.clear(point, approach if arrow else point) and (not arrow or index.clear(point, tip)):
            return point, approach, tip, direction
    raise RuntimeError("No unobstructed port for highlighted compound reaction")


def _route_overview_highlights(root, drawing, metadata, nodes, edges, groups):
    """Keep native unobstructed curves; route colliding highlight bundles only."""
    circles = {}
    for key, record in nodes.items():
        circle = groups[record["svg_id"]].find(f"{{{SVG}}}circle")
        if circle is None:
            raise RuntimeError("Overview compound node has no circular geometry")
        circles[key] = (float(circle.get("cx")), float(circle.get("cy")), float(circle.get("r")))
    selected = [edge for edge in edges if edge["candidate_ranks"]]
    rerouted, family_count = [], defaultdict(int)
    for record in selected:
        source, target = record["source"], record["target"]
        group = groups[record["svg_id"]]
        original = _CircleIndex([circle for cid, circle in circles.items() if cid not in {source, target}])
        family = tuple(sorted((source, target)))
        lane = family_count[family]
        family_count[family] += 1
        if _highlight_clear(group, original):
            record["geometry_rerouted"] = False
            continue
        paths = group.findall(f"{{{SVG}}}path")
        if not paths:
            raise RuntimeError("Highlighted reaction has no SVG path")
        flattened = [_svg_path_lines(path.get("d", "")) for path in paths]
        start = tuple(sum(lines[0][0][i] for lines in flattened) / len(paths) for i in (0, 1))
        end = tuple(sum(lines[-1][1][i] for lines in flattened) / len(paths) for i in (0, 1))
        width = max(float(path.get("stroke-width", "1")) for path in paths)
        spacing = max(2.6, width + 0.4)
        half = (len(paths) - 1) * spacing / 2
        arrow_half = max(5.0, half + width / 2 + 0.5)
        # The envelope covers both all parallel stroke centres and the arrow.
        envelope = max(half + width / 2, arrow_half + width / 2)
        padding = _ROUTE_CLEARANCE + envelope + lane * 3.5 + 0.15
        obstacles = _CircleIndex([(x, y, radius + padding) for x, y, radius in original.circles])
        start, _, _, _ = _free_port(circles[source], start, obstacles)
        body_end, approach, tip, direction = _free_port(circles[target], end, obstacles, arrow=True)
        all_points = list(circles.values())
        bounds = (min(x - radius for x, _, radius in all_points) - 100,
                  min(y - radius for _, y, radius in all_points) - 100,
                  max(x + radius for x, _, radius in all_points) + 100,
                  max(y + radius for _, y, radius in all_points) + 100)
        road = _route_grid(start, approach, obstacles, bounds)
        if math.dist(road[-1], body_end) > 1e-8:
            road.append(body_end)
        for index, path in enumerate(paths):
            points = _offset_polyline(road, (index - (len(paths) - 1) / 2) * spacing)
            path.set("d", "M" + " L".join(f"{x:.5f},{y:.5f}" for x, y in points))
            path.set("stroke-linejoin", "bevel")
        polygons = group.findall(f"{{{SVG}}}polygon")
        arrow = polygons[0] if polygons else ET.SubElement(group, f"{{{SVG}}}polygon")
        normal = -direction[1], direction[0]
        triangle = [(body_end[0] + normal[0] * arrow_half, body_end[1] + normal[1] * arrow_half),
                    tip, (body_end[0] - normal[0] * arrow_half, body_end[1] - normal[1] * arrow_half)]
        arrow.set("points", " ".join(f"{x:.5f},{y:.5f}" for x, y in [*triangle, triangle[0]]))
        colour = paths[0].get("stroke", _INK)
        arrow.set("fill", colour)
        arrow.set("stroke", colour)
        arrow.set("stroke-width", str(width))
        for extra in polygons[1:]:
            group.remove(extra)
        if not _highlight_clear(group, original):
            raise RuntimeError(f"Highlighted reaction {record['svg_id']} still intersects a non-endpoint compound")
        group.set("data-geometry-rerouted", "true")
        record["geometry_rerouted"] = True
        rerouted.append(record["svg_id"])
    metadata["geometry_routing"] = {"method": "bounded-cubic-check-and-grid-astar", "clearance": _ROUTE_CLEARANCE,
                                    "curve_error_bound": _GEOMETRY_ERROR, "checked_highlighted_edges": len(selected),
                                    "rerouted_edge_ids": rerouted, "remaining_collisions": 0}
    # Keep any outward detour inside the SVG without moving nodes or the legend.
    if rerouted:
        translate = re.search(r"translate\(([-+\d.eE]+)[ ,]+([-+\d.eE]+)\)", drawing.get("transform", ""))
        tx, ty = map(float, translate.groups()) if translate else (0.0, 0.0)
        x, y, width, height = map(float, root.get("viewBox").split())
        points = []
        for record in selected:
            group = groups[record["svg_id"]]
            for path in group.findall(f"{{{SVG}}}path"):
                for a, b in _svg_path_lines(path.get("d", "")):
                    points.extend((a, b))
            for polygon in group.findall(f"{{{SVG}}}polygon"):
                points.extend(_polygon_points(polygon))
        left, top = min(x, min(point[0] + tx - 16 for point in points)), min(y, min(point[1] + ty - 16 for point in points))
        right, bottom = max(x + width, max(point[0] + tx + 16 for point in points)), max(y + height, max(point[1] + ty + 16 for point in points))
        root.set("viewBox", f"{left:g} {top:g} {right - left:g} {bottom - top:g}")
        root.set("width", f"{right - left:g}pt")
        root.set("height", f"{bottom - top:g}pt")


def _overview_legend(root, metadata):
    """Put a colour key outside the layout without adding graph entities."""
    bounds = [float(value) for value in root.get("viewBox", "0 0 360 80").split()]
    x, y, width, height = bounds
    routes = metadata["routes"]
    scheme_mode = metadata.get("candidate_unit") == "extra_ec_set"
    rows = len(routes) if routes else 1
    extra = 50 + rows * 24 + (24 if not metadata["nodes"] else 0)
    width = max(width, 360)
    root.set("viewBox", f"{x:g} {y:g} {width:g} {height + extra:g}")
    root.set("width", f"{width:g}pt")
    root.set("height", f"{height + extra:g}pt")
    legend = ET.SubElement(root, f"{{{SVG}}}g", {
        "id": "pathfiller-route-legend", "class": "route-legend",
        "transform": f"translate({x + 18:g} {y + height + 20:g})",
        "font-family": "Arial", "font-size": "13", "fill": _INK,
    })
    ET.SubElement(legend, f"{{{SVG}}}text", {"x": "0", "y": "0"}).text = "EC supplementation schemes" if scheme_mode else "Candidate routes"
    for index, route in enumerate(routes, 1):
        group = ET.SubElement(legend, f"{{{SVG}}}g", {
            "id": f"route-legend-{route['rank']}", "data-candidate-rank": str(route["rank"]),
            "data-route-color": route["color"],
        })
        if scheme_mode:
            group.set("data-extra-ecs", _json_text(route["extra_ecs"]))
            _set_title(group, f"Scheme {route['rank']}\n" + _supplement_label(route["extra_ecs"]))
        baseline = index * 24
        ET.SubElement(group, f"{{{SVG}}}line", {
            "x1": "0", "x2": "28", "y1": str(baseline - 4), "y2": str(baseline - 4),
            "stroke": route["color"], "stroke-width": "3",
        })
        steps = route["reaction_steps"]
        label = f"{'Scheme' if scheme_mode else 'Candidate'} {route['rank']}: {steps} steps" + (" (target only)" if steps == 0 else "")
        ET.SubElement(group, f"{{{SVG}}}text", {"x": "38", "y": str(baseline)}).text = label
    if not routes:
        ET.SubElement(legend, f"{{{SVG}}}text", {"x": "0", "y": "24"}).text = "No supplementation schemes" if scheme_mode else "No candidate routes"
    if not metadata["nodes"]:
        ET.SubElement(legend, f"{{{SVG}}}text", {"x": "0", "y": str((rows + 1) * 24)}).text = "No network compounds available"


def _set_title(group, text):
    title = group.find(f"{{{SVG}}}title")
    if title is None:
        title = ET.Element(f"{{{SVG}}}title")
        group.insert(0, title)
    title.text = text


def export_network_svg(graph: nx.MultiDiGraph, output_path, *, title: str | None = None,
                       overview: bool = False, highlighted_paths=None) -> Path:
    """Export all network nodes and edges as SVG without modifying the input graph.

    Graphviz's ``dot`` executable must be installed and available on PATH.
    Nodes use circular C labels and show their complete names below when
    Entity.name is available. Layout-only labels do not add business nodes or
    edges. Reactions show every EC annotation, and map
    links are dashed for legacy map-state networks. A graph marked with
    ``representation='compound'`` is flat: each C identifier is drawn once,
    with original map memberships and merged edge origins retained in metadata.
    The SVG includes ``metadata#pathfiller-network`` JSON and
    stable node/edge groups with ``data-*`` attributes for later applications.
    overview=True uses sfdp with compact labels and native hover titles;
    highlighted_paths maps ranks to ordered (source, target, merged key) edges.
    """
    source, metadata, nodes, edges, heading = _prepare_export(
        graph, title, overview=overview, highlighted_paths=highlighted_paths,
    )
    logger.debug("开始导出 SVG：%d 个节点、%d 条边，正在计算 Graphviz 布局", len(nodes), len(edges))
    engine = "sfdp" if overview else "dot"
    executable = shutil.which(engine)
    if executable is None:
        raise RuntimeError(f"SVG export requires Graphviz '{engine}' on PATH; install Graphviz first")
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
    logger.debug("SVG 导出完成：%s", output.resolve())
    return output
