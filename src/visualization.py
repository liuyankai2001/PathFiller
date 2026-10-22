"""Display a KEGG compound network with NetworkX and Matplotlib."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from matplotlib.figure import Figure


_INK = "#303030"
_BLUE = "#1F497D"
_EDGE_COLORS = {True: "#26804A", False: "#B83D3D", None: _INK}


def draw_network(
    graph: nx.MultiDiGraph,
    *,
    title: str | None = None,
    figsize: tuple[float, float] = (16, 11),
    seed: int = 22,
) -> Figure:
    """Return a figure without modifying the graph or calling ``plt.show()``.

    Node data must contain ``entity`` and edge data ``connection``, as returned
    by ``src.kegg.parse_map``. Parallel edges and both directed edges of a
    reversible reaction are drawn separately. Positions are computed for this
    preview; they do not reproduce the original KEGG map layout.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    if not isinstance(graph, nx.MultiDiGraph):
        raise TypeError("graph must be a networkx.MultiDiGraph")
    # Keep NetworkX's lazily populated drawing/layout caches off the caller's
    # graph as well. Entity and Connection objects are only read.
    graph = graph.copy()

    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 13,
        "figure.dpi": 150,
        "svg.fonttype": "none",
    }):
        figure, ax = plt.subplots(figsize=figsize)
        figure.subplots_adjust(left=0.02, right=0.98, top=0.87, bottom=0.16)
        heading = title or graph.graph.get("name") or "Compound network"
        figure.suptitle(heading, fontsize=21, fontweight="bold", y=0.97)
        map_id = graph.graph.get("map_id", "KEGG")
        figure.text(
            0.5, 0.916,
            f"{map_id}  |  {graph.number_of_nodes()} compounds  |  "
            f"{graph.number_of_edges()} directed connections",
            ha="center", color=_INK, fontsize=13,
        )

        nodes = list(graph)
        # A simple undirected copy prevents parallel reactions from pulling
        # their compound labels on top of one another during layout.
        positions = nx.spring_layout(
            nx.Graph(graph), seed=seed, iterations=250, weight=None, k=0.85,
        )
        node_size = 3300
        tca_nodes = [node for node in nodes if graph.nodes[node]["entity"].is_tca]
        other_nodes = [node for node in nodes if not graph.nodes[node]["entity"].is_tca]

        # Curve parallel arrows separately. Reverse arrows curve on the opposite
        # side because their source and target have been exchanged.
        grouped_edges = defaultdict(list)
        for source, target, key, data in graph.edges(keys=True, data=True):
            grouped_edges[source, target].append((key, data["connection"]))
        for (source, target), edges in grouped_edges.items():
            has_reverse = source != target and (target, source) in grouped_edges
            for index, (key, connection) in enumerate(edges):
                radius = (
                    0.10 + 0.10 * index if has_reverse
                    else 0.14 * (index - (len(edges) - 1) / 2)
                )
                nx.draw_networkx_edges(
                    graph, positions, ax=ax, edgelist=[(source, target, key)],
                    node_size=node_size, nodelist=nodes,
                    edge_color=_EDGE_COLORS[connection.is_endogenous],
                    width=1.5, arrows=True, arrowstyle="-|>", arrowsize=17,
                    connectionstyle=f"arc3,rad={radius}",
                    min_source_margin=6, min_target_margin=6,
                )

        for selected, fill, text_color in (
            (tca_nodes, _BLUE, "white"), (other_nodes, "white", _INK),
        ):
            if not selected:
                continue
            nx.draw_networkx_nodes(
                graph, positions, nodelist=selected, node_size=node_size,
                node_color=fill, edgecolors=_BLUE, linewidths=1.5, ax=ax,
            )
            nx.draw_networkx_labels(
                graph, positions, labels={node: str(node) for node in selected},
                font_size=13, font_color=text_color, font_family="DejaVu Sans",
                ax=ax,
            )

        if not nodes:
            ax.text(0.5, 0.5, "No compounds to display", ha="center", va="center")
        ax.margins(0.12)
        ax.set_axis_off()
        legend = [
            Line2D([], [], color=color, lw=2, label=label)
            for color, label in (
                (_EDGE_COLORS[True], "Endogenous"),
                (_EDGE_COLORS[False], "No host EC match"),
                (_EDGE_COLORS[None], "Not assessed"),
            )
        ] + [
            Line2D([], [], marker="o", linestyle="none", markersize=12,
                   markerfacecolor=fill, markeredgecolor=_BLUE, label=label)
            for fill, label in (
                (_BLUE, "TCA stop compound"), ("white", "Other compound"),
            )
        ]
        figure.legend(
            handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.045),
            ncol=3, frameon=False, fontsize=13,
        )
        figure.text(
            0.5, 0.021,
            "Arrows follow KGML. Layout is computed for this preview.",
            ha="center", fontsize=13, color=_INK,
        )
        return figure


def main() -> None:
    """Preview TCA by default; use --no-show to export without opening a window."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", default="map00020", dest="map_id")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()
    if args.no_show:
        import matplotlib

        matplotlib.use("Agg")

    from .kegg import get_tca_graph, parse_map

    graph = get_tca_graph() if args.map_id == "map00020" else parse_map(args.map_id)
    figure = draw_network(graph)
    output_dir = Path(__file__).resolve().parent.parent / "figures" / "tca_network"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = "tca_network" if graph.graph["map_id"] == "map00020" else f"{graph.graph['map_id']}_network"
    for suffix in ("png", "svg"):
        output = output_dir / f"{stem}.{suffix}"
        figure.savefig(output, dpi=180, bbox_inches="tight", facecolor="white")
        print(output)

    import matplotlib.pyplot as plt

    if not args.no_show:
        plt.show()
    plt.close(figure)


if __name__ == "__main__":
    main()
