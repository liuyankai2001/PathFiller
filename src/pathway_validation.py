"""单条化合物路线内部的依赖检查，不判断全部共底物的外部供应。"""

from .schemas import PathWay


def analyze_pathway_dependencies(pathway: PathWay) -> dict:
    """起点和较早反应的全部产物可用；提前使用后续化合物构成内部循环。

    其他未知共底物只列出，本轮不判断其供应。结果无内部循环也不等于路线已生物学验证。
    """
    if pathway is None:
        raise ValueError("没有候选路线，无法检查路线依赖")
    connections = pathway.connections
    available = {connections[0].source.id} if connections else set()
    future = {connection.target.id for connection in connections}
    visited = set(available)
    cycles, background = [], []
    for index, connection in enumerate(connections, 1):
        required = set(connection.substrate_ids) - available
        cyclic = required & future
        if cyclic:
            cycles.append({"kind": "downstream_substrate", "step": index,
                           "compound_ids": sorted(cyclic), "reaction_ids": list(connection.reaction_ids)})
        if connection.target.id in visited:
            cycles.append({"kind": "compound_cycle", "step": index,
                           "compound_ids": [connection.target.id], "reaction_ids": list(connection.reaction_ids)})
        other = required - future
        if other:
            background.append({"step": index, "compound_ids": sorted(other)})
        available.update(connection.product_ids)
        available.add(connection.target.id)
        visited.add(connection.target.id)
        future.discard(connection.target.id)
    return {
        "status": "circular_dependency" if cycles else "no_internal_cycle_detected",
        "circular_dependencies": cycles, "unverified_substrates": background,
        "scope": "仅检查单条路线内部的循环依赖；未检查其他共底物供应、计量、碳流和宿主底物适用性。",
    }
