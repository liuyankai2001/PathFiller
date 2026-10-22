from src.search import build_pathway_network

result = build_pathway_network("C05432")
graph = result.graph  # 合并后的 NetworkX MultiDiGraph

print("已加载的 map：", list(result.graphs))
print("化合物数：", graph.number_of_nodes())
print("反应边数：", graph.number_of_edges())
print("数据是否完整：", result.is_complete)