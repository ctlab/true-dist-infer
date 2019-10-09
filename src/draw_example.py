from src.graphs.linear_genome_graph import LinearGenomeGraph

g = LinearGenomeGraph(5, 1)

for i in range(5):
    g.do_k2_break()

g.save_pygraphviz("example.svg")
