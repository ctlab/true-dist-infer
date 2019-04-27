from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph

filename = "example.svg"
g = LinearGenomeGraph(5, 1)

for i in range(5):
    g.do_k2_break()

g.save_pygraphviz(filename)
