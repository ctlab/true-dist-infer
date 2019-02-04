from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph

filename = "file_%d.svg"
g = LinearGenomeGraph(5, 2)

for i in range(5):
    g.save_pygraphviz(filename % i)
    g.do_k2_break()
