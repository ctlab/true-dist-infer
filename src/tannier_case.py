import argparse

from src.graphs.real_data_graph import RealDataGraph
from src.real_data_est_common import print_graph_stats

g = RealDataGraph()

file1, file2 = '/Users/alexey/Downloads/conica_chr_X.gen', '/Users/alexey/Downloads/conica_chr_Y.gen'
alpha = 0.33
g.build_grimm(file1, file2)

print_graph_stats(g, alpha, True)

g.save_pygraphviz("example.pdf")