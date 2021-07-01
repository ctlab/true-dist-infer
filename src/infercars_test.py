import argparse

from src.graphs.real_data_graph import RealDataGraph

from src.utils.parsers import parse_to_df
from src.real_data_est_common import get_dist_param, print_graph_stats



g = RealDataGraph()

# file = '/Users/alexey/Downloads/data/Conserved.Segments'
# file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/WORKING 50k.test/Conserved.Segments'
# file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/WORKING 50k.maml/Conserved.Segments'
# file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/WORKING 50k.hg38.mm10_not_alt/Conserved.Segments'
file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/WORKING 50k.maml.ucsc/Conserved.Segments'
# file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/TEST.50K/SFs/Conserved.Segments'
# file = '/Users/alexey/PycharmProjects/true-dist-infer-model/real_data/MGRA_mammals/mammals.infercars'
sp1 = 'hg38'
sp2 = 'mm10'

df = parse_to_df(file)

alpha = 0.33

g.infercars(df, sp1, sp2)

print_graph_stats(g, alpha, True)
