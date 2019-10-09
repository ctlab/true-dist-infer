import argparse

from src.graphs.real_data_graph import RealDataGraph
from src.real_data_est_common import print_graph_stats

parser = argparse.ArgumentParser(
    description='Construct a breakpoint graph with real data provided in grimm format. '
                'Then calculates the necessary statistics and estimate true evolutionary distance using different approaches.',
)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('--file1', '-f1', required=True,
                      help='Path to first file in grimm format')

required.add_argument('--file2', '-f2', required=True,
                      help='Path to second file in grimm format')

optional.add_argument('--alpha', '-a', type=float, default=0.33,
                      help='Set alpha to a specific value. Float, default value is 0.33')

optional.add_argument('--stats', '-s', action='store_true',
                      help='Show main statistics of the graph: b, d, p_even, p_odd and paths and cycles count')

args = parser.parse_args()
d = vars(args)

g = RealDataGraph()

file1, file2 = d['file1'], d['file2']
alpha = d['alpha']
g.build_grimm(file1, file2)

print_graph_stats(g, alpha, d['stats'])
