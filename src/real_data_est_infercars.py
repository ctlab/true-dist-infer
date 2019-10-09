import argparse

from src.graphs.real_data_graph import RealDataGraph

from src.utils.parsers import parse_to_df
from src.real_data_est_common import get_dist_param, print_graph_stats

parser = argparse.ArgumentParser(
    description='Construct a breakpoint graph with real data provided in infercars format. '
                'Then calculates the necessary statistics and estimate true evolutionary distance using different approaches.',
)

parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('--file', '-f', required=True,
                      help='Path to file in infercars format')

required.add_argument('--species1', '-s1', required=True,
                      help='Name of first species')

required.add_argument('--species2', '-s2', required=True,
                      help='Name of second species')

optional.add_argument('--fit_alpha', '-fa', action='store_true',
                      help='Fit alpha')

optional.add_argument('--alpha', '-a', type=float, default=0.33,
                      help='Set alpha to a specific value. Float value, default is 0.33')

optional.add_argument('--stats', '-s', action='store_true',
                      help='Show main statistics of the graph: b, d, p_even, p_odd and paths and cycles count')

args = parser.parse_args()
d = vars(args)

g = RealDataGraph()

file, sp1, sp2 = d['file'], d['species1'], d['species2']
df = parse_to_df(file)

alpha = d['alpha']
if d['fit_alpha']:
    alpha = get_dist_param(df, sp1, sp2)
    print(f'Alpha fitted, alpha={alpha}')

g.infercars(df, sp1, sp2)

print_graph_stats(g, alpha, d['stats'])
