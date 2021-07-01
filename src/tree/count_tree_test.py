from src.tree.fit_to_tree import count_tree_errors, count_correlation_pairs_alignment, count_alignment_pairs_alignment
from src.utils.parsers import parse_to_df
from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator


import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import chain
from itertools import combinations


tree_file = 'real_data/procars/tree.nwk'
infercars_file = 'real_data/procars/orthology_blocks.txt'
#
# errors1 = count_tree_errors(tree_file, infercars_file)
# errors2 = count_tree_errors(tree_file, infercars_file, error_func=lambda y, y_tree: math.sqrt(np.sum(np.power((y - y_tree), 2))))
# print(errors1)
# print(errors2)

class ParsimonyEstimator:
    def predict(self, g):
        return 0, g.d()

etimators = {
    'Parsimony': ParsimonyEstimator(),
    'Uniform': UniformDBEstimator(),
    'Dirichlet flat': FlatDirichletDBEstimator(),
    'Dirichlet 1/3': DirichletDBEstimator(1 / 3),
    # 'Corr dirichlet 1/3': CorrectedDirichletDBEstimator(1 / 3)
}

# res = 1000
# tree_file = 'real_data/Yersinia_pestis/tree.nwk'
# infercars_file = f'real_data/Yersinia_pestis/{res}/blocks_unique_coords.infercars'

df = parse_to_df(infercars_file)
species = df.species.unique()
sp_pairs = list(combinations(species, 2))

d_bs = []
# for _, (i, (sp1, sp2)) in zip(range(10), enumerate(sp_pairs)):
for i, (sp1, sp2) in enumerate(sp_pairs):
    print(f'{i} of {len(sp_pairs)}: {sp1}–{sp2}')
    g = RealDataGraph()
    g.infercars(df, sp1, sp2, cyclic=True)

    if g.b() != 0:
        d_bs.append(g.d() / g.b())

sns.set_style('whitegrid')
sns.histplot(d_bs, bins=np.arange(0.5, 1, 0.05))

plt.xlabel('d/b')

plt.tight_layout()
plt.savefig('y_procars_d_b_distribution.pdf')
plt.show()