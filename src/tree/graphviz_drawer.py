from src.tree.fit_to_tree import fit_to_tree_graphviz

import numpy as np
from ete3 import Tree


# tree_file = '/Users/alexey/PycharmProjects/true-dist-infer/src/tree/hg38.100way.tree'
# res = '510k'
# genomes = '7c'
# infercars_file = f'/Users/alexey/PycharmProjects/true-dist-infer/real_data/ucsc_{genomes}_diff_res/{res}/Conserved.Segments'

res = 3000
suffix = '_5_left'
tree_file = f'real_data/Yersinia_pestis/tree{suffix}.nwk'
infercars_file = f'real_data/Yersinia_pestis/{res}/blocks_unique_coords{suffix}.infercars'

etimators = [
    'Parsimony', 'Uniform', 'Dirichlet flat', 'Dirichlet 1/3'
]

est = etimators[0]

fit_to_tree_graphviz(tree_file, infercars_file, est,
                     f'yersinia_tree_{est.replace("/", "_")}_{res}{suffix}.pdf',
                     f'yersinia_{est.replace("/", "_")}_{res}{suffix}.csv')