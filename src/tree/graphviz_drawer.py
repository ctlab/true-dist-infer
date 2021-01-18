from src.tree.fit_to_tree import fit_to_tree_graphviz

import numpy as np
from ete3 import Tree


tree_file = '/Users/alexey/PycharmProjects/true-dist-infer/src/tree/hg38.100way.tree'
res = '510k'
genomes = '7c'
infercars_file = f'/Users/alexey/PycharmProjects/true-dist-infer/real_data/ucsc_{genomes}_diff_res/{res}/Conserved.Segments'

etimators = [
    'Parsimony', 'Uniform', 'Dirichlet flat', 'Dirichlet 1/3'
]

est = etimators[0]

fit_to_tree_graphviz(tree_file, infercars_file, est, f'ucsc_{genomes}_{est}_{res}.pdf', f'ucsc_{genomes}_{est}_{res}.csv')