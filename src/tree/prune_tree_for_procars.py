from src.utils.parsers import parse_to_df
from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator
from src.utils.parsers import parse_to_df


from ete3 import Tree, TextFace, TreeStyle
from itertools import combinations
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.optimize
import os.path
import pickle

tree_file = 'hg38.100way.tree'
species = '11'
infer_file = f'/Users/alexey/PycharmProjects/true-dist-infer/real_data/ucsc_{species}_diff_res/500k/Conserved.Segments'
# infer_file = f'/Users/alexey/PycharmProjects/true-dist-infer-model/real_data/procars/orthology_blocks.txt'

t = Tree(tree_file)

df = parse_to_df(infer_file)

df['len'] = df.chr_end - df.chr_beg
print(df.describe())

t.prune(df.species.unique())

out_tree_file = f'for_procars_{species}.tree'
t.write(outfile=out_tree_file, format=9)
