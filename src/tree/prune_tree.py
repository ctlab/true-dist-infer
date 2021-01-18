from src.utils.parsers import parse_to_df
from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator

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

t = Tree(tree_file)

ucsc_species = {'hg38': ('Human', '2013'), #0
                'papAnu4': ('Baboon', '2017'),
                'felCat9': ('Cat', '2017'), #1
                'galGal6': ('Chicken', '2018'),
                'panTro6': ('Chimp', '2018'),
                'bosTau9': ('Cow', '2018'),
                'macFas5': ('Crab-eating macaque', '2013'),
                'canFam4': ('Dog', '2020'), #0
                'gorGor6': ('Gorilla', '2019'),
                'chlSab2': ('Green Monkey', '2014'),
                'equCab3': ('Horse', '2018'), #1
                'calJac3': ('Marmoset', '2010'),
                'mm10': ('Mouse', '2012'), #0
                'monDom5': ('Opossum', '2007'),
                'ponAbe3': ('Orangutan', '2018'),
                'susScr11': ('Pig', '2017'),
                'ornAna2': ('Platypus', '2007'),
                'oryCun2': ('Rabbit', '2009'),
                'rn6': ('Rat', '2014'), #0
                'rheMac10': ('Rhesus', '2019'), #1
                'oviAri4': ('Sheep', '2015'),
                'danRer10': ('Zebrafish', '2014')}

t.prune(ucsc_species.keys())

for node in t.traverse():
    # Hide node circles
    node.img_style['size'] = 0

    if node.is_leaf():
        name_face = TextFace(ucsc_species[node.name][0] + ', ' + ucsc_species[node.name][1])
        node.add_face(name_face, column=0)


ts = TreeStyle()
ts.mode = 'r'
ts.scale = 500
# Disable the default tip names config
ts.show_leaf_name = False
# ts.show_branch_support = False
t.render('pruned_tree.pdf', w=1000, tree_style=ts)