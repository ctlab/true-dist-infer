from src.utils.parsers import parse_to_df
from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator
from src.real_data_est_common import print_graph_stats

from ete3 import Tree
from itertools import combinations
from collections import defaultdict
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from itertools import chain

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.optimize
import os.path
import pickle


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


def construct_A(sp_pairs, t):
    A = np.zeros((len(sp_pairs), len(t.get_edges())))
    for i, (sp1, sp2) in enumerate(sp_pairs):
        n1 = t & sp1
        n2 = t & sp2
        for j, (v, u) in enumerate(t.iter_edges()):
            if (n1 in v and n2 in u) or (n2 in v and n1 in u):
                A[i, j] = 1
    return A


def count_distances_bp(sp_pairs, df, etimators, cache_file, cyclic):
    if os.path.isfile(cache_file):
        with open(cache_file, 'rb') as f:
            print('reading distances from cache')
            distances = pickle.load(f)
    else:
        distances = {}
        for i, (sp1, sp2) in enumerate(sp_pairs):
            print(f'{i} of {len(sp_pairs)}: {sp1}–{sp2}')
            g = RealDataGraph()
            g.infercars(df, sp1, sp2, cyclic=cyclic)

            for est_name, est in etimators.items():
                distances[(sp1, sp2, est_name)] = est.predict(g)[1]
        with open(cache_file, 'wb') as f:
            print('writing distances to cache')
            pickle.dump(distances, f)

    return distances


def construct_y(sp_pairs, distances, est):
    return np.array([distances[(sp1, sp2, est)] for sp1, sp2 in sp_pairs])


def filter_alt_(df):
    allowed_blocks = set()

    for block, df_block in df.groupby('block'):
        if all(['alt' not in chr for chr in df_block['chr']]):
            allowed_blocks.add(block)

    return df.loc[df['block'].isin(allowed_blocks)].copy()


def __count_species_distances(infercars_file, cyclic, filter_alt=True):
    df = parse_to_df(infercars_file)
    if filter_alt: df = filter_alt_(df)
    species = df.species.unique()
    sp_pairs = list(combinations(species, 2))

    cache_file = infercars_file + '.cache.pkl'
    distances = count_distances_bp(sp_pairs, df, etimators, cache_file, cyclic=cyclic)

    return distances, sp_pairs


def __fit_to_tree(sp_pairs, distances, est_name, tree_file, unroot):
    t = Tree(tree_file)

    species = [sp for pair in sp_pairs for sp in pair]
    t.prune(species)

    if unroot: t.unroot()
    # print('Tree:', t.write(format=5))

    # y = Aw, finding w
    A = construct_A(sp_pairs, t)
    y = construct_y(sp_pairs, distances, est_name)

    w, error = scipy.optimize.nnls(A, y)
    y_tree = A.dot(w)

    return t, w, y, y_tree


def fit_to_tree_graphviz(tree_file, infercars_file, est_name, graphviz_file, table_file, filter_alt=True, round=2,
                         unroot=True, cyclic=False):
    import pygraphviz as pgv

    def rec_convert_t_to_graphviz(n, g, edge_to_w, acc_name):
        g.add_node(acc_name, label=(n.name if n.is_leaf() else ''), shape=('ellipse' if n.is_leaf() else 'point'))
        for child, direction in zip(n.children, ['->0', '->1', '->2']):
            child_name = acc_name + direction
            w = edge_to_w[frozenset(child.get_leaves())]
            g.add_edge(acc_name, child_name, label=w.round(round))
            rec_convert_t_to_graphviz(child, g, edge_to_w, child_name)

    distances, sp_pairs = __count_species_distances(infercars_file, cyclic, filter_alt)
    t, ws, y, y_tree = __fit_to_tree(sp_pairs, distances, est_name, tree_file, unroot)

    edge_to_w = {frozenset(e[0]): w for e, w in zip(t.iter_edges(), ws)}
    G = pgv.AGraph(strict=False, directed=True)
    rec_convert_t_to_graphviz(t, G, edge_to_w, 'root')

    G.draw(graphviz_file, prog='dot')
    with open(table_file, 'w') as f:
        print('sp1,sp2,y,y_tree', file=f)
        for (sp1, sp2), y_, y_tree_ in zip(sp_pairs, y, y_tree):
            print(sp1, sp2, y_, y_tree_.round(round), sep=',', file=f)


def count_tree_errors(tree_file, infercars_file, dh=None, alignment_type=None, filter_alt=True,
                      error_func=lambda y, y_tree: np.sqrt(np.sum(np.power((y - y_tree) / y, 2))), unroot=True,
                      cyclic=False):
    distances, sp_pairs = __count_species_distances(infercars_file, cyclic, filter_alt)

    if alignment_type != None:
        for sp1, sp2 in sp_pairs:
            distances[(sp1, sp2, alignment_type)] = dh.get_dist(sp1, sp2)

    errors = {}
    ls = etimators if alignment_type == None else chain(etimators, [alignment_type])
    for est in ls:
        _1, _2, y, y_tree = __fit_to_tree(sp_pairs, distances, est, tree_file, unroot)
        errors[est] = error_func(y, y_tree)

    return errors


def count_alignment_pairs_alignment(infercars_file, distance_holder, filter_alt=True, cyclic=False):
    distances, sp_pairs = __count_species_distances(infercars_file, cyclic, filter_alt)

    al_ds = [distance_holder.get_dist(sp1, sp2) for sp1, sp2 in sp_pairs]
    est_ds = {est: [distances[(sp1, sp2, est)] for sp1, sp2 in sp_pairs]for est in etimators}

    return est_ds, al_ds


def count_correlation_score(xs, ys, fit_intercept=True):
    xs = np.array(xs)
    xs = xs.reshape(-1, 1)
    reg = LinearRegression(fit_intercept=fit_intercept).fit(xs, ys)
    return reg.score(xs, ys)


def count_correlation_pairs_alignment(infercars_file, distance_holder, filter_alt=True, cyclic=False):
    distances, sp_pairs = __count_species_distances(infercars_file, cyclic, filter_alt)

    tree_ds = [distance_holder.get_dist(sp1, sp2) for sp1, sp2 in sp_pairs]
    corrs = {}
    for est in etimators:
        ds = [distances[(sp1, sp2, est)] for sp1, sp2 in sp_pairs]
        corrs[est] = pearsonr(tree_ds, ds)[0]
        # corrs[est] = count_correlation_score(tree_ds, ds)

    return corrs