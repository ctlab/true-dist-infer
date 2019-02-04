import networkx as nx
import random
import numpy as np
from operator import itemgetter
from src.graphs.common_utils import get_probabilities_by_distribution
from src.graphs.abstract_graph import AbstractGraph


class CyclicGenomeGraph(AbstractGraph):
    def __init__(self, n, distribution=None, params=[]):
        super().__init__()
        self._n = n
        ps = get_probabilities_by_distribution(distribution, params, size=n)
        for i in range(n):
            self.add_edge(str(i) + 'h', str((i + 1) % n) + 't', probability=ps[i]/np.sum(ps), label='red')

    def do_k2_break(self):
        self.do_k_break(2)

    def do_k_break(self, k):
        def generate_new_permutation(n):
            def check_permutation(p):
                return not any(
                    (p[i] % 2 == 0 and p[i] == p[i + 1] - 1)
                    or (p[i + 1] % 2 == 0 and p[i + 1] == p[i] - 1)
                    for i in range(0, len(p), 2)
                )

            permutation = np.random.permutation(n)
            while not check_permutation(permutation):
                permutation = np.random.permutation(n)
            return permutation

        def generate_new_probabilities(ps):
            def generate_normed_exp(n):
                xs = [random.expovariate(1) for _ in range(n)]
                return xs / np.sum(xs)

            r = np.array([generate_normed_exp(len(ps)) for _ in range(len(ps))])
            return np.dot(r.transpose(), ps)

        all_old_ps = list(map(lambda x: x[2]["probability"], self.edges.data()))
        edges_indexes = np.random.choice(self._n, k, replace=False, p=all_old_ps)

        edges = itemgetter(*edges_indexes)(list(self.edges))
        self.remove_edges_from(edges)

        old_ps = itemgetter(*edges_indexes)(all_old_ps)
        new_ps = generate_new_probabilities(old_ps)
        new_order = generate_new_permutation(k * 2)

        vertexes_in_edges = [
            item for sublist in map(lambda x: x[:2], edges) for item in sublist
        ]
        for i, new_w in zip(range(0, len(new_order), 2), new_ps):
            self.add_edge(
                vertexes_in_edges[new_order[i]],
                vertexes_in_edges[new_order[i + 1]],
                probability=new_w,
                label='red'
            )

    def construct_bp_graph(self):
        g = nx.MultiGraph(self)

        for i in range(self._n):
            g.add_edge(str(i) + 'h', str((i + 1) % self._n) + 't', label='black')

        return g