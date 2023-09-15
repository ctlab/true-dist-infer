import random
import networkx as nx
import numpy as np

from src.graphs.common_utils import get_probabilities_by_distribution
from src.graphs.abstract_graph import AbstractGraph


class LinearGenomeGraph(AbstractGraph):
    def __init__(self, n, chrs, distribution=None, params=[]):
        super().__init__()
        self._n, self.chrs = n, chrs
        self.ns = list(range(n))
        self.add_nodes_from([str(i) + 't' for i in self.ns])
        self.add_nodes_from([str(i) + 'h' for i in self.ns])
        [self.add_edge(str(i) + 't', str(i) + 'h', label='block') for i in self.ns]

        ps = get_probabilities_by_distribution(distribution, params, size=(n + chrs))
        i = 0
        for xs in np.array_split(self.ns, chrs):
            self.add_edge(str(xs[0]) + 't', 'telomere', label='adj-red', probability=ps[i] / sum(ps))
            self.add_edge(str(xs[-1]) + 'h', 'telomere', label='adj-red', probability=ps[i + 1] / sum(ps))
            i += 2
            for x, y in zip(xs[:-1], xs[1:]):
                self.add_edge(str(x) + 'h', str(y) + 't', label='adj-red', probability=ps[i] / sum(ps))
                i += 1

    def construct_bp_graph(self):
        g = nx.MultiGraph()

        g.add_nodes_from(map(lambda x: str(x) + 't', self.ns))
        g.add_nodes_from(map(lambda x: str(x) + 'h', self.ns))

        for xs in np.array_split(self.ns, self.chrs):
            for x, y in zip(xs[:-1], xs[1:]):
                g.add_edge(str(x) + 'h', str(y) + 't', label='black')

        # adding bp edges
        for e in self.label_edges_with_probability('adj-red'):
            if e[1] != 'telomere':
                g.add_edge(e[0], e[1], label='red')

        return g

    def label_edges_with_probability(self, label):
        # return list(self.edges(data=True))
        return list(map(lambda e: (e[0], e[1], e[2]),
                        filter(lambda e: e[2]['label'] == label, self.edges(data=True))))

    def do_k2_break(self):
        def insert_edge(e, label='adj-red'):
            self.add_edge(e[0], e[1], probability=e[2], label=label)

        def delete_edge(e):
            self.remove_edge(e[0], e[1])

        def new_probabilities(w1, w2=0):
            r1, r2 = random.random(), random.random()
            return w1 * r1 + w2 * r2, w1 * (1 - r1) + w2 * (1 - r2)

        break_is_done = False
        while not break_is_done:
            es = self.label_edges_with_probability('adj-red')
            all_ps = list(map(lambda e: e[2]['probability'], es))
            i1, i2 = np.random.choice(len(es), size=2, replace=True, p=all_ps)
            e1, e2 = es[i1], es[i2]

            if e1[0] == 'telomere' or e2[0] == 'telomere':
                print("WARNING")

            # Breaking chromosome
            if i1 == i2:
                if e1[1] != 'telomere':
                    # print("________Breaking_chromosome")
                    delete_edge(e1)
                    p1, p2 = new_probabilities(e1[2]['probability'])
                    insert_edge((e1[0], 'telomere', p1))
                    insert_edge((e1[1], 'telomere', p2))
                    break
                else:
                    continue

            # Merging chromosome
            if e1[1] == 'telomere' and e2[1] == 'telomere':
                # print("________Merging_chromosome")
                delete_edge(e1)
                delete_edge(e2)
                insert_edge((e1[0], e2[0], e1[2]['probability'] + e2[2]['probability']), label='adj-red')

                if len(nx.cycle_basis(self.without_telomere())) > 0:
                    delete_edge((e1[0], e2[0]))
                    insert_edge((e1[0], e1[1], e1[2]['probability']))
                    insert_edge((e2[0], e2[1], e2[2]['probability']))
                    continue
                break

            # Translocation
            if i1 == i2:
                continue

            e1 = (e1[0], e1[1], e1[2]['probability'])
            e2 = (e2[0], e2[1], e2[2]['probability'])
            p1, p2 = new_probabilities(e1[2], e2[2])

            delete_edge(e1)
            delete_edge(e2)

            if random.random() > 0.5:
                ne1 = (e1[0], e2[1], p1)
                ne2 = (e1[1], e2[0], p2)
            else:
                ne1 = (e1[0], e2[0], p1)
                ne2 = (e1[1], e2[1], p2)

            if self.has_edge(ne1[0], ne1[1]) or self.has_edge(ne2[0], ne2[1]):
                insert_edge(e1)
                insert_edge(e2)
                continue

            insert_edge(ne1)
            insert_edge(ne2)

            if len(nx.cycle_basis(self.without_telomere())) > 0:
                delete_edge(ne1)
                delete_edge(ne2)
                insert_edge(e1)
                insert_edge(e2)
                continue

            break_is_done = True

    def without_telomere(self):
        g = nx.Graph(self)
        g.remove_node('telomere')
        return g
