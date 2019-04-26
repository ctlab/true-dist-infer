import math
import networkx as nx
import pygraphviz as pgv

ind_by_node = lambda n: int(n[:-1]) * 2 + (1 if n[-1] == 'h' else 0)


class AbstractGraph(nx.MultiGraph):
    def __init__(self):
        super().__init__()

    def construct_bp_graph(self):
        raise NotImplementedError('subclasses must override construct_bp_graph!')

    def count_components(self, predicate):
        return len(list(filter(predicate, nx.connected_component_subgraphs(self.construct_bp_graph()))))

    def p_even(self):
        return self.count_components(lambda c: len(c) != len(c.edges) and len(c.edges) % 2 == 0)

    def p_odd(self):
        return self.count_components(lambda c: len(c) != len(c.edges) and len(c.edges) % 2 == 1)

    def chr(self):
        return (self.p_even() + self.p_odd()) / 2

    def p_m(self, m):
        return self.count_components(lambda c: len(c) != len(c.edges) and len(c.edges) == m)

    def c(self):
        return self.count_components(lambda c: len(c) == len(c.edges))

    def c_m(self, m):
        return self.count_components(lambda c: len(c) == len(c.edges) and len(c) == m * 2)

    def n(self):
        return len(self) // 2

    def d(self):
        return self.n() - self.c() - self.p_even() // 2

    def b(self):
        b = 0.0
        for component in nx.connected_component_subgraphs(self.construct_bp_graph()):
            if len(component) > 2 or (len(component) == 2 and len(component.edges) == 1):
                b += len(component) / 2
        return b

    def save_pygraphviz(self, filename, prog='neato'):
        def cacl_pos(node):
            radius = self.n() ** (6 / 7) / 2
            ind = ind_by_node(node)
            return "%f,%f!" % (math.cos(math.pi / self.n() * ind) * radius,
                               math.sin(math.pi / self.n() * ind) * radius)

        bpg = self.construct_bp_graph()
        a = pgv.AGraph(strict=False)

        for node in bpg.nodes:
            a.add_node(node, shape="circle", pos=cacl_pos(node))

        for edge in bpg.edges(data='label'):
            a.add_edge(edge[0], edge[1], color=edge[2])

        a.draw(filename, prog=prog)
