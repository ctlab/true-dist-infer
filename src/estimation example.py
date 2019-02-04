from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator

filename = "file_%d.svg"
g = LinearGenomeGraph(1000, 10, "gamma", [1])

dir_est = DirichletDBEstimator()
uni_est = UniformDBEstimator()
gam_est = DirichletDBEstimator(1)

for i in range(1000):
    print(i)
    g.do_k2_break()

print(dir_est.predict_k(g.d(), g.b()))
print(uni_est.predict_k(g.d(), g.b()))
print(gam_est.predict_k(g.d(), g.b()))

