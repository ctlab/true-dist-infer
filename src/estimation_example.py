from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator

g = LinearGenomeGraph(1000, 10, "gamma", [1])

uni_est = UniformDBEstimator()
flat_dir_est = FlatDirichletDBEstimator()
alpha1_3_dir_est = DirichletDBEstimator(1/3)

for i in range(300):
    g.do_k2_break()

print(uni_est.predict_k(g))
print(flat_dir_est.predict_k(g))
print(alpha1_3_dir_est.predict_k(g))

