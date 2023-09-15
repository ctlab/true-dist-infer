from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator

import pandas as pd

from tqdm import tqdm

n = 1000
k_max = 5000
tries = 200

flat_dir_est = FlatDirichletDBEstimator()
corr_dir_est = CorrectedDirichletDBEstimator(1)

columns = ['x', 'k_real', 'k_flat_dir_est', 'c_max']
rows = []

for try_number in tqdm(range(tries)):
    # print('try', try_number, 'of', tries)
    g = CyclicGenomeGraph(n, "gamma", [1])
    
    for k_real in tqdm(range(1, k_max)):
        g.do_k2_break()

        x = 2 * k_real / n
        k_flat_dir_est = flat_dir_est.predict_k(g, quiet=True)

        max_cycle_length = max(g.count_paths_and_cycles()[0].keys())

        rows.append([x, k_real, k_flat_dir_est, max_cycle_length])
    
# pd.DataFrame(rows, columns=columns).to_csv('tannier_simulations.csv', sep=',', index=False)
