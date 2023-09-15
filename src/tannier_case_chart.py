from src.graphs.cyclic_genome_graph import CyclicGenomeGraph
from src.graphs.linear_genome_graph import LinearGenomeGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('tannier_simulations.csv')
# print(df)

sns.set_style('whitegrid')
# sns.lineplot(data=df, x='x', y='k_real / k_uni_est')
sns.lineplot(data=df, x='x', y='c_max')
# sns.lineplot(data=df, x='x', y='k_flat_dir_est')
# sns.lineplot(data=df, x='x', y='k_real / k_alpha_1_3_dir_est')

plt.tight_layout()
plt.savefig('threshold_simulations_lineplot_c_max.pdf')
plt.show()