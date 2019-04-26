import scipy.stats
import sys

from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator
from src.utils.parsers import parse_to_df
from src.utils.block_stats import dist_between_blocks

# grimm_flag = True
# file1 = "data/H.gen"
# file2 = "data/M.gen"
# file = "data/Conserved.Segments"
# sp1 = "hg19"
# sp2 = "mm10"
# fit_alpha = True
alpha = 1 / 3


def get_dist_param(df, sp1, sp2):
    def df_sp(sp):
        return df.loc[df['species'] == sp].sort_values(by=['chr', 'chr_beg'])

    dist = getattr(scipy.stats, "gamma")
    df_sp1 = df_sp(sp1)
    df_sp2 = df_sp(sp2)

    ls = dist_between_blocks(df_sp1) + dist_between_blocks(df_sp2)
    param = dist.fit(ls)[0]
    return param


g = RealDataGraph()
if sys.argv[1] == 'grimm':
    file1, file2 = sys.argv[2], sys.argv[3]
    g.build_grimm(file1, file2)
    if len(sys.argv) > 4:
        alpha = float(sys.argv[4])
else:
    if sys.argv[1] == 'infercars':
        file, sp1, sp2 = sys.argv[2], sys.argv[3], sys.argv[4]
        df = parse_to_df(file)
        g = RealDataGraph()
        if len(sys.argv) > 5:
            if sys.argv[5] == 'fit':
                alpha = get_dist_param(df, sp1, sp2)
                print(f'Alpha fitted, alpha={alpha}')
            else:
                alpha = float(sys.argv[5])
        g.infercars(df, sp1, sp2)
    else:
        print('Undefined format')
        exit(0)

flat_dir_est = FlatDirichletDBEstimator()
uni_est = UniformDBEstimator()
alpha1_3_dir_est = DirichletDBEstimator(alpha)
corr_alpha1_3_dir_est = CorrectedDirichletDBEstimator(alpha)

print(f"d={g.d()}, b={g.b()}")
print(f"Uniform estimator, k: {uni_est.predict_k(g)}")
print(f"Flat Dirichlet estimator, k: {flat_dir_est.predict_k(g)}")
print(f"Dirichlet estimator with alpha={alpha}, k: {alpha1_3_dir_est.predict_k(g)}")
print(f"Corrected dirichlet estimator with alpha={alpha}, k: {corr_alpha1_3_dir_est.predict_k(g)}")
