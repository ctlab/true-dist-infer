import scipy.stats

from src.utils.block_stats import dist_between_blocks
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator


def get_dist_param(df, sp1, sp2):
    def df_sp(sp):
        return df.loc[df['species'] == sp].sort_values(by=['chr', 'chr_beg'])

    dist = getattr(scipy.stats, "gamma")
    df_sp1 = df_sp(sp1)
    df_sp2 = df_sp(sp2)

    ls = dist_between_blocks(df_sp1) + dist_between_blocks(df_sp2)
    param = dist.fit(ls)[0]
    return param


def print_graph_stats(g, alpha, stats):
    flat_dir_est = FlatDirichletDBEstimator()
    uni_est = UniformDBEstimator()
    alpha1_3_dir_est = DirichletDBEstimator(alpha)
    corr_alpha1_3_dir_est = CorrectedDirichletDBEstimator(alpha)

    print(f"Uniform estimator, k: {uni_est.predict_k(g)}")
    print(f"Flat Dirichlet estimator, k: {flat_dir_est.predict_k(g)}")
    print(f"Dirichlet estimator with alpha={alpha}, k: {alpha1_3_dir_est.predict_k(g)}")
    print(f"Corrected dirichlet estimator with alpha={alpha}, k: {corr_alpha1_3_dir_est.predict_k(g)}")

    if stats:
        print()
        print(f'd = {g.d()}', f'b = {g.b()}', f'p_even = {g.p_even()}', f'p_odd = {g.p_odd()}', f'c = {g.c()}', sep='\n')
        cycles, paths = g.count_paths_and_cycles()

        print('paths:')
        for key in sorted(paths.keys()):
            print(f'  s_{key} = {paths[key]}')

        print('cycles:')
        for key in sorted(cycles.keys()):
            print(f'  c_{key} = {cycles[key]}')