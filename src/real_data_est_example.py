import scipy.stats

from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator, CorrectedDirichletDBEstimator
from src.utils.parsers import parse_to_df
from src.utils.block_stats import dist_between_blocks


grimm_flag = False
file1 = "data/H.gen"
file2 = "data/M.gen"
file = "data/Conserved.Segments"
sp1 = "hg19"
sp2 = "mm10"
fit_alpha = True
alpha = 1/3


def get_dist_param(df, sp1, sp2):
    def df_sp(sp):
        return df.loc[df['species'] == sp].sort_values(by=['chr', 'chr_beg'])

    dist = getattr(scipy.stats, "gamma")
    df_sp1 = df_sp(sp1)
    df_sp2 = df_sp(sp2)

    ls = dist_between_blocks(df_sp1) + dist_between_blocks(df_sp2)
    param = dist.fit(ls)[0]
    return param


if grimm_flag:
    g = RealDataGraph()
    g.build_grimm(file1, file2)
else:
    df = parse_to_df(file)
    g = RealDataGraph()
    if fit_alpha:
        alpha = get_dist_param(df, sp1, sp2)
    print(alpha)
    g.infercars(df, sp1, sp2)

flat_dir_est = FlatDirichletDBEstimator()
uni_est = UniformDBEstimator()
alpha1_3_dir_est = DirichletDBEstimator(1/3)
corr_alpha1_3_dir_est = CorrectedDirichletDBEstimator(1/3)
print((g.p_even() + g.p_odd()) / 2)

print("d", g.d())
print("b", g.b())
print(g.d(), g.b(), g.c(), len(g.edges))
print(uni_est.predict_k(g))
print(flat_dir_est.predict_k(g))
print(alpha1_3_dir_est.predict_k(g))
print(corr_alpha1_3_dir_est.predict_k(g))
