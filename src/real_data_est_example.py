import scipy.stats

from src.graphs.real_data_graph import RealDataGraph
from src.estimators.flat_dirichlet_estimator import FlatDirichletDBEstimator
from src.estimators.uniform_db_estimator import UniformDBEstimator
from src.estimators.dirichlet_db_estimator import DirichletDBEstimator
from src.utils.parsers import parse_to_df
from src.utils.block_stats import dist_between_blocks

grimm_flag = False
file1 = "data/H.gen"
file2 = "data/M.gen"
file = "data/Conserved.Segments"
sp1 = "hg19"
sp2 = "mm10"
fit_alpha = True

norm_and_sort = lambda xs: sorted(list(filter(lambda x: x != 0, map(lambda x: x * len(xs) * 0.3 / sum(xs), xs))))


def get_dist_param(df, sp):
    dist = getattr(scipy.stats, "gamma")
    df_sp = df.loc[df['species'] == sp].sort_values(by=['chr', 'chr_beg'])
    ls = dist_between_blocks(df_sp, True)

    ls = list(sorted(filter(lambda l: l != 0, map(lambda l: l / sum(ls), ls))))
    # ls = norm_and_sort(ls)
    print(sorted(ls))
    print(sum(ls))
    param = dist.fit(ls, floc=0, fscale=1)
    print(param)


if grimm_flag:
    g = RealDataGraph()
    g.build_grimm(file1, file2)
else:
    df = parse_to_df(file)
    g = RealDataGraph()
    if fit_alpha:
        get_dist_param(df, sp1)
        get_dist_param(df, sp2)

    g.build_not_grimm(df, sp1, sp2)

dir_est = FlatDirichletDBEstimator()
uni_est = UniformDBEstimator()
gam_est = DirichletDBEstimator(0)

print("d", g.d())
print("b", g.b())
print(g.d(), g.b(), g.c(), len(g.edges))
print(uni_est.predict_k(g.d(), g.b()))
print(dir_est.predict_k(g.d(), g.b()))
print(gam_est.predict_k(g.d(), g.b()))
