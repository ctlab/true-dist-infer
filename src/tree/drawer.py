import math

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from itertools import combinations
from collections import defaultdict
from sklearn import linear_model
from sklearn.metrics import r2_score
from scipy.stats import pearsonr, spearmanr

from src.tree.fit_to_tree import count_tree_errors, count_correlation_pairs_alignment, count_alignment_pairs_alignment
from src.tree.distance_holder import DistanceHolder

tree_file = '/Users/alexey/PycharmProjects/true-dist-infer/src/tree/hg38.100way.tree'


def diff_ks_ucsc_7():
    # infercars_file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/maml.ucsc_7/%dk/Conserved.Segments'
    # infercars_file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/chicken/%dk/Conserved.Segments'
    genomes = "7"
    infercars_file = f'real_data/ucsc_{genomes}_diff_res/%dk/Conserved.Segments'

    alignment_type = 'nuc_alpha_05'
    distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'
    dh = DistanceHolder(distances_file)

    ys = defaultdict(list)
    ks = list(range(100, 1010, 500))
    # ks = list(range(50, 1001, 10))

    for k in ks:
        print(k)
        # errors = count_tree_errors(tree_file, infercars_file % k, dh, alignment_type, error_func=lambda y, y_tree: math.sqrt(np.sum(np.power((y - y_tree), 2))))
        errors = count_tree_errors(tree_file, infercars_file % k, dh, alignment_type)
        # print(k, errors[alignment_type] > errors['Parsimony'], errors[alignment_type] - errors['Parsimony'])
        print(k, errors['Uniform'] > 0.15, errors['Uniform'])
        for est, error in errors.items():
            ys[est].append(error)

    sns.set(style="whitegrid", font="serif")

    for est, ys in ys.items():
        if est == alignment_type: continue
        plt.plot(ks, ys, label=f'{est}', ls=('--' if est == alignment_type else '-'))

    # plt.title(f'Tree error')
    plt.xlabel('minimal rearrangement size, kilobase')
    plt.ylabel('tree error')
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(f'ucsc_{genomes}_ks_tree_rel.png')
    # plt.savefig(f'ucsc_{genomes}_ks_tree.png')
    plt.show()


def diff_ks_ucsc_7_correlation_pairs():
    genomes = "7"
    infercars_file = f'real_data/ucsc_{genomes}_diff_res/%dk/Conserved.Segments'
    alignment_type = 'nuc_alpha_05'
    distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'

    dh = DistanceHolder(distances_file)
    ys = defaultdict(list)
    # ks = list(range(50, 2001, 300))
    ks = list(range(50, 1001, 10))

    for k in ks:
        print(k)
        corrs = count_correlation_pairs_alignment(infercars_file % k, dh)
        for est, corr in corrs.items():
            ys[est].append(corr)

    sns.set(style="whitegrid", font="serif")

    for est, ys in ys.items():
        plt.plot(ks, ys, label=f'{est}')

    # plt.title(f'Tree corr')
    plt.xlabel('minimal rearrangement size, kilobase')
    plt.ylabel('tree corr for pairs')
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(f'ucsc_{genomes}_ks_tree_corr_{alignment_type}.png')
    plt.show()


# def scatter_correlation_pairs():
#     # genomes = "7c"
#     # res = 500
#     # infercars_file = f'real_data/ucsc_{genomes}_diff_res/{res}k/Conserved.Segments'
#     # alignment_type = 'nuc_alpha_05'
#     # distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'
#
#     infercars_file = f'real_data/procars/orthology_blocks.txt'
#     alignment_type = 'nuc_05'
#     distances_file = f'real_data/align_distances/distances_procars_{alignment_type}.csv'
#
#     dh = DistanceHolder(distances_file)
#
#     ests_ds, al_ds = count_alignment_pairs_alignment(infercars_file, dh)
#
#     for est, est_ds in ests_ds.items():
#         plt.figure()
#         sns.set(style="whitegrid", font="serif")
#         sns.regplot(x=est_ds, y=al_ds, label=f'{est}', scatter_kws={'s': 3}, line_kws={'lw': 0.5})
#
#         # plt.annotate("$R^2$ = {:.2f}, pearson = {:.2f}".format(r2_score(est_ds,al_ds), pearsonr(est_ds, al_ds)[0]), (50, 0.75))
#         plt.annotate("$R^2$ = {:.2f}, pearson = {:.2f}".format(r2_score(np.array(al_ds) / max(al_ds), np.array(est_ds) / max(est_ds)), pearsonr(est_ds, al_ds)[0]), (50, 0.75))
#
#         plt.xlabel('Estimated k')
#         plt.ylabel('Distance in alignment')
#         # plt.legend(loc='best')
#         plt.title(est)
#
#         plt.xlim(xmin=0)
#         plt.ylim(ymin=0)
#
#         plt.tight_layout()
#         plt.savefig(f'procars_correlation_{alignment_type}_{est.replace(" ", "_").replace("/", "_")}.pdf')
#         plt.show()

def scatter_correlation_pairs():
    # genomes = "7c"
    # res = 500
    # infercars_file = f'real_data/ucsc_{genomes}_diff_res/{res}k/Conserved.Segments'
    # alignment_type = 'nuc_alpha_05'
    # distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'

    infercars_file = f'real_data/procars/orthology_blocks.txt'
    alignment_type = 'nuc_05'
    distances_file = f'real_data/align_distances/distances_procars_{alignment_type}.csv'

    dh = DistanceHolder(distances_file)

    ests_ds, al_ds = count_alignment_pairs_alignment(infercars_file, dh)

    for est, est_ds in ests_ds.items():
        plt.figure()
        sns.set(style="whitegrid", font="serif")

        X = np.array(est_ds).reshape(-1, 1)
        y = np.array(al_ds).reshape(-1, 1)

        ols = linear_model.LinearRegression(fit_intercept=False)
        model = ols.fit(X, y)
        response = model.predict(X)

        r2 = model.score(X, y)

        plt.plot(X, response, color='k', label='Regression model')
        plt.scatter(X, y, edgecolor='k', facecolor='grey', alpha=0.7, label='Sample data')

        plt.annotate("$R^2$ = {:.3f}".format(r2), (50, 0.7))

        plt.xlabel('Estimated k')
        plt.ylabel('Distance in alignment')
        # plt.legend(loc='best')
        plt.title(est)

        plt.xlim(xmin=0)
        plt.ylim(ymin=0)

        plt.tight_layout()
        plt.savefig(f'procars_correlation_{alignment_type}_{est.replace(" ", "_").replace("/", "_")}_no_intercept.pdf')
        plt.show()


def diff_combinations_500k():
    infercars_folder = '/Users/alexey/Downloads/DESCHRAMBLER/examples/500k.maml.ucsc_7_%s/Conserved.Segments'

    files = ['simple', 'no_horse', 'no_dog', 'no_cat', 'baboon', 'pig_cow', 'pig_sheep', 'cow_sheep', 'pig_cow_sheep',
             'chicken']
    ys = defaultdict(list)

    for file in files:
        print(file)
        infercars_file = infercars_folder % file
        errors = count_tree_errors(tree_file, infercars_file,
                                   error_func=lambda y, y_tree: math.sqrt(np.sum(np.power((y - y_tree), 2))))
        for est, error in errors.items():
            if est == 'Corr dirichlet 1/3': continue
            ys[est].append(error)

    sns.set(style="whitegrid", font="serif")
    for est, ys in ys.items():
        plt.plot(files, ys, label=f'{est}')

    # plt.title(f'Tree error')
    plt.xlabel('set of genomes')
    plt.ylabel('tree error')
    plt.legend(loc=2)
    plt.xticks(rotation=20)

    plt.tight_layout()
    plt.savefig(f'ucsc_diff_comb.pdf')
    plt.show()


def boxplot_errors():
    genomes = "procars"
    tree_file = 'real_data/procars/tree.nwk'
    infercars_file = 'real_data/procars/orthology_blocks.txt'

    # genomes = "7c"
    # infercars_file = f'real_data/ucsc_{genomes}_diff_res/410k/Conserved.Segments'
    #

    # suffix = '_5_left'

    # res = 5000
    # tree_file = f'real_data/Yersinia_pestis/tree{suffix}.nwk'
    # infercars_file = f'real_data/Yersinia_pestis/{res}/blocks_unique_coords{suffix}.infercars'

    errors = count_tree_errors(tree_file, infercars_file, cyclic=False,
                               error_func=lambda y, y_tree: (y - y_tree) / y)
                               # error_func=lambda y, y_tree: np.power(y - y_tree, 2) / y)

    errs_2d = [[est, er] for est, ers in errors.items() for er in ers if er == er and er != np.inf]

    for est, errs in errors.items():
        print(est)
        # filt_errs = list(filter(lambda v: v == v and v != np.inf, errs))
        print('median:', np.mean(errs))
        print('mean:', np.mean(errs))
        print('sum:', np.sum(errs))
        print()

    df = pd.DataFrame(errs_2d, columns=['Estimator', 'Relative Error, $(y - y_{tree}) / y$'])

    sns.set(style="whitegrid", font="serif")
    sns.boxplot(data=df, x='Estimator', y='Relative Error, $(y - y_{tree}) / y$')

    plt.tight_layout()
    plt.savefig(f'boxplot_procars_y_tree_relative.pdf')
    plt.show()

# scatter_correlation_pairs()
boxplot_errors()
