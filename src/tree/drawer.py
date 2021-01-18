import math

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from itertools import combinations
from collections import defaultdict


from src.tree.fit_to_tree import count_tree_errors, count_correlation_pairs_alignment, count_alignment_pairs_alignment
from src.tree.distance_holder import DistanceHolder
tree_file = '/Users/alexey/PycharmProjects/true-dist-infer/src/tree/hg38.100way.tree'


def diff_ks_ucsc_7():
    # infercars_file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/maml.ucsc_7/%dk/Conserved.Segments'
    # infercars_file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/chicken/%dk/Conserved.Segments'
    genomes = "7c"
    infercars_file = f'real_data/ucsc_{genomes}_diff_res/%dk/Conserved.Segments'

    alignment_type = 'nuc_alpha_05'
    distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'
    dh = DistanceHolder(distances_file)

    ys = defaultdict(list)
    # ks = list(range(100, 1010, 25))
    ks = list(range(50, 1001, 10))

    for k in ks:
        print(k)
        # errors = count_tree_errors(tree_file, infercars_file % k, error_func=lambda y, y_tree: math.sqrt(np.sum(np.power((y - y_tree), 2))))
        errors = count_tree_errors(tree_file, infercars_file % k, dh, alignment_type)
        # print(k, errors[alignment_type] > errors['Parsimony'], errors[alignment_type] - errors['Parsimony'])
        print(k, errors['Uniform'] > 0.15, errors['Uniform'])
        for est, error in errors.items():
            ys[est].append(error)

    sns.set(style="whitegrid", font="serif")

    for est, ys in ys.items():
        plt.plot(ks, ys, label=f'{est}', ls=('--' if est == alignment_type else '-'))

    # plt.title(f'Tree error')
    plt.xlabel('minimal rearrangement size, kilobase')
    plt.ylabel('tree error')
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(f'ucsc_{genomes}_ks_tree_rel.pdf')
    plt.show()


def diff_ks_ucsc_7_correlation_pairs():
    genomes = "7c"
    infercars_file = f'real_data/ucsc_{genomes}_diff_res/%dk/Conserved.Segments'
    alignment_type = 'nuc_alpha_05_logdet'
    distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'

    dh = DistanceHolder(distances_file)
    ys = defaultdict(list)
    # ks = list(range(50, 2001, 300))
    ks = list(range(100, 1001, 100))

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
    plt.savefig(f'ucsc_{genomes}_ks_tree_corr_{alignment_type}.pdf')
    plt.show()


def scatter_correlation_pairs():
    genomes = "7cz"
    res = 500
    infercars_file = f'real_data/ucsc_{genomes}_diff_res/{res}k/Conserved.Segments'
    alignment_type = 'nuc_alpha_05'
    distances_file = f'real_data/align_distances/distances_{alignment_type}.csv'

    dh = DistanceHolder(distances_file)

    ests_ds, al_ds = count_alignment_pairs_alignment(infercars_file, dh)

    for est, est_ds in ests_ds.items():
        plt.figure()
        sns.set(style="whitegrid", font="serif")
        sns.regplot(est_ds, al_ds, label=f'{est}', scatter_kws={'s':3}, line_kws={'lw': 0.5})

        plt.xlabel('Estimated k')
        plt.ylabel('Distance in alignment')
        # plt.legend(loc='best')

        plt.tight_layout()
        plt.savefig(f'scatter_{genomes}_{res}k_{alignment_type}_{est.replace(" ", "_").replace("/", "_")}.pdf')
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


diff_ks_ucsc_7()

