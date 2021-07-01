from src.utils.parsers import parse_to_df, export_df_to_infercars


from ete3 import Tree
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict

import numpy as np

tree_file = 'real_data/Yersinia_pestis/tree.nwk'
tree_file_out = 'real_data/Yersinia_pestis/tree_%s_left.nwk'

res = 1000
infercars_file = f'real_data/Yersinia_pestis/{res}/blocks_unique_coords.infercars'
infercars_file_out = f'real_data/Yersinia_pestis/{res}/blocks_unique_coords_%s_left.infercars'

t = Tree(tree_file)

n = len(t.get_leaves())

m = np.zeros((n, n))
for i, leaf1 in enumerate(t.get_leaves()):
    for j, leaf2 in enumerate(t.get_leaves()):
        m[i, j] = t.get_distance(leaf1, leaf2)

cls = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average',
                              distance_threshold=0.0004).fit_predict(m)

print(cls)
print(np.unique(cls))

used = defaultdict(bool)
survivors = []

for cl, leaf in zip(cls, t.get_leaves()):
    if not used[cl]:
        used[cl] = True
        survivors.append(leaf)

t.prune(survivors)
print(len(survivors))
t.write(outfile=tree_file_out % len(survivors), format=5)

df = parse_to_df(infercars_file)
print('df before', len(df))

df = df.loc[df['species'].isin([s.name for s in survivors])].copy()
print('df after', len(df))

export_df_to_infercars(df, infercars_file_out % len(survivors))