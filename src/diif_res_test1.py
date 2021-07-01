import matplotlib.pyplot as plt
import seaborn as sns

from itertools import combinations
from collections import defaultdict

from src.graphs.real_data_graph import RealDataGraph
from src.utils.parsers import parse_to_df
from src.tree.fit_to_tree import filter_alt_


file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/maml.ucsc_7/%dk/Conserved.Segments'


ys = defaultdict(list)
ks = list(range(200, 1010, 25))

for k in ks:
    df = parse_to_df(file % k)
    df = filter_alt_(df)

    species = df.species.unique()
    print(k)
    for sp1, sp2 in combinations(species, 2):
        g = RealDataGraph()
        g.infercars(df, sp1, sp2)
        print('   ', sp1, sp2, g.b(), g.d())

        ys[(sp1, sp2)].append(g.d() / g.b())

sns.set(style="whitegrid", font="serif", font_scale=0.7)


for (sp1, sp2), ys in ys.items():
    plt.plot(ks, ys, label=f'{sp1}-{sp2}')

plt.title(f'All d/b for different minimal block size')
plt.xlabel('minimal rearrangement size, kilobase')
plt.ylabel('d/b')
plt.legend(loc=4)

plt.savefig(f'ucsc_7-d_b.pdf')
plt.show()



