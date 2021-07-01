import matplotlib.pyplot as plt
import seaborn as sns

from itertools import combinations
from collections import defaultdict

from src.graphs.real_data_graph import RealDataGraph
from src.utils.parsers import parse_to_df


file = '/Users/alexey/Downloads/DESCHRAMBLER/examples/maml.ucsc/%dk/Conserved.Segments'
species = ['hg38', 'mm10', 'rn6', 'canFam4']

ks = list(range(100, 510, 10))

ys = []

for k in ks:
    df = parse_to_df(file % k)
    ys.append(len(df.block.unique()))

sns.set(style="whitegrid", font="serif")

plt.plot(ks, ys)

plt.title(f'All blocks count for different minimal block size')
plt.xlabel('minimal rearrangement size, kilobase')
plt.ylabel('blocks count')
# plt.legend(loc=4)

plt.savefig(f'all.new-len.pdf')
plt.show()



