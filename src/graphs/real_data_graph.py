from src.graphs.abstract_graph import AbstractGraph
from src.utils.parsers import parse_grimm_file

uniq_predicate = lambda ls: lambda b: ls.count(b) == 1
sp_blocks_from_df = lambda df, sp: df.loc[df['species'] == sp]['block'].tolist()


class RealDataGraph(AbstractGraph):
    def infercars(self, df, sp1, sp2, cyclic=False):
        def add_one_sp(label, sp):
            df_sp = df.loc[df['species'] == sp].sort_values(by=['chr', 'chr_beg'])

            for chr in df_sp['chr'].unique():
                df_sp_chr = df_sp.loc[df_sp['chr'] == chr]
                bs = [row['block'] if row['orientation'] == "+" else -row['block'] for _, row in df_sp_chr.iterrows()]
                self.add_edges_from_list(list(filter(lambda b: abs(b) in allowed_blocks, bs)), label, cyclic)

        super().__init__()
        allowed_blocks = sp_blocks_from_df(df, sp1)
        allowed_blocks = list(filter(uniq_predicate((sp_blocks_from_df(df, sp2))), allowed_blocks))
        add_one_sp("red", sp1)
        add_one_sp("black", sp2)

    def build_grimm(self, file1, file2, cyclic=False):
        def add_from_bss(label, bss):
            for bs in bss:
                self.add_edges_from_list(list(filter(lambda b: abs(b) in allowed_blocks, bs)), label, cyclic)

        bss1 = parse_grimm_file(file1)
        bss2 = parse_grimm_file(file2)

        bss1_flat = list([abs(b) for bs in bss1 for b in bs])
        bss2_flat = list([abs(b) for bs in bss2 for b in bs])

        allowed_blocks = range(min(bss1_flat), max(bss1_flat))
        allowed_blocks = list(filter(lambda b: bss1_flat.count(b) == 1, allowed_blocks))
        allowed_blocks = list(filter(lambda b: bss2_flat.count(b) == 1, allowed_blocks))

        add_from_bss("red", bss1)
        add_from_bss("black", bss2)

    def add_edges_from_list(self, ls, label, cyclic=False):
        def add_block_edge(i, j):
            self.add_edge(str(abs(i)) + ("h" if i > 0 else "t"), str(abs(j)) + ("t" if j > 0 else "h"), label=label)

        self.add_nodes_from(map(lambda x: str(abs(x)) + "h", ls))
        self.add_nodes_from(map(lambda x: str(abs(x)) + "t", ls))
        for i, j in zip(ls, ls[1:] + ls[:1]) if cyclic else zip(ls[:-1], ls[1:]):
            add_block_edge(i, j)

    def construct_bp_graph(self):
        return self
