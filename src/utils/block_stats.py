from bisect import bisect_left, bisect_right

every_chr = lambda df, chr_f: [r for chr in df['chr'].unique() for r in chr_f(df.loc[df['chr'] == chr])]


def dist_between_blocks(df, with_begs=False):
    def chr_f(df):
        ss = sorted(df['chr_beg'].tolist())
        es = sorted(df['chr_end'].tolist())
        ans = [ss[0]] if with_begs else []
        for e in es:
            l = bisect_left(ss, e)
            r = bisect_right(ss, e)
            if r != len(ss):  # and ss[l] - e != 0:
                ans.append(ss[r] - e)
        return ans
    return every_chr(df, chr_f)