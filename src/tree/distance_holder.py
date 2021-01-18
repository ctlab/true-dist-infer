import pandas as pd

from string import digits

remove_digits = str.maketrans('', '', digits)


class DistanceHolder:
    def __init__(self, file):
        self.df = pd.read_csv(file, index_col=0)
        self.df.fillna(0, inplace=True)

        self.actual_assembly = {a.translate(remove_digits): a for a in self.df.columns}

    def get_dist(self, sp1, sp2):
        def get_column_name(sp):
            return self.actual_assembly[sp.translate(remove_digits)]

        cl1, cl2 = get_column_name(sp1), get_column_name(sp2)
        return self.df[cl1][cl2] + self.df[cl2][cl1]
