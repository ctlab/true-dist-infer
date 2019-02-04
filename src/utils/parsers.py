import re
import numpy as np
import pandas as pd

p = "([A-Za-z0-9]+)\.([A-Za-z0-9_]+):(\d+)-(\d+) ([+|-]).*"
pattern = re.compile(p)
columns = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]


def parse_to_df(file_name):
    def find_indices(lst, condition):
        return [i for i, elem in enumerate(lst) if condition(elem)]

    with open(file_name) as f:
        lines = f.readlines()

    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []

    for i, b in enumerate(bs):
        b = b[1:-1]
        for oc in b:
            m = pattern.match(oc)
            temp.append([i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])

    return pd.DataFrame(temp, columns=columns)


def parse_grimm_file(file):
    lss = []
    with open(file, "r") as f:
        for ch_index, line in enumerate(filter(lambda l: len(l) > 1 and (l[-2] == "$" or l[-2] == "@"), f)):
            split = list(map(lambda x: int(x), line.split()[:-1]))
            lss.append(split)
    return lss