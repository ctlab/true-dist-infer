import numpy as np

from collections import defaultdict

from src.estimators.abstract_db_estimator import AbstractDBEstimator

cache = defaultdict(lambda: defaultdict(lambda: {}))


def c_m(m, t):
    def f(x):
        if m in cache and t in cache[m] and x in cache[m][t]:
            return cache[m][t][x]
        if m == 1:
            return t ** t / (x + t) ** t
        else:
            prev = c_m(m - 1, t)(x)
            pr = 1
            for j in range(m - 2):
                pr *= (m * t + m + j) / (m * t - t + m - 1 + j)
            pr *= ((m - 1) * t + 2 * m - 4) / m
            ret = prev * x * t ** (t + 1) / (x + t) ** (t + 2) * pr
            if not (m in cache and t in cache[m] and x in cache[m][t]):
                cache[m][t][x] = ret
            return ret

    return f


b_ = lambda t: lambda x: 1 - c_m(1, t)(x)
d_ = lambda t, mx: lambda x: 1 - np.sum(c_m(m, t)(x) for m in range(1, mx))


class DirichletDBEstimator(AbstractDBEstimator):
    def __init__(self, alpha, mx=50):
        self.t = alpha
        self.mx = mx

    def d_over_n(self, x):
        return d_(self.t, self.mx)(x)

    def b_over_n(self, x):
        return b_(self.t)(x)
