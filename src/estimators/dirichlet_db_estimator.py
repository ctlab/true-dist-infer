import numpy as np

from collections import defaultdict
from scipy import optimize

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
d_b = lambda t, mx: lambda x: d_(t, mx)(x) / b_(t)(x)


class DirichletDBEstimator(AbstractDBEstimator):
    def __init__(self, alpha, mx=50):
        self.t = alpha
        self.mx = mx

    def d_over_n(self, x):
        return d_(self.t, self.mx)(x)

    def b_over_n(self, x):
        return b_(self.t)(x)


class CorrectedDirichletDBEstimator(DirichletDBEstimator):
    def __init__(self, t, mx=50):
        super().__init__(t, mx)

    def predict(self, g):
        addings = lambda chr, b: lambda x: (-0.33 + ((x - 0.4) ** (-14 / 11)) * 5 / 100) * chr / b * x / (1 + x)
        chr_ = g.chr()
        d = g.d()
        b = g.b()
        d_over_b_corrected = lambda r: lambda x: d_b(self.t, self.mx)(x) - r \
                                                 + (0 if x <= 0.5 else  # chr / b * x / (1 + x)
                                                   x / 1000 - 0.0045
                                                   + (- 0.0015 if chr_ > 0 else 0)
                                                   + (- 0.002 if self.t <= 0.7 else 0)
                                                   + (- 0.002 if self.t <= 0.5 else 0)
                                                   + addings(chr_, b)(x))

        x = optimize.bisect(d_over_b_corrected(d / b), 1e-6, 3, xtol=1e-4)

        b_n = self.b_over_n(x)
        n = b / b_n
        return n, n * x / 2
