import math

from src.estimators.abstract_db_estimator import AbstractDBEstimator


class UniformDBEstimator(AbstractDBEstimator):
    def d_over_n(self, x):
        c_m = lambda m: lambda x: math.exp(- m * x) * (x ** (m - 1)) * (m ** (m - 2)) / math.factorial(m)
        return 1 - sum(c_m(m)(x) for m in range(1, 100))

    def b_over_n(self, x):
        return 1 - math.exp(- x)
