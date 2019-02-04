from scipy.special import hyp2f1

from src.estimators.abstract_db_estimator import AbstractDBEstimator


class FlatDirichletDBEstimator(AbstractDBEstimator):
    def d_over_n(self, x):
        return 1 - (1 + x) ** 2 * (hyp2f1(-2 / 3, -1 / 3, 1 / 2, 27 * x / (4 * (1 + x) ** 3)) - 1) / (3 * x)

    def b_over_n(self, x):
        return x / (1 + x)
