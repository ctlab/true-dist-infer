import numpy as np
from scipy import optimize


class TannierEstimator:
    @staticmethod
    def predict(d, b, c2):
        def e_c1(n, k):
            return n * sum(
                np.prod([(- 2 * k / (n + u)) for u in range(0, l)])
                for l in range(50))

        def dc1_dn(n, k):
            return e_c1(n, k) / n - n * sum(np.prod([(- 2 * k / (n + u))
                                                     for u in range(l)])
                                            * sum(1 / (n + u) for u in range(l))
                                            for l in range(50))

        def dc1_dk(n, k):
            return n * sum(
                l * np.prod([(- 2 * k / (n + u)) for u in range(l)]) / k
                for l in range(50))

        def e_c2(n, k):
            return k * n * n * sum(
                (l + 1) * (m + 1) /
                (4 * (k - 1) ** 2 * np.prod([(n + u) / (- 2 * (k - 1)) for u in range(l + m + 2)]))
                for l, m in np.ndindex((50, 50)))

        def dc2_dn(n, k):
            return e_c2(n, k) / n * 2 - k * n * n * sum(
                (l + 1) * (m + 1) /
                (4 * (k - 1) ** 2 * np.prod([(n + u) / (- 2 * (k - 1)) for u in range(l + m + 2)])) *
                sum(1 / (n + u) for u in range(l + m + 2))
                for l, m in np.ndindex((50, 50)))

        def dc2_dk(n, k):
            return e_c2(n, k) / k + k * n * n * sum(
                (l + 1) * (m + 1) * (l + m) /
                (4 * (k - 1) ** 3 * np.prod([(n + u) / (- 2 * (k - 1)) for u in range(l + m + 2)]))
                for l, m in np.ndindex((50, 50)))

        def create_fun(real_b, real_c2):
            def fun(x):
                return [x[0] - e_c1(x[0], x[1]) - real_b,
                        e_c2(x[0], x[1]) - real_c2]

            return fun

        def jac(x):
            return np.array([[1 - dc1_dn(x[0], x[1]),
                              - dc1_dk(x[0], x[1])],
                             [dc2_dn(x[0], x[1]),
                              dc2_dk(x[0], x[1])]])

        fun = create_fun(b, c2)
        prediction = optimize.root(fun, np.array([3 * b, d]), jac=jac, method='hybr')

        # print(prediction)
        return prediction.x[0], prediction.x[1]

    @staticmethod
    def predict_k(d, b, c2):
        return TannierEstimator.predict(d, b, c2)[1]

    @staticmethod
    def predict_n(d, b, c2):
        return TannierEstimator.predict(d, b, c2)[0]
