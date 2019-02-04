from scipy import optimize


class AbstractDBEstimator:
    def d_over_n(self, x):
        raise NotImplementedError('subclasses must override d_over_n!')

    def b_over_n(self, x):
        raise NotImplementedError('subclasses must override b_over_n!')

    def predict(self, d, b):
        d_over_b = lambda r: lambda x: self.d_over_n(x) / self.b_over_n(x) - r

        x = optimize.bisect(d_over_b(d / b), 1e-6, 3, xtol=1e-4)

        b_n = self.b_over_n(x)
        n = b / b_n
        return n, n * x / 2

    def predict_k(self, d, b):
        return self.predict(d, b)[1]

    def predict_n(self, d, b):
        return self.predict(d, b)[0]