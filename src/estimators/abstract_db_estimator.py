from scipy import optimize


class AbstractDBEstimator:
    def d_over_n(self, x):
        raise NotImplementedError('subclasses must override d_over_n!')

    def b_over_n(self, x):
        raise NotImplementedError('subclasses must override b_over_n!')

    def predict(self, g):
        d_over_b = lambda r: lambda x: self.d_over_n(x) / self.b_over_n(x) - r

        d = g.d()
        b = g.b()
        if d >= b - 1e-6: d = d - 1e-1

        print('d, b:', d, b)
        if d > 0 and b > 0 and d / b > 1/2:
            x = optimize.bisect(d_over_b(d / b), 1e-6, 3, xtol=1e-4)
            print('d/b:', d/b)
        else:
            x = 1e-6

        b_n = self.b_over_n(x)
        n = b / b_n

        print('estimated:', max(d, n * x / 2))
        print()
        return n, max(d, n * x / 2)

    def predict_k(self, g):
        return self.predict(g)[1]

    def predict_n(self, g):
        return self.predict(g)[0]