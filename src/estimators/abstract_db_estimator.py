from scipy import optimize


class AbstractDBEstimator:
    def d_over_n(self, x):
        raise NotImplementedError('subclasses must override d_over_n!')

    def b_over_n(self, x):
        raise NotImplementedError('subclasses must override b_over_n!')

    def predict(self, g, quiet=False):
        d_over_b = lambda r: lambda x: self.d_over_n(x) / self.b_over_n(x) - r

        d = g.d()
        b = g.b()
        if d >= b - 1e-6: d = d - 1e-1

        is_nan = False

        if not quiet: print('d, b:', d, b)
        if d > 0 and b > 0 and d / b > 1/2:
            try:
                x = optimize.bisect(d_over_b(d / b), 1e-6, 3, xtol=1e-4)
            except ValueError:
                # print("NANING")
                is_nan = True
                x = float('nan')
            if not quiet: print('d/b:', d/b)
        else:
            x = 1e-6

        b_n = self.b_over_n(x)
        n = b / b_n

        if not quiet: print('estimated:', max(d, n * x / 2))
        if not quiet: print()

        if is_nan: return x, x
        return n, max(d, n * x / 2)

    def predict_k(self, g, quiet=False):
        return self.predict(g, quiet=quiet)[1]

    def predict_n(self, g, quiet=False):
        return self.predict(g, quiet=quiet)[0]