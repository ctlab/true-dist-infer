import scipy.stats


def get_probabilities_by_distribution(distribution, params, size):
    if distribution is None:
        return [1 / size for _ in range(size)]
    else:
        dist = getattr(scipy.stats, "gamma")
        return dist.rvs(*params, size=size)