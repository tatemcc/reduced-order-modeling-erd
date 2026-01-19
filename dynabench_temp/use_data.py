from dynabench.dataset import DynabenchIterator

advection_iterator = DynabenchIterator(equation='advection',
                                       structure='cloud',
                                       resolution='low',
                                       lookback=4,
                                       rollout=16)

for sample in advection_iterator:
    x, y, points = sample

    # x is the input data with shape (lookback, n_points, n_features)
    # y is the target data with shape (rollout, n_points, n_features)
    # points are the observation points with shape (n_points, dim)
    # for the advection equation n_features=1 and dim=2
