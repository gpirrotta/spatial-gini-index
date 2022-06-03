import numpy as np
from sgi.weights import RangeDistanceBand



def trapezoidify_area(bases, heights):
    
    bases = np.asarray(bases)
    bases = bases/ (bases.sum())
    cumulative_bases = bases.cumsum()
    total_bases_sum = np.r_[cumulative_bases[0],np.resize(cumulative_bases[1:],cumulative_bases.shape) + cumulative_bases][:-1]

    heights = np.asarray(heights)
    heights = heights / heights.sum()
    cumulative_heights = heights.cumsum()
    total_heights_differences = np.diff(np.r_[[0], cumulative_heights])

    area = (np.sum(total_heights_differences * total_bases_sum) / 2)

    return area, total_bases_sum, total_heights_differences


def get_lagged_target(sgi):
    distance_matrix = sgi.distance_matrix.copy()
    np.fill_diagonal(distance_matrix, np.inf)
    min_threshold=0
    max_threshold = sgi.h_distances[0]
    w = RangeDistanceBand.create_W(distance_matrix, min_threshold, max_threshold)
    y_lagged = np.ones(distance_matrix.shape[0])
    for row, columns in w.neighbors.items():
        tmp = []
        for column in columns:
            tmp.append(sgi.target[column])
        y_lagged[row] = np.mean(tmp)

    return y_lagged

