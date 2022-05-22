import numpy as np
from tqdm import tqdm

class DistanceMetric():

    @classmethod
    def haversine(cls, points_coordinates):
        data = np.deg2rad(points_coordinates)                     
        lat = data[:,0]                     
        lng = data[:,1]         
        diff_lat = lat[:,None] - lat
        diff_lng = lng[:,None] - lng
        d = np.sin(diff_lat/2)**2 + np.cos(lat[:,None])*np.cos(lat) * np.sin(diff_lng/2)**2

        return 2 * 6371 * np.arcsin(np.sqrt(d))


    @classmethod
    def euclidean(cls, x, y):
        """ Computing pairwise distances using memory-efficient
        vectorization.

        Parameters
        ----------
        x : numpy.ndarray, shape=(M, D)
        y : numpy.ndarray, shape=(N, D)

        Returns
        -------
        numpy.ndarray, shape=(M, N)
            The Euclidean distance between each pair of
            rows between `x` and `y`."""
        sqr_dists = -2 * np.matmul(x, y.T)
        sqr_dists +=  np.sum(x**2, axis=1)[:, np.newaxis]
        sqr_dists += np.sum(y**2, axis=1)
        return  np.sqrt(np.clip(sqr_dists, a_min=0, a_max=None))
