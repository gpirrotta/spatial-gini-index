import dataclasses
from tqdm import tqdm
from scipy.sparse import csr_matrix
import numpy as np
from sgi.distance import DistanceMetric
from sgi.weights import RangeDistanceBand
from sgi.reports import Report                

import pandas as pd
import scipy.stats
import timeit

np.random.seed(42)

#example of decay function
def decay_function(distance, h_distance, s, phi=5, alpha=2):
    #return (phi * (distance[s[0,0],s[:,1]]/h_distance)**alpha).tolist()
    return (distance[s[0,0],s[:,1]]**-2).tolist()

class SGI:
    '''
    Spatial Gini Index 
    (Mucciardi et al.)
    '''
    def __init__(self, points, y, step='max-min', labels=None, permutations=0):
        self._metric = 'SGI'
        self._points = points
        self._y = y
        self._step = step
        self._w = None
        self._links = []
        self._variabilities = []
        self._h_distances = []
        self._sgi = None
        self._report = None
        self._report_pvalue = None
        self._labels = labels
        self._pvalue_z = None
        self._pvalue_mc = None

        self._n = points.shape[0]
        self._distance_matrix = None
        self._sgis = []
        self._permutations = permutations
        self._time_sgi = 0.0
        self._time_pvalue = 0.0
        
        self._build()

    
    def _build(self):
        start_time = timeit.default_timer()

        self._distance_matrix = DistanceMetric.euclidean(self._points, self._points)
        self._target_matrix = self._build_target_matrix(self._y)
        self._h_distances, self._links, self._variabilities = self._build_weights()
        self._sgi = self._build_sgi(self._variabilities)

        self._time_sgi= timeit.default_timer() - start_time

        if self._permutations > 0:
            start_time = timeit.default_timer()
            self._pvalue_z, self._pvalue_mc = self._build_pvalue()
            self._time_pvalue= timeit.default_timer() - start_time
        


    def _build_weights(self):

        if self._step == 'max-min':
            return self._build_weights_maxmin()
        elif self._step == 'constant':
            return self._build_weights_constant()
        
        

    
    def _build_weights_maxmin(self):
        
        distance_matrix = self._distance_matrix.copy()

        np.fill_diagonal(distance_matrix, np.inf)
        num_links = np.inf
        dim = distance_matrix.shape[0]
        total_links = dim * (dim - 1)
        progressive_num_links = 0
        min_threshold=0
        disable = False

        h_distances = []
        links = []
        variabilities = []


        with tqdm(total=total_links, disable=disable) as pbar:            
            while progressive_num_links < total_links:
                pbar.set_description(f"Processing {progressive_num_links} links")
                
                minimums = np.min(distance_matrix, axis=0)
                minimus = minimums[minimums != np.inf]
                max_threshold = np.max(minimus)
                
                filter_distance = distance_matrix <= max_threshold
                num_links = filter_distance.sum().sum()
                
                variability = self._spatial_lag_total(distance_matrix, min_threshold, max_threshold)

                distance_matrix[filter_distance] = np.inf
                progressive_num_links += num_links
                min_threshold = max_threshold

                h_distances.append(max_threshold)
                links.append(num_links)
                variabilities.append(variability)
                pbar.update(num_links)
            
        return h_distances, links, variabilities

        
    def _build_weights_constant(self):
        distance_matrix = self._distance_matrix.copy()
        np.fill_diagonal(distance_matrix, np.inf)

        spatial_lag = 0
        num_links = np.inf
        dim = distance_matrix.shape[0]
        total_links = dim * (dim - 1)
        progressive_num_links = 0
        min_threshold = 0

        h_distances = []
        links = []
        variabilities = []

        with tqdm(total=total_links) as pbar:
            while progressive_num_links < total_links:
                pbar.set_description(f"Processing {progressive_num_links} links")
                spatial_lag += 1
                max_threshold = 0
                if spatial_lag == 1:
                    self._delta = 0
                    minimums = np.min(distance_matrix, axis=0)
                    minimus = minimums[minimums != np.inf]
                    self._delta = np.max(minimus)
                    max_threshold = self._delta
                else:
                    max_threshold = self._delta * spatial_lag

                filter_distance = distance_matrix <= max_threshold
                num_links = filter_distance.sum().sum()

                variability = self._spatial_lag_total(distance_matrix, min_threshold, max_threshold)

                distance_matrix[filter_distance] = np.inf
                progressive_num_links += num_links
                min_threshold = max_threshold

                h_distances.append(max_threshold)
                links.append(num_links)
                variabilities.append(variability)
                pbar.update(num_links)
            
        return h_distances, links, variabilities


    def _build_pvalue(self):
        h_distances  = [0]+self._h_distances
        ids = np.arange(self._n)
        sgis = np.zeros((self._permutations, ))
        disable = False

        with tqdm(total=self._permutations, disable=disable) as pbar:            
            for perm in range(self._permutations):
                pbar.set_description(f"Processing permutation {perm}")
                np.random.shuffle(ids)
                distance_matrix = self._distance_matrix.copy()
                np.fill_diagonal(distance_matrix, np.inf)
                y_shuffled = self._y[ids]
                #y_shuffled = self._y 
                variabilities = []
                for min_threshold,max_threshold in zip(h_distances[0:], h_distances[1:]):
                    variability = 0
                    
                    filter_distance = distance_matrix <= max_threshold
                    w = RangeDistanceBand.create_W(distance_matrix, min_threshold, max_threshold)
                    distance_matrix[filter_distance] = np.inf

                    for row, columns in w.neighbors.items():
                        for column in columns:
                            variability+=(y_shuffled[row]-y_shuffled[column])**2

                    variabilities.append(variability)
                
                sgis[perm] = self._build_sgi(variabilities)
                pbar.update(1)

        sgis = np.asarray(sgis)
        self._sgis = sgis
        sqm = np.std(sgis)
        z_cal = (self._sgi-0.5)/(sqm/np.sqrt(self._permutations))        
        self._z_cal = z_cal
        pvalue_z = scipy.stats.norm.sf(abs(z_cal))*2

    
        above = sgis >= self._sgi
        larger = above.sum()
        if (self._permutations - larger) < larger:
            larger = self._permutations - larger

        pvalue_mc = (larger + 1.) / (self._permutations + 1.)

        return pvalue_z, pvalue_mc

    
    def _spatial_lag_total(self, distance_matrix, min_threshold, max_threshold, target_matrix=None):
        
        w = RangeDistanceBand.create_W(distance_matrix, min_threshold, max_threshold)
        if target_matrix is None:
            target_matrix= self._target_matrix

        shape = distance_matrix.shape
        flatted_weights = np.concatenate(list(w.weights.values()))
        sparse_matrix = csr_matrix((flatted_weights, w.neighbors_indexes), shape=shape)
        variability = sparse_matrix.multiply(target_matrix).sum()

        return variability
    
        

    def _build_sgi(self, variabilities):            
        links = np.asarray(self._links)
        links = links / links.sum()
        cumulative_links = links.cumsum()
        links_total_differences = np.diff(np.r_[[0], cumulative_links])

        variabilities = np.asarray(variabilities)
        variabilities = variabilities/ (variabilities.sum())
        cumulative_variabilities = variabilities.cumsum()
        J_total = np.r_[cumulative_variabilities[0],np.resize(cumulative_variabilities[1:],cumulative_variabilities.shape) + cumulative_variabilities][:-1]
                
        sgi = 1 - (np.sum(links_total_differences * J_total) / 2)

        return sgi



    def _build_target_matrix(self, y):
        return (y[:, np.newaxis] - y) ** 2


    @property
    def metric(self):
        return self._metric

    @property
    def h_distances(self):
        return self._h_distances

    @property
    def links(self):
        return self._links

    @property
    def variabilities(self):
        return self._variabilities

    @property
    def sgi(self):
        return self._sgi

    @property
    def distance_matrix(self):
        return self._distance_matrix
            
    @property
    def report(self):
        report = Report(self)
    
        return report.report

    @property
    def report_pvalue(self):
        report = Report(self)
        return report.report_pvalue

    
    @property
    def pvalue_z(self):
        if self._pvalue_z is None:
            print("PValue Z not calculated")
        else:
            return self._pvalue_z


    @property
    def pvalue_mc(self):
        if self._pvalue_mc is None:
            print("PValue MC not calculated")
        else:
            return self._pvalue_mc




class SGD(SGI):
    '''
    Spatial Gini Decomposition
    (Rey et al.)
    '''
    def __init__(self, points, y, step='max-min', labels=None, permutations = 0):
        SGI.__init__(self,points, y, step, labels, permutations)
        self._metric = 'SGD'
        

    def _target_matrix(self):
        return np.abs(self._y[:, np.newaxis] - self._y)

    def _build_sgi(self, variabilities):
        # In Rey Decomposition SGI is not calculated
        pass

