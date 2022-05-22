import numpy as np

class W:
    def __init__(self, neighbors, weights):
        self._neighbors = neighbors
        self._weights = weights
        self._neighbors_indexes = None
        self._s0 = np.sum([len(v) for v in weights.values()])
        self._build_neighbors_indexes()
        self._transform = 'B'

    @property
    def neighbors(self):
        return self._neighbors

    @property
    def weights(self):
        return self._weights

    @property
    def s0(self):
        return self._s0

    @property
    def neighbors_indexes(self):
        return self._neighbors_indexes


    @property
    def transform(self):
        return self._transform

    @transform.setter
    def transform(self, transformation):
        if self._transform == transformation:
            pass
        elif transformation == 'R':
            self._transform='R'
            for row, columns in self._weights.items():
                cols = np.array(columns)/len(columns)
                self._weights[row] = cols.tolist()
        elif transformation == 'B':
            self._transform='B'
            for row, columns in self._weights.items():
                cols = np.ones(len(columns)).tolist()
                self._weights[row] = cols
        
                
            

    def _build_neighbors_indexes(self):
        rows = []
        cols = []
        for row, columns in self._neighbors.items():
            for col in columns:
                rows.append(row)
                cols.append(col)
        self._neighbors_indexes = tuple([np.array(rows), np.array(cols)])
        

class RangeDistanceBand:

    @classmethod
    def create_W(cls, distance_matrix, min_threshold, max_threshold, binary=True, decay=None):
        np.fill_diagonal(distance_matrix, np.inf)
        filter = np.logical_and(distance_matrix > min_threshold, distance_matrix <= max_threshold)
       
        neighbors_indexes =  np.where(filter)       

        tmp = np.asarray(neighbors_indexes).T
        sort_arr = tmp[tmp[:, 0].argsort(), :]
        split_arr = np.split(sort_arr, np.where(np.diff(sort_arr[:,0]))[0] + 1)

        neighbors = {s[0,0]: s[:,1].tolist() for s in split_arr}

        if binary:
            weights = {s[0,0]: np.ones(s[:,1].shape[0]).tolist() for s in split_arr}
        else:
            weights = {s[0,0]: decay(distance_matrix, max_threshold, s) for s in split_arr}

        return W(neighbors, weights)

