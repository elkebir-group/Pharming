from dataclasses import dataclass
import numpy as np
import pandas as pd



'''
N = number of cells
M = number of SNVs
G = number of segments
'''

@dataclass
class Data:
    var : np.array  # N x M matrix with the number of variant read counts
    total : np.array #  N x M matrix with the number of total read counts
    copy_numbers: np.array# N x G matrix with the copy number profile of each cell in each segment
    snv_to_seg: dict  # a dictionary that maps each snv to seg
    seg_to_snvs: dict #a ditionary that maps each seg to a list of snvs in that segment
    # cell_lookup : pd.Series # an N length series mapping internal cell index to input cell label
    # mut_lookup : pd.Series #an M length series mapping interal SNV index to input SNV label

    def __post_init__(self):
        self.nseg = len(self.seg_to_snvs)
        self.N, self.M = self.var.shape
        self.segments = list(self.seg_to_snvs.keys())
        self.segments.sort()
        self.cells = np.arange(self.N)
        self.muts = np.arange(self.M)
    
    def __str__(self):
        mystr= f"Input data contains: {self.N} cells, {self.M} SNVs and {self.nseg} segments"
        return mystr

    def compute_vafs(self, cells=None, snvs=None):
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return var.sum(axis=0)/total.sum(axis=0)
    
    def count_marginals(self, seg):
        cell_map = self.cells_by_cn(seg)
        snvs = self.seg_to_snvs[seg]
        alt = {cn: self.var[np.ix_(cells, snvs)].sum(axis=0) for cn,cells in cell_map.items()}
        total = {cn: self.total[np.ix_(cells, snvs)].sum(axis=0) for cn,cells in cell_map.items()}

        return snvs, alt, total

    
    def compute_likelihood(self):
        pass 

    def cells_by_cn(self, seg):
        arr= self.copy_numbers[:,seg]
        mapping = {}
        unique_values, inverse_indices, value_counts = np.unique(arr, return_inverse=True, return_counts=True)
        indices_per_value = [np.where(inverse_indices == i)[0] for i in range(len(unique_values))]
        for value,  indices in zip(unique_values, indices_per_value):
            mapping[value] = list(indices)
        return mapping

            # print(f"{value}: count={count}, indices={list(indices)}")



