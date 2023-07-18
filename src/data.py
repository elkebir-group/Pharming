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
    seg_to_snvs: dict #a dictionary that maps each seg to a list of snvs in that segment
    # cell_lookup : pd.Series # an N length series mapping internal cell index to input cell label
    # mut_lookup : pd.Series #an M length series mapping interal SNV index to input SNV label

    def __post_init__(self):
        self.nseg = len(self.seg_to_snvs)
        self.N, self.M = self.var.shape
        self.segments = list(self.seg_to_snvs.keys())
        self.segments.sort()
        self.cells = np.arange(self.N)
        self.muts = np.arange(self.M)


        # Create the multi-index
        multi_index = pd.MultiIndex.from_product([self.cells, self.muts], names=['cell', 'snv'])


        # Flatten the array to create the Series data
        self.total_series = pd.Series(self.total.flatten(), index=multi_index)
        self.total_series = self.total_series[self.total_series > 0]


    
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

    def count_cells_by_snv(self, seg):
        snvs = self.seg_to_snvs[seg]
        cells_by_snvs = {j: [] for j in snvs}
        filtered_series = self.total_series.loc[pd.IndexSlice[:, snvs]]
        for i,j in filtered_series.index:
              cells_by_snvs[j].append(i)

        # indices = np.argwhere(self.total > 0)


        # for i,j in indices:
        #     if j not in snvs:
        #         continue 
        #     print(f"cell{i}, snv: {j}  {self.total[i,j]}")
          

        return np.count_nonzero(self.total[:,snvs],axis=0), cells_by_snvs
    
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
    
    def cn_states_by_seg(self, seg):
        cn_states = []
        total_cn= np.unique(self.copy_numbers[:,seg])
        for cn in total_cn:
            cn_states.append((cn-1,1))
        return cn_states



            # print(f"{value}: count={count}, indices={list(indices)}")



