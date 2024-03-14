from dataclasses import dataclass
import numpy as np
import pandas as pd
import pickle
import argparse
from collections import defaultdict, Counter
from scipy.stats import binom

'''
N = number of cells
M = number of SNVs
G = number of segments
'''
#we expect some observed VAFs to be np.NaN due to d_{ij} =0 in ulc coverage data
np.seterr(invalid='ignore')

@dataclass
class Data:
    var : np.array  # N x M matrix with the number of variant read counts
    total : np.array #  N x M matrix with the number of total read counts
    # copy_numbers: pd.DataFrame# N rows with columns x, y  and MultIndex (seg_id,cell)  where (x,y) is the copy number state of each cell in each segment
    copy_x: np.array #n x k matrix with the copy number of the maternal allele
    copy_y: np.array #n x k matrix with the copy number of the paternal allele

    #use copy_numbers.loc[seg_id] to subset rows based on segment
    #use copy_numbers.xs(cell, level='cell') to subset rows based on cell
    #use copy_numbers.loc[(seg_id, cell)] to access tuple (x,y)
    snv_to_seg: dict  # a dictionary that maps each snv to seg
    seg_to_snvs: dict #a dictionary that maps each seg to a list of snvs in that segment
    cell_lookup : pd.Series # an n length series mapping internal cell index to input cell label
    mut_lookup : pd.Series #an m length series mapping interal SNV index to input SNV label

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


    def binomial_likelihood(self, cells, snvs, vaf, alpha=0.001, axis=1):
  
        if isinstance(vaf, float) or isinstance(vaf, int):
            vaf = np.array(vaf)

        var =  self.var[np.ix_(cells, snvs)]
        total =self.total[np.ix_(cells, snvs)]
     

        adj_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
        adj_vaf = adj_vaf.reshape(1, -1)
        cellprobs= -1*binom.logpmf(var, total, p=adj_vaf)
            
        cellprobs = np.nansum(cellprobs, axis=axis)

        return cellprobs
    
    def compute_vafs(self, cells=None, snvs=None):
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return np.sum(var, axis=0)/np.sum(total, axis=0)
    

    
    def obs_vafs(self, cells=None, snvs=None):
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return var/total
    
    def nonzero(self, axis):

        varcount = np.count_nonzero(self.var, axis=axis)
        totcount = np.count_nonzero(self.total, axis=axis)
        return varcount, totcount
    

    def export_mut_lookup(self, fname):
        df = pd.DataFrame({"index": self.mut_lookup.index, "label": self.mut_lookup.values}).reset_index(drop=True)
        df.to_csv(fname, index=False)
    
    def export_cell_lookup(self, fname):
        df = pd.DataFrame({"index": self.cell_lookup.index, "label": self.cell_lookup.values}).reset_index(drop=True)
        df.to_csv(fname, index=False)


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

        return np.count_nonzero(self.total[:,snvs],axis=0), cells_by_snvs
    
    def copy_profiles_by_seg(self, segments, cells=None):
        '''
        returns a pandas dataframe  with cols ['x','y' ]subsetted to specified segment and cell set
        '''
        if cells is None:
            return self.copy_x[:,segments], self.copy_y[:,segments]
        else:
            return self.copy_x[np.ix_(cells,segments)], self.copy_y[np.ix_(cells,segments)]
        # if cells is None:
        #     return self.copy_numbers.loc[segment]
        # else:
        #     return self.copy_numbers.loc[(segment, cells), :]
    

    def cells_by_cn(self, seg):
        copy_x_segment = self.copy_x[:, seg]
        copy_y_segment = self.copy_y[:, seg]
        combined_copy = np.column_stack((copy_x_segment, copy_y_segment))
        copy_states_dict = defaultdict(list)
        unique_copy_states = set(map(tuple, combined_copy))

        # Iterate through unique copy states
        for state in unique_copy_states:
            # Find the indices where the state occurs in combined_copy
            indices = np.where(np.all(combined_copy == state, axis=1))
            
            # Append the indices to the dictionary
            copy_states_dict[state] = list(indices[0])

        # Convert the defaultdict to a regular dictionary
        return dict(copy_states_dict)

    
    def cn_states_by_seg(self, seg):

        x = self.copy_x[:,seg]
        y = self.copy_y[:,seg]
        # cn_states = set([(x, y) for x, y in zip(x,y)])
        cn_states = [(x, y) for x, y in zip(x,y)]
        item_counts = Counter(cn_states)

        return set(cn_states), item_counts 

    def num_cn_states(self, seg, root_state=None):
        cn_states, counts = self.cn_states_by_seg(seg)
        # cn_states = set(root_state).intersection(cn_states)
        return len(cn_states)

    def cn_proportions(self, seg):
        cn_props = {}
        states, counts = self.cn_states_by_seg(seg)
        for cn in states:
            cn_props[cn] = counts[cn]/self.N
        return cn_props
        

    def get_largest_segments(self,n =2,  min_cn_states=3):
        cand_segments = [ell for ell in self.segments if self.num_cn_states(ell) >= min_cn_states ]
        if len(cand_segments) < n:
            return cand_segments
        cand_segments = sorted(cand_segments, reverse=True, key=lambda x: len(self.seg_to_snvs[x]))
        return cand_segments[:n]


    def save(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)




