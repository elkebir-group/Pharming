from dataclasses import dataclass
import numpy as np
import pandas as pd
import pickle
from collections import defaultdict, Counter
from scipy.stats import binom
import utils 
from scipy.sparse import lil_matrix

'''
N = number of cells
M = number of SNVs
G = number of segments
'''
#we expect some observed VAFs to be np.NaN due to d_{ij} =0 in ulc coverage data
np.seterr(invalid='ignore')
from scipy.sparse import csr_matrix

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
    seg_lookup: pd.Series
    alpha : float = 0.001 #per base sequencing error rate 

    def __post_init__(self):
        self.nseg = len(self.seg_to_snvs)
        self.N, self.M = self.var.shape
        self.segments = list(self.seg_to_snvs.keys())
        self.segments.sort()
        self.cells = np.arange(self.N)
        self.muts = np.arange(self.M)
        self.likelihood_dict = self.precompute_likelihood(self.alpha)


        # Create the multi-index
        multi_index = pd.MultiIndex.from_product([self.cells, self.muts], names=['cell', 'snv'])


        # Flatten the array to create the Series data
        self.total_series = pd.Series(self.total.flatten(), index=multi_index)
        self.total_series = self.total_series[self.total_series > 0]


    
    def __str__(self):
        mystr= f"Input data contains: {self.N} cells, {self.M} SNVs and {self.nseg} segments"
        return mystr

    def to_sparse(self):
 

        # Assuming var and total are numpy arrays
        # Find the indices where total is nonzero
        nonzero_indices = np.nonzero(self.total)

        # Extract the non-zero elements and their corresponding row and column indices
        total_nonzero = self.total[nonzero_indices]
        var_nonzero = self.var[nonzero_indices]

        # Create a sparse matrix using csr_matrix
        self.total_sparse = csr_matrix((total_nonzero, nonzero_indices), shape=self.total.shape)

        # Repeat the same for the var matrix if needed
        self.var_sparse = csr_matrix((var_nonzero, nonzero_indices), shape=self.var.shape)
    

    def unique_vafs(self):
        all_vals  = self.copy_x + self.copy_y
        copy_numbers = np.unique(all_vals)
        vafs = []
        for total in copy_numbers[copy_numbers > 0]:
            for alt in range(total+1):
                vafs.append(alt/total)
   
        self.vafs = set(vafs)
    
    def compute_cell_likelihoods(self, vaf_to_snvs, cells=None):
        if cells is None:
            cells = self.cells
        
        cell_likes = np.zeros(len(cells))
        for v, snvs in vaf_to_snvs.items():
            cell_likes += self.likelihood_dict[v][cells[:, np.newaxis], snvs].sum(axis=1)

        return cell_likes
    
    def compute_snv_likelihoods(self, v, snvs, cells=None):
        if cells is None:
            return self.likelihood_dict[v][:,snvs].sum(axis=0)
        
        return self.likelihood_dict[v][cells[:, np.newaxis],snvs].sum(axis=0)
      
    
    def compute_cmb(self, cells, snvs):
        var = self.var[cells,:][:,snvs]
        nom = np.count_nonzero(var, axis=1)
        total = self.total[cells,:][:,snvs]
        denom = np.count_nonzero(total, axis=1)
      
        return nom/denom


    def precompute_likelihood(self, alpha=0.001):
        self.to_sparse()
        self.unique_vafs()

        # Initialize a dictionary to store likelihood matrices
        self.likelihood_dict = {}
        nonzero_indices = self.total_sparse.nonzero()

        # Convert VAF array to numpy array
        vaf = np.array([v for v in self.vafs])
        adj_vaf = vaf * (1 - alpha) + (1 - vaf) * (alpha / 3)
        for v in self.vafs:
            self.likelihood_dict[v] = lil_matrix(self.total_sparse.shape, dtype=float)
        var, total = [], []
        for i, j in zip(*self.total_sparse.nonzero()):
            var.append(self.var[i,j])
            total.append(self.total[i,j])
        
        
    
            # Compute the logpmf for the current entry

        for i,v in enumerate(adj_vaf):
            self.likelihood_dict[vaf[i]][nonzero_indices] =  -1 * binom.logpmf( var, total, p=v)
            self.likelihood_dict[vaf[i]] = self.likelihood_dict[vaf[i]].toarray()

        return self.likelihood_dict
        

        #     # Update the likelihood matrices for each VAF
        #     for k, v in enumerate(self.vafs):
        #         likelihood_dict[v][i, j] = logpmf_entries_per_vaf[k]


        
         
  

# # If needed, convert the result back to a sparse matrix
# logpmf_sparse = csr_matrix((logpmf_values, nonzero_indices), shape=total.shape)
            
#             # Compute logpmf for each entry (var[i, j], total[i, j])
#             logpmf_entry = binom.logpmf(self.var_sparse, self.total_sparse, p=adj_vaf)
            
        #     # Append the logpmf values for the current vaf_value
        #     logpmf_values.append(logpmf_entry)

        # # Convert the list of logpmf values to a numpy array
        # logpmf_values_array = np.array(logpmf_values)
        # print("done")





 
    # def binomial_likelihood2(self, snvs, vaf, alpha=0.001):
    #      vaf = np.atleast_1d(vaf)
    #      var = self.var[:, snvs]
    #      total = self.total[:,snvs]
    #      var = var[:, np.newaxis, :]
    #      total = total[:, np.newaxis, :]
    #      adj_vaf = vaf * (1 - alpha) + (1 - vaf) * (alpha / 3)
    #      logpmf_values= -binom.logpmf(var, total, p=adj_vaf)
    #      cell_probs = np.sum(logpmf_values, axis=2)
    #      return cell_probs


    def binomial_likelihood(self, cells, snvs, vaf, alpha=0.001, axis=1):
        vaf = np.atleast_1d(vaf)

        cells = np.asarray(cells)  # Convert to NumPy array
        if cells.ndim != 1:
            raise ValueError("Cells must be a 1-dimensional array")
        var = self.var[cells[:, np.newaxis], snvs]
        total = self.total[cells[:, np.newaxis], snvs]

        adj_vaf = vaf * (1 - alpha) + (1 - vaf) * (alpha / 3)
        adj_vaf = adj_vaf.reshape(1, -1)
        cellprobs = -binom.logpmf(var, total, p=adj_vaf)

        cellprobs = np.nansum(cellprobs, axis=axis)
 

        return cellprobs
    # def binomial_likelihood(self, cells, snvs, vaf, alpha=0.001, axis=1):
  
    #     # if isinstance(vaf, float) or isinstance(vaf, int):
    #     #     vaf = np.array(vaf)

    #     vaf = np.atleast_1d(vaf)

    #     var =  self.var[np.ix_(cells, snvs)]
    #     total =self.total[np.ix_(cells, snvs)]
    #     # var = self.var[cells[:, None], snvs]
    #     # total = self.total[cells[:, None], snvs]
     

    #     adj_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
    #     adj_vaf = adj_vaf.reshape(1, -1)
    #     cellprobs= -1*binom.logpmf(var, total, p=adj_vaf)
            
    #     cellprobs = np.nansum(cellprobs, axis=axis)

    #     return cellprobs
    
    def var_marg(self, snvs):
        return self.var[:,snvs].sum(axis=0)
    
    def total_marg(self, snvs):
        return self.total[:,snvs].sum(axis=0)
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
    
    def num_snvs(self, ell):
        return len(self.seg_to_snvs[ell])

    def count_cells_by_snv(self, seg):
        snvs = self.seg_to_snvs[seg]
        cells_by_snvs = {j: [] for j in snvs}
        filtered_series = self.total_series.loc[pd.IndexSlice[:, snvs]]
        for i,j in filtered_series.index:
              cells_by_snvs[j].append(i)

        return np.count_nonzero(self.total[:,snvs],axis=0), cells_by_snvs
    
    def copy_profiles_by_seg(self, segments, cells=None):
        if cells is None:
            return self.copy_x[:, segments], self.copy_y[:, segments]
        else:
            segments = np.array(segments)
   
            return self.copy_x[cells[:, np.newaxis], segments], self.copy_y[cells[:, np.newaxis], segments]
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
    
    def thresholded_cn_prop(self, seg, thresh=0, start_state=(1,1) ):
        cn_props = self.cn_proportions(seg)
        cn_states = [state for state, prop in cn_props.items() if prop >= thresh or state==start_state]
        
        norm = sum(cn_props[s] for s in cn_states)
        norm_cn_props = {state: cn_props[state]/norm for state in cn_states}
        
        if start_state not in norm_cn_props:
            norm_cn_props[start_state] =0 
        return norm_cn_props

    
    def cn_states_by_seg(self, seg):

        x = self.copy_x[:,seg]
        y = self.copy_y[:,seg]
        # cn_states = set([(x, y) for x, y in zip(x,y)])
        cn_states = [(x, y) for x, y in zip(x,y)]
        item_counts = Counter(cn_states)

        return set(cn_states), item_counts 

    def num_cn_states(self, seg, thresh=0):
        cn_prop =  self.cn_proportions(seg)

        return len([s for s in cn_prop if cn_prop[s] >= thresh])

    def cn_proportions(self, seg):
        cn_props = {}
        states, counts = self.cn_states_by_seg(seg)
        for cn in states:
            cn_props[cn] = counts[cn]/self.N
        return cn_props
        

    # def get_largest_segments(self,n =2,  min_cn_states=3):
    #     cand_segments = [ell for ell in self.segments if self.num_cn_states(ell) >= min_cn_states ]
    #     if len(cand_segments) < n:
    #         return cand_segments
    #     cand_segments = sorted(cand_segments, reverse=True, key=lambda x: len(self.seg_to_snvs[x]))
    #     return cand_segments[:n]


    def save(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)




