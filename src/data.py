"""A module to model the input data for Pharming."""

from collections import Counter
from dataclasses import dataclass
import pickle
import numpy as np
import pandas as pd
from scipy.stats import binom
from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix



#we expect some observed VAFs to be np.NaN due to d_{ij} =0 in ulc coverage data
np.seterr(invalid='ignore')


@dataclass
class Data:
    """The Pharming data class.

    This class is used to store the input data for Pharming. It contains the following attributes:
    Attributes
    ----------
    var : np.array
        N x M matrix with the number of variant read counts
    total : np.array
        N x M matrix with the number of total read counts
    copy_x : np.array
        N x G matrix with the copy number of the maternal allele
    copy_y : np.array
        N x G matrix with the copy number of the paternal allele
    snv_to_seg : dict
        A dictionary that maps each SNV to a segment    
    seg_to_snvs : dict
        A dictionary that maps each segment to a list of SNVs in that segment
    cell_lookup : pd.Series
        An N length series mapping internal cell index to input cell label    
    mut_lookup : pd.Series
        An M length series mapping internal SNV index to input SNV label
    seg_lookup : pd.DataFrame
        A G length dataframe mapping internal segment index to input segment label
    alpha : float
        Per base sequencing error rate
    seg_weights : np.ndarray
        The weight of each segment for the objective function  
    n : int
        The number of cells
    m : int
        The number of SNVs
    nseg : int
        The number of segments
    segments : list
        A list of segments
    cells : np.array
        An array of cell indices
    muts : np.array
        An array of SNV indices
    vafs : set
        A set of unique VAFs in the data
    likelihood_dict : dict  
        A dictionary that maps each VAF to a likelihood matrix

    """ 

    var : np.ndarray
    total : np.ndarray
    copy_x: np.ndarray
    copy_y: np.ndarray
    snv_to_seg: dict
    seg_to_snvs: dict
    cell_lookup : pd.Series
    mut_lookup : pd.Series
    seg_lookup: pd.DataFrame
    alpha : float = 0.001
    seg_weights : np.ndarray = None


    def __post_init__(self):
        """Set key attributes of the data object from the input data."""
        self.nseg = len(self.seg_to_snvs)
        self.N, self.m = self.var.shape
        self.segments = list(self.seg_to_snvs.keys())
        self.segments.sort()
        self.cells = np.arange(self.N)
        self.muts = np.arange(self.m)
        if self.seg_weights is None:
            self.seg_weights = np.ones(self.nseg)

         # Find the indices where total is nonzero
        nonzero_indices = np.nonzero(self.total)

        # Extract the non-zero elements and their corresponding row and column indices
        total_nonzero = self.total[nonzero_indices]

        # Create a sparse matrix using csr_matrix
        total_sparse = csr_matrix((total_nonzero, nonzero_indices), shape=self.total.shape)

        copy_numbers  = np.unique(self.copy_x + self.copy_y)
        vafs = []
        for total in copy_numbers[copy_numbers > 0]:
            for alt in range(total+1):
                vafs.append(alt/total)
   
        self.vafs = set(vafs)

          # Initialize a dictionary to store likelihood matrices
        self.likelihood_dict = {}
        nonzero_indices = total_sparse.nonzero()

        # Convert VAF array to numpy array
        vaf = np.array([v for v in self.vafs])
        adj_vaf = vaf * (1 - self.alpha) + (1 - vaf) * (self.alpha / 3)
        for v in self.vafs:
            self.likelihood_dict[v] = lil_matrix(total_sparse.shape, dtype=float)
        var, total = [], []
        for i, j in zip(*total_sparse.nonzero()):
            var.append(self.var[i,j])
            total.append(self.total[i,j])
        
        
        # Compute the logpmf for the current entry
        for i,v in enumerate(adj_vaf):
            self.likelihood_dict[vaf[i]][nonzero_indices] =  -1 * binom.logpmf( var, total, p=v)
            self.likelihood_dict[vaf[i]] = self.likelihood_dict[vaf[i]].toarray()

        # Create the multi-index
        # multi_index = pd.MultiIndex.from_product([self.cells, self.muts], names=['cell', 'snv'])

        # Flatten the array to create the Series data
        # self.total_series = pd.Series(self.total.flatten(), index=multi_index)
        # self.total_series = self.total_series[self.total_series > 0]

    
    def __str__(self):
        """Return a string representation of the data object."""
        mystr= f"Input data contains: {self.N} cells, {self.m} SNVs and {self.nseg} segments"
        return mystr
   
    def compute_cell_likelihoods(self, vaf_to_snvs, cells=None):
        """Compute the likelihood of each cell given the VAFs and SNVs.

        Parameters
        ----------
        vaf_to_snvs : dict
            A dictionary that maps each VAF to a list of SNVs
        cells : np.array, optional
            An array of cell indices to compute the likelihood for
        
        Returns
        -------
        cell_likes : np.array
            An array of likelihood values for each cell
        """
        if cells is None:
            cells = self.cells
        
        cell_likes = np.zeros(len(cells))
        for v, snvs in vaf_to_snvs.items():
            cell_likes += self.likelihood_dict[v][cells[:, np.newaxis], snvs].sum(axis=1)

        return cell_likes
    
    def compute_snv_likelihoods(self, v, snvs, cells=None):
        """Compute the likelihood of each SNV given the VAF and a listlike object of cells."""
        if cells is None:
            return self.likelihood_dict[v][:,snvs].sum(axis=0)
        
        return self.likelihood_dict[v][cells[:, np.newaxis],snvs].sum(axis=0)
      
    def compute_cmb(self, cells, snvs):
        """Compute the per cell CMB for a given listlike object of SNVs."""
        var = self.var[cells,:][:,snvs]
        nom = np.count_nonzero(var, axis=1)
        total = self.total[cells,:][:,snvs]
        denom = np.count_nonzero(total, axis=1)
      
        return nom/denom

    def max_alleles(self, ell):
        """Return the maximum number of alleles in a given list of segments."""

  
        x = self.copy_x[:,ell].max()
        y = self.copy_y[:,ell].max()
            
        return max(x,y)
    
    def sample_segments(self, k, rng=None,  max_cn_states=1, min_snvs =0, thresh=0, max_cn=10):
        """Randomly sample k segments with at least min_snvs and at most max_cn_states."""
        candidates = [ell for ell in self.seg_to_snvs if self.num_snvs(ell) >= min_snvs and 
                      self.num_cn_states(ell, thresh) <= max_cn_states and self.max_alleles(ell) <= max_cn]
        
        print(f"Total candidates: {len(candidates)}")
        if rng is None:
            rng = np.random.default_rng()
        samples = list(rng.choice(candidates, k, replace=False))
        return samples


    def binomial_likelihood(self, cells, snvs, vaf, alpha=0.001, axis=1):
        """Compute the binomial likelihood of the data."""
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
   
    
    def var_marg(self, snvs):
        """Compute the marginal variant count for a given listlike object of SNVs."""
        return self.var[:,snvs].sum(axis=0)
    
    def total_marg(self, snvs):
        """Compute the marginal total count for a given listlike object of SNVs."""
        return self.total[:,snvs].sum(axis=0)
    
    def compute_vafs(self, cells=None, snvs=None):
        """Compute the VAFs for a given subset of cells and SNVs."""
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return np.sum(var, axis=0)/np.sum(total, axis=0)
    
    def obs_vafs(self, cells=None, snvs=None):
        """Compute the observed VAFs for a given subset of cells and SNVs."""
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return var/total
    
    def export_mut_lookup(self, fname):
        """Export the mutation lookup table to a CSV file."""
        df = pd.DataFrame({"index": self.mut_lookup.index, "label": self.mut_lookup.values}).reset_index(drop=True)
        df.to_csv(fname, index=False)
    
    def export_cell_lookup(self, fname):
        """Export the cell lookup table to a CSV file."""
        df = pd.DataFrame({"index": self.cell_lookup.index, "label": self.cell_lookup.values}).reset_index(drop=True)
        df.to_csv(fname, index=False)

    def export_segment_lookup(self, fname):
        """Export the segment lookup table to a CSV file."""
        print(self.seg_lookup.head())
        # df = pd.DataFrame({"index": self.seg_lookup.index, "label": self.seg_lookup.values}).reset_index(drop=True)
        self.seg_lookup.to_csv(fname, index=False)

    def num_snvs(self, ell):
        if ell not in self.seg_to_snvs:
            return 0
        """Return the number of SNVs in a given segment."""
        return len(self.seg_to_snvs[ell])
    
    def copy_profiles_by_seg(self, segments, cells=None):
        """Return the copy number profiles for a given list of segments."""
        if cells is None:
            return self.copy_x[:, segments], self.copy_y[:, segments]
        else:
            segments = np.array(segments)
   
            return self.copy_x[cells[:, np.newaxis], segments], self.copy_y[cells[:, np.newaxis], segments]

    def save(self, fname):
        """Save the data object to a pickle file."""
        with open(fname, 'wb') as file:
            pickle.dump(self, file)
    
    def thresholded_cn_prop(self, seg, thresh:float=0, start_state=(1,1), include_start_state=True ):
        """Return the proportion of cells in each copy number state for a given segment above a threshold."""
        cn_props = self.cn_proportions(seg)
        cn_states = [state for state, prop in cn_props.items() if prop >= thresh]
        
        norm = sum(cn_props[s] for s in cn_states)
        norm_cn_props = {state: cn_props[state]/norm for state in cn_states}
        
        if include_start_state and start_state not in norm_cn_props:
            norm_cn_props[start_state] =0 
        return norm_cn_props
    
    def cn_states_by_seg(self, seg):
        """Return the copy number states for a given segment."""
        x = self.copy_x[:,seg]
        y = self.copy_y[:,seg]
        cn_states = [(x, y) for x, y in zip(x,y)]
        item_counts = Counter(cn_states)

        return set(cn_states), item_counts 

    def num_cn_states(self, seg, thresh=0, start_state=(1,1), include_start_state=False):
        """Return the number of copy number states for a given segment."""
        cn_prop =  self.thresholded_cn_prop(seg, thresh, start_state, include_start_state)

        return len(cn_prop)

    def cn_proportions(self, seg):
        """Return the proportion of cells in each copy number state for a given segment."""
        cn_props = {}
        states, counts = self.cn_states_by_seg(seg)
        for cn in states:
            cn_props[cn] = counts[cn]/self.N
        return cn_props
    
    def consensus_profile(self, cells, ell):
        """Return the consensus copy number profile for a given segment."""
        x = self.copy_x[cells, ell]
        y = self.copy_y[cells, ell]
        states = list(zip(x, y))
        cntr = Counter(states)
        consensus_state = cntr.most_common(1)
        return consensus_state[0][0]

    
    # def count_marginals(self, seg):
    #     """Compute the marginal variant and total counts for a given segment."""
    #     cell_map = self.cells_by_cn(seg)
    #     snvs = self.seg_to_snvs[seg]
    #     alt = {cn: self.var[np.ix_(cells, snvs)].sum(axis=0) for cn,cells in cell_map.items()}
    #     total = {cn: self.total[np.ix_(cells, snvs)].sum(axis=0) for cn,cells in cell_map.items()}

    #     return snvs, alt, total
    

    # def count_cells_by_snv(self, seg):
    #     snvs = self.seg_to_snvs[seg]
    #     cells_by_snvs = {j: [] for j in snvs}
    #     filtered_series = self.total_series.loc[pd.IndexSlice[:, snvs]]
    #     for i,j in filtered_series.index:
    #           cells_by_snvs[j].append(i)

    #     return np.count_nonzero(self.total[:,snvs],axis=0), cells_by_snvs

        # def nonzero(self, axis):
    #     """Compute the number of non-zero entries in the data matrix."""
    #     varcount = np.count_nonzero(self.var, axis=axis)
    #     totcount = np.count_nonzero(self.total, axis=axis)
    #     return varcount, totcount
    

        # if cells is None:
        #     return self.copy_numbers.loc[segment]
        # else:
        #     return self.copy_numbers.loc[(segment, cells), :]
    

    # def cells_by_cn(self, seg):
    #     copy_x_segment = self.copy_x[:, seg]
    #     copy_y_segment = self.copy_y[:, seg]
    #     combined_copy = np.column_stack((copy_x_segment, copy_y_segment))
    #     copy_states_dict = defaultdict(list)
    #     unique_copy_states = set(map(tuple, combined_copy))

    #     # Iterate through unique copy states
    #     for state in unique_copy_states:
    #         # Find the indices where the state occurs in combined_copy
    #         indices = np.where(np.all(combined_copy == state, axis=1))
            
    #         # Append the indices to the dictionary
    #         copy_states_dict[state] = list(indices[0])

    #     # Convert the defaultdict to a regular dictionary
    #     return dict(copy_states_dict)
    


    # def get_largest_segments(self,n =2,  min_cn_states=3):
    #     cand_segments = [ell for ell in self.segments if self.num_cn_states(ell) >= min_cn_states ]
    #     if len(cand_segments) < n:
    #         return cand_segments
    #     cand_segments = sorted(cand_segments, reverse=True, key=lambda x: len(self.seg_to_snvs[x]))
    #     return cand_segments[:n]






