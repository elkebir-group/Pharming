from dataclasses import dataclass
import numpy as np
import pandas as pd
import pickle
from segment_genome import Segment
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


    def binomial_likelihood(self, cells, snvs, vaf, alpha=0.001):
  

        var =  self.total[np.ix_(cells, snvs)]
        total =self.total[np.ix_(cells, snvs)]
     

        adj_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
        adj_vaf = adj_vaf.reshape(1, -1)
        cellprobs= -1*binom.logpmf(var, total, p=adj_vaf)
            
        cellprobs = np.nansum(cellprobs, axis=1)

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


    def cn_proportions(self, seg):
        cn_props = {}
        states, counts = self.cn_states_by_seg(seg)
        for cn in states:
            cn_props[cn] = counts[cn]/self.N
        return cn_props
        



    def save(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)

def load_from_files(var_fname, copy_fname ):
    col_names = ['segment', 'mutation_label', 'cell_label','var', 'total']

   
    read_counts = pd.read_table(
        var_fname, header=None, names=col_names, skiprows=[0])
    


    copy_numbers = pd.read_csv(copy_fname, header=None,names=["segment", "cell_label", "x", "y"], skiprows=[0])

    cell_labels_not_in_df1 = copy_numbers.loc[~copy_numbers['cell_label'].isin(read_counts['cell_label']), 'cell_label']
    excluded_labels = cell_labels_not_in_df1.unique()
    print(f"Excluding {len(excluded_labels)} cells due to lack of read counts:")
    for c in excluded_labels:
        print(c)

    #exclude any cells that have no read counts at any position
    copy_numbers = copy_numbers[copy_numbers['cell_label'].isin(read_counts['cell_label'])]
  
    return load(read_counts, copy_numbers)

def load(read_counts,copy_numbers ):


    cell_labels = np.sort(read_counts['cell_label'].unique())

    mut_labels = np.sort(read_counts['mutation_label'].unique())



    #create indexed series of mapping of cell index to label
    cell_lookup = pd.Series(data=cell_labels, name="cell_label").rename_axis("cell")     
    mut_lookup = pd.Series(data=mut_labels, name="mutation_label").rename_axis("mut")

    
    read_counts = pd.merge(read_counts, cell_lookup.reset_index(), on='cell_label', how='left')
    read_counts = pd.merge(read_counts, mut_lookup.reset_index(), on='mutation_label', how='left')

    #in long format
    segs = copy_numbers.loc[:, ["segment"]].drop_duplicates()
    seg_labels = np.sort(segs['segment'].unique())
    seg_lookup = pd.Series(data=seg_labels, name="segment").rename_axis("seg_id")

    copy_numbers = pd.merge(copy_numbers, cell_lookup.reset_index(), on='cell_label', how='left').drop("cell_label", axis=1)
    copy_numbers= pd.merge(copy_numbers, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)

   

    read_counts = pd.merge(read_counts, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)
    seg_to_mut_mapping = read_counts.loc[:, ["seg_id", "mut"]].drop_duplicates()
    snv_to_seg = seg_to_mut_mapping.set_index("mut")["seg_id"].to_dict()
    seg_to_snvs =  {value: [k for k, v in snv_to_seg.items() if v == value] for value in set(snv_to_seg.values())}


    read_counts= read_counts.set_index(["cell", "mut"])
    var = read_counts["var"].unstack(level="mut", fill_value=0).to_numpy()
    total = read_counts["total"].unstack(level="mut", fill_value=0).to_numpy()



    copy_numbers["cell"] = copy_numbers["cell"].astype(int)
    copy_numbers = copy_numbers.set_index([ "cell", "seg_id"])

    copy_x = copy_numbers["x"].unstack(level="seg_id").to_numpy()
    copy_y = copy_numbers["y"].unstack(level="seg_id").to_numpy()


    return Data(var, total, copy_x, copy_y, snv_to_seg, seg_to_snvs, cell_lookup, mut_lookup)

def segment(cn_profiles, tol=0.0, pseudo_counts=1e-6):
    '''
    cn_profiles: an numpy array with shape N x B (cell by bin) matrix with the total copy number profiles of for each cell in each bin
    
    returns: a pd.DataFrame with Segment ID as index and columns start and end bin (inclusive)
    '''
    return Segment(tol=tol, pseudo_counts=pseudo_counts).fit(cn_profiles)


def load_from_pickle(fname):
    return pd.read_pickle(fname)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-m", "--mut-lookup", type=str,
        help="filename of mutation label to export")
    parser.add_argument("-l", "--cell-lookup", type=str,
        help="filename of cell label to export")
    # parser.add_argument("-i", "--index", type=str, default="cell",
    #     help="index of copy number profiles")
    # parser.add_argument("-t", "--tolerance",  type=float, default=0.0,
    #     help="lower bound KL-divergence to establish a new segment")
    # parser.add_argument("-p", "--pseudo", required=False, type=float, default=1e-6,
    #     help="pseudo counts to use for computing KL- divergence")
    parser.add_argument("-D", "--data",  type=str, 
        help="filename of pickled data object")
    
    args = parser.parse_args()
    # folder = "s12_m1000_k25_l7"
    # cell_folder= "n5000_c0.01_e0.15"
    # pth = f"simulation_study/input/{folder}/{cell_folder}"
    # args = parser.parse_args([
    #     "-f", f"{pth}/sparse.p0",
    #     "-c",f"{pth}/cells.p0",
    #     "-D", f"{pth}/data.pickle"
    #     ])

    dat = load_from_files(args.file, args.profiles)
    if args.data is not None:
        dat.save(args.data)
    
    if args.mut_lookup is not None:
        dat.export_mut_lookup(args.mut_lookup)
    
    if args.cell_lookup is not None:
        dat.export_cell_lookup(args.cell_lookup)
