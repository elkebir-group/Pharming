from dataclasses import dataclass
import numpy as np
import pandas as pd
import pickle
from segment_genome import Segment
import argparse
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
    copy_numbers: np.array# N x G matrix with the copy number profile (x,y) of each cell in each segment
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

    def compute_vafs(self, cells=None, snvs=None):
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return var.sum(axis=0)/total.sum(axis=0)
    
    def obs_vafs(self, cells=None, snvs=None):
        if cells is None:
            cells = self.cells 
        if snvs is None:
            snvs = self.muts 
        
        var =self.var[np.ix_(cells, snvs)]
        total = self.total[np.ix_(cells, snvs)]
        return var/total
    

    
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
    
    def compute_likelihood(self):
        pass 

    def cells_by_cn(self, seg):
        arr= self.copy_numbers[:,seg]
        mapping = {}
        unique_values, inverse_indices, value_counts = np.unique(arr, return_inverse=True, return_counts=True)
        indices_per_value = [np.where(inverse_indices == i)[0] for i in range(len(unique_values))]
        for value,  indices in zip(unique_values, indices_per_value):
            state = tuple(value)
            mapping[state] = list(indices)
        return mapping
    
    def cn_states_by_seg(self, seg):
        all_states = self.copy_numbers[:,seg]
        all_states =[tuple(state) for state in all_states]
    
        # df = self.copy_numbers[self.copy_numbers["seg_id"]==seg]
        # unique_rows =df.drop_duplicates(subset=['x', 'y'])
        # cn_states = [(x, y) for x, y in zip(unique_rows['x'], unique_rows['y'])]

        return set(all_states)






    def save(self, fname):
        with open(fname, 'wb') as file:
            pickle.dump(self, file)

def load_from_files(read_counts_fname,copy_numbers_fname ):
        col_names = ['chr', 'segment', 'mutation_label', 'cell_label','var', 'total']


    
        read_counts = pd.read_table(
            read_counts_fname, sep="\t", header=None, names=col_names, skiprows=[0])
        
    

        copy_numbers = pd.read_table(copy_numbers_fname, names=["segment", "cell_label", "cn"])
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


    dtype = np.dtype([('x', int), ('y', int)])
    copy_numbers = copy_numbers.set_index(["seg_id", "cell"])
    copy_numbers= copy_numbers.unstack(level="seg_id").to_numpy(dtype)

    return Data(var, total, copy_numbers, snv_to_seg, seg_to_snvs, cell_lookup, mut_lookup)

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
    # parser.add_argument("-i", "--index", type=str, default="cell",
    #     help="index of copy number profiles")
    # parser.add_argument("-t", "--tolerance",  type=float, default=0.0,
    #     help="lower bound KL-divergence to establish a new segment")
    # parser.add_argument("-p", "--pseudo", required=False, type=float, default=1e-6,
    #     help="pseudo counts to use for computing KL- divergence")
    parser.add_argument("-D", "--data",  type=str, 
        help="filename of pickled data object")
    
    args = parser.parse_args()

    dat = load(args.file, args.profiles)
    if args.data is not None:
        dat.save(args.data)
