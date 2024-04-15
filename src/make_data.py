import pandas as pd 
import numpy as np
import argparse
from data import Data


def load_from_files(var_fname, copy_fname ):
    col_names = ['segment', 'mutation_label', 'cell_label','var', 'total']


    read_counts = pd.read_table(
        var_fname, header=None, names=col_names, skiprows=[0], sep=None)
    print(read_counts.head())
    
    copy_numbers = pd.read_csv(copy_fname, header=None,names=["segment", "cell_label", "x", "y"], skiprows=[0])
    print(copy_numbers.head())

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

# def segment(cn_profiles, tol=0.0, pseudo_counts=1e-6):
#     '''
#     cn_profiles: an numpy array with shape N x B (cell by bin) matrix with the total copy number profiles of for each cell in each bin
    
#     returns: a pd.DataFrame with Segment ID as index and columns start and end bin (inclusive)
#     '''
#     return Segment(tol=tol, pseudo_counts=pseudo_counts).fit(cn_profiles)



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


    # instance = "s11_m5000_k25_l5_d2"
    # # instance = "s12_m5000_k25_l7"
    # folder = "n2000_c0.25_e0" 
    # pth = f"simulation_study/sims"

    # opth = "test"

    # args = parser.parse_args([

    #     "-f", f"{pth}/{instance}/{folder}/sparse.p0",
    #     "-c", f"{pth}/{instance}/{folder}/cells.p0",
    #     "-D", f"{opth}/data.pkl",

    # ])



    dat = load_from_files(args.file, args.profiles)
    if args.data is not None:
        dat.save(args.data)
    
    if args.mut_lookup is not None:
        dat.export_mut_lookup(args.mut_lookup)
    
    if args.cell_lookup is not None:
        dat.export_cell_lookup(args.cell_lookup)
