"""Module to prepare the input data for Pharming."""

import argparse
import pandas as pd 
import numpy as np
from data import Data


def load_from_files(var_fname, copy_fname, alpha=0.001 ):
    """Read in the data files and prepare the Pharming data object.

    Parameters
    ----------
    var_fname : str
        Filename of the variant data file. Columns should be [segment label, snv id, cell id, variant read counts, total read counts].
    copy_fname : str
        Filename of the copy number data file. Columns should be either [segment label, cell id, materinal copies, paternal copies] or
        [chromosome, start locus, stop locus, segment label, cell id, materinal copies, paternal copies].
    alpha : float, optional
        The per base sequencing error rate, by default 0.001.
    
    Returns
    -------
    Data
        A Pharming data object.
    """
    col_names = ['segment', 'mutation_label', 'cell_label','var', 'total']


    read_counts = pd.read_table(
        var_fname, header=None, names=col_names, skiprows=[0], sep=None)
    print(read_counts.head())
    
    copy_numbers = pd.read_csv(copy_fname, header=None, skiprows=[0])


    if copy_numbers.shape[1]== 7:
        copy_numbers.columns =["chr", "start", "stop", "segment", "cell_label", "x", "y"]
    elif copy_numbers.shape[1] == 4:
        copy_numbers.columns =["segment", "cell_label", "x", "y"]
    else:
        raise ValueError("Copy number file should have 4 or 7 columns")
    
    print(copy_numbers.head())

    cell_labels_not_in_df1 = copy_numbers.loc[~copy_numbers['cell_label'].isin(read_counts['cell_label']), 'cell_label']
    excluded_labels = cell_labels_not_in_df1.unique()
    print(f"Excluding {len(excluded_labels)} cells due to lack of read counts:")
    for c in excluded_labels:
        print(c)

    #exclude any cells that have no read counts at any position
    copy_numbers = copy_numbers[copy_numbers['cell_label'].isin(read_counts['cell_label'])]
  
    return load(read_counts, copy_numbers, alpha)

def load(read_counts, copy_numbers, alpha=0.001):
    """Prepare the Pharming data object from the read counts and copy number dataframes.

    Parameters
    ----------
    read_counts : pd.DataFrame
        A dataframe with the SNV variant and total read counts. Columns should be labeled [segment, snv, cell, variant, total].
    copy_numbers : pd.DataFrame
        A dataframe with the allele-specific copy numbers. Columns should either be named  [segment, cell, x, y] or
        [chr, start, stop, segment, cell, x, y].
    alpha : float, optional
        The per base sequencing error rate, by default 0.001.
    
    Returns
    -------
    Data
        A Pharming data object.
    """
    cell_labels = np.sort(read_counts['cell_label'].unique())
    mut_labels = np.sort(read_counts['mutation_label'].unique())


    #create indexed series of mapping of cell index to label
    cell_lookup = pd.Series(data=cell_labels, name="cell_label").rename_axis("cell")     
    mut_lookup = pd.Series(data=mut_labels, name="mutation_label").rename_axis("mut")

    read_counts = pd.merge(read_counts, cell_lookup.reset_index(), on='cell_label', how='left')
    read_counts = pd.merge(read_counts, mut_lookup.reset_index(), on='mutation_label', how='left')


    #get distinct segments
    if copy_numbers.shape[1] == 7:
        seg_lookup = copy_numbers.loc[:, ["segment", "chr", "start", "stop"]].drop_duplicates().sort_values(by="segment")
        seg_lookup["length"] = seg_lookup["stop"] - seg_lookup["start"] + 1
        seg_lookup["weight"] = seg_lookup["length"] / seg_lookup["length"].sum()
    elif copy_numbers.shape[1] == 4:
        seg_lookup = copy_numbers.loc[:, ["segment"]].drop_duplicates().sort_values(by="segment")
        seg_lookup["weight"] = 1.0

    seg_lookup = seg_lookup.rename_axis("seg_id").reset_index()

    seg_weights = seg_lookup["weight"].to_numpy()
    
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


    return Data(var,
                total,
                copy_x,
                copy_y,
                snv_to_seg,
                seg_to_snvs,
                cell_lookup,
                mut_lookup,
                seg_lookup=seg_lookup,
                alpha = alpha,
                seg_weights=seg_weights)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabeled columns: [chr segment snv cell var total]")
    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-m", "--mut-lookup", type=str,
        help="filename of mutation label to export")
    parser.add_argument("-l", "--cell-lookup", type=str,
        help="filename of cell label to export")
    parser.add_argument("-s", "--segment-lookup", type=str,
        help="filename of segment label to export")
    parser.add_argument("-a", "--alpha", type=float, default=0.001,
        help="sequencing error rate")
    parser.add_argument("-D", "--data",  type=str,
        help="filename of pickled data object")
    
    args = parser.parse_args()



    dat = load_from_files(args.file, args.profiles, args.alpha)
    if args.data is not None:
        dat.save(args.data)
    
    if args.mut_lookup is not None:
        dat.export_mut_lookup(args.mut_lookup)
    
    if args.cell_lookup is not None:
        dat.export_cell_lookup(args.cell_lookup)

    if args.segment_lookup is not None:
        print("saving segments...")
        dat.export_segment_lookup(args.segment_lookup)