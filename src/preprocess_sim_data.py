import argparse
import numpy as np
import pandas as pd 
import itertools 
from clonal_tree import ClonalTree
import networkx as nx
from data import load
import os 


def read_tree(fname):
    tree = nx.DiGraph()
    with open(fname, "r+") as file:
        for line in file:
            if "#leaves" in line:
                break 
            if "#edges" in line:
                continue
            edge_string = line.strip().split(" ")
            edge = [int(e) for e in edge_string]
            tree.add_edge(*edge)
    return tree


def dataframe_to_mapping(df, col="cell"):

    df = df.rename(columns={col: 'value', 'cluster':'key'})

    mapping = {key : [] for key in df['key'] }
    for index,row in df.iterrows():
        mapping[row['key']].append(row['value'])
    
    return mapping 




        
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-s", "--segments", type=str,
                        help="filename of precomputed segmentation")
    parser.add_argument( "--bin-mapping", type=str,
                        help="filename of snv to bin map")
    parser.add_argument( "--mut-mapping", type=str,
                        help="filename of snv clone assignments")
    parser.add_argument( "--cell-mapping", type=str,
                        help="filename of cell clone assignments")
    parser.add_argument( "--loss-mapping", type=str,
                        help="filename of snv loss clone assignments")
    parser.add_argument( "-t", "--tree", type=str,
                        help="filename of ground truth clonal tree")
    parser.add_argument( "-T", "--clonal-tree", type=str,
                        help="filename of pickled output clonal tree")
    parser.add_argument( "-F", "--read-counts", type=str,
                        help="filename of read counts in pharming input format")
    parser.add_argument( "-C", "--cell-cn-profiles", type=str,
                        help="filename of cell copy number profiles by segment in pharming input format")   
    parser.add_argument("-D", "--data",  type=str, 
        help="filename of pickled data object")
    parser.add_argument( "--draw",  type=str, 
        help="filename of clonal tree drawing")
    parser.add_argument("--segtrees", type=str,
        help = "directory where the segment trees should be saved")
    parser.add_argument("-L", "--likelihoods", type=str,
        help = "filename of marginal likelihood")
    # parser.add_argument("-i", "--index", type=str, default="cell",
    #     help="index of copy number profiles")
    # parser.add_argument("-t", "--tolerance",  type=float, default=0.0,
    #     help="lower bound KL-divergence to establish a new segment")
    # parser.add_argument("-p", "--pseudo", required=False, type=float, default=1e-6,
    #     help="pseudo counts to use for computing KL- divergence")

    # inpath = "/scratch/data/leah/pharming/sim_study/input/s13_n1000_m15000_c5_p0.1_l0"

    args = parser.parse_args()

    # args = parser.parse_args([  
    #     "-f", f"{inpath}/dataframe.tsv",
    #     "-c",  f"{inpath}/copy_number_profiles.csv",
    #     "-s", f"{inpath}/segementation.csv",
    #     "-t", f"{inpath}/true_tree.txt",
    #     "--bin-mapping", f"{inpath}/snv_bin_mapping.csv",
    #     "--mut-mapping", f"{inpath}/mutclust_gt.csv",
    #     "--cell-mapping", f"{inpath}/cellclust_gt.csv",
    #     "--loss-mapping", f"{inpath}/mut_loss_clust_gt.csv",
    #     "--draw", f"{inpath}/true_tree.png",
    #     "--segtree", f"{inpath}/SegTrees",
    #     "-L", f"{inpath}/segment_likelihoods.csv"
    #     ])


    







    
    segments = pd.read_csv(args.segments)




    cnames = ["chr", "mut", "cell", "base", "var", "total", "cn"]


    read_counts = pd.read_table(args.file, names=cnames)
    #assume mutation labels are unique indicies from 0 to M-1
    read_counts['chr_mutation'] = read_counts['mut']
    cell_labels = np.sort(read_counts['cell'].unique())
    mut_labels = np.sort(read_counts['chr_mutation'].unique())

 

    # cell_to_index = {c: i for i,c in enumerate(cell_labels)}
    # index_to_cell = {i: c for c,i in cell_to_index.items()}
    # mut_to_index = {c: i for i,c in enumerate(mut_labels)}
    # index_to_mut = {i: c for c,i in mut_to_index.items()}


    cell_gt = pd.read_csv(args.cell_mapping)
    # cell_gt = cell_gt.rename(columns={'cell' : 'cell_label'})
    # cell_gt['cell'] =cell_gt["cell_label"].map(cell_to_index)

    cell_mapping = dataframe_to_mapping(cell_gt, "cell")

  
  
    mut_gt  = pd.read_csv(args.mut_mapping)
    mut_mapping  = dataframe_to_mapping(mut_gt, "mutation")
    if args.loss_mapping is not None:
        loss_gt = pd.read_csv(args.loss_mapping)
        loss_mapping = dataframe_to_mapping(loss_gt, "mutation")
    else:
        loss_mapping = {}

    tree = read_tree(args.tree)
    ct = ClonalTree(0,tree, mut_mapping, loss_mapping, cell_mapping)
    if args.draw is not None:
        ct.draw(args.draw)
    
    if args.clonal_tree is not None:
        ct.save(args.clonal_tree)



    geno_cn = pd.read_csv(args.profiles)
    geno_cn["cluster"] = geno_cn["genotype"].str.replace('genotype', '').astype(int)
    geno_cn.drop("genotype", axis=1, inplace=True)
   

    geno_long =  pd.melt(geno_cn, id_vars= ['cluster'], var_name='bin', value_name='cn')
    geno_long['bin'] = geno_long["bin"].str.replace("bin", "").astype(int)
   
 

    
    bins= geno_long['bin'].unique()
    bin_to_segment = {b: segments[(b >= segments["start"]) & (b <= segments["end"])]["segment_id"].values[0] for b in bins}
    
    geno_long["segment"] = geno_long['bin'].map(bin_to_segment)

  
    geno_long = geno_long.drop("bin", axis=1).drop_duplicates()
    geno_wide =  geno_long.pivot(index='cluster', columns='segment', values='cn')
    # print(geno_wide.head())

    cn_profiles = cell_gt.merge(geno_long, on="cluster")

    # pharming input format without column headers [segment cell totalCN]"
    cn_profiles = cn_profiles[["segment", "cell", "cn"]]
    if args.cell_cn_profiles is not None:
        cn_profiles.to_csv(args.cell_cn_profiles, sep="\t", index=False, header=False)
    

  

    snv_bin_map = pd.read_csv(args.bin_mapping, names=["mut", "chr", "arm", "bin"])
    snv_bin_map["segment"] = snv_bin_map["bin"].map(bin_to_segment)
    snv_bin_map= snv_bin_map.drop(["bin", "arm"], axis=1)
   
    read_counts = read_counts.merge(snv_bin_map, on=["mut", "chr"])
    
    # pharming input format [chr segment snv cell var total]
    read_counts = read_counts[["chr", "segment", "mut", "cell", "var", "total"]]

    if args.read_counts is not None:
        read_counts.to_csv(args.read_counts, sep="\t", index=False, header=False)

    read_counts.columns =   ['chr', 'segment', 'mutation_label', 'cell_label','var', 'total']

    cn_profiles.columns = ["segment", "cell_label", "cn"]
    dat = load(read_counts, cn_profiles)

    if args.data is not None:
            dat.save(args.data)
    
       



    
  
    print("Done reading input....")

    #now we need to build the ground truth segment tree for each segment
    if args.segtrees is not None:
        if not os.path.exists(args.segtrees):
            os.makedirs(args.segtrees)
            print(f"Directory '{args.segtrees}' created successfully.\n")
        else:
            print(f"Directory '{args.segtrees}' already exists.\n")

            

        likelihoods = {}   
        for g in dat.segments:
            snvs = dat.seg_to_snvs[g]
            df= mut_gt[mut_gt['mutation'].isin(snvs)]
            mut_map= dataframe_to_mapping(df, "mutation")
            clone_map = geno_wide.iloc[:,g].to_dict()
            seg_tree = tree.copy()
            for k, cn in clone_map.items():
                seg_tree.nodes[k]["genotype"] = (cn-1, 1)
            mut_copies = {s: 1 for s in snvs}
            T_Seg= ClonalTree(g,seg_tree, mut_map, loss_mapping, cell_mapping, mut_copies)
            
            like = T_Seg.compute_likelihood(dat,g)
            likelihoods[g] =like
            print(f"Segment {g} log-likelihood: {like}")
            pred_cell, pred_mut = T_Seg.generate_results(dat.cell_lookup, dat.mut_lookup)

            # pred_cell.to_csv(f"{args.segtrees}/gt_cell_g{g}.csv", index=False)
            pred_mut = pred_mut[["mutation", "cluster"]]
            pred_mut.to_csv(f"{args.segtrees}/gt_mut_g{g}.csv",index=False)
            T_Seg.draw(f"{args.segtrees}/tree_g{g}.png")
            T_Seg.save(f"{args.segtrees}/tree_g{g}.pickle")
        
        if args.likelihoods is not None:
            print("Saving segment likelihoods...")
            with open(args.likelihoods, "w+") as file:
                file.write("segment,likelihood\n")
                for g, like in likelihoods.items():
                    file.write(f"{g},{like}\n")
     
        print("\nDone building segment trees....")
    print("Data processing complete!")
        

