

# Created by: L.L. Weber
# Created on: 2024-05-06 09:20:11


from clonal_tree import ClonalTree
import utils 
import pandas as pd 
import numpy as np 
import networkx as nx 
import itertools
import argparse 






def clade_cmb(n, ca, ct, dat):
   
    clade_nodes =  nx.descendants(ct.tree, n) | {n}
    outside_nodes =  ct.clones() - clade_nodes 

    snvs = ct.mut_mapping[n]

    within_clade_cells = list(itertools.chain.from_iterable([ca.get_cells(u) for u in clade_nodes]))
    outside_clade_cells =list(itertools.chain.from_iterable([ca.get_cells(u) for u in outside_nodes]))

    not_lost_snvs = []
    for j in snvs:
        not_lost = []
        for u in clade_nodes:
            geno = ct.genotypes[u][j]
            not_lost.append(max(geno[2], geno[3]) > 0)
        if all(not_lost):
            not_lost_snvs.append(j)


    #returns a vector of length outside_clade_cells
    outside_cmb = dat.compute_cmb(outside_clade_cells, not_lost_snvs)
    within_cmb = dat.compute_cmb(within_clade_cells, not_lost_snvs)


    outside_df = pd.DataFrame({'cell': outside_clade_cells, 'cmb': outside_cmb})
    outside_df["within_clade"] = 0

    within_df = pd.DataFrame({'cell': within_clade_cells, 'cmb': within_cmb})
    within_df["within_clade"] = 1

    df = pd.concat([outside_df, within_df])
    df["clade_snvs"] = len(not_lost_snvs)
    df["clade"] = n

    return df 

   
def cmb_all(ct, ca, dat):
    ct.update_mappings()
    cmb_list = []
    for n in ct.mut_mapping:
        if len(ct.mut_mapping[n]) > 0:
            cmb_list.append(clade_cmb(n, ca, ct, dat))

    cmb_df = pd.concat(cmb_list)
    cmb_df["cell"] = cmb_df["cell"].astype(int)
  
    return cmb_df 

def main(args):
    sol = pd.read_pickle(args.solution)
    ct = sol[0].ct
    ca = sol[0].phi
    dat  = pd.read_pickle(args.data)

    cmb_df = cmb_all(ct, ca, dat)
    cmb_df.to_csv(args.out, index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-s", "--solution", required=False,
                        help = "pickled solution file")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="filename where cmb values will be saved")

    pth = "/scratch/data/leah/Pharming/dlp/SA1090_v0"
    pth2 = "pharming/decifer/isegs3_tm4_top5_lamb1000"
    args = parser.parse_args()
    
    args = parser.parse_args([
        "-d", f"{pth}/input/data.pkl",
        "-s", f"{pth}/{pth2}/solutions.pkl",
        "-o", f"{pth}/{pth2}/cmb.csv"
    ])


    # tpath = "/scratch/data/leah/Pharming/simulation_study/sims/s10_m10000_k25_l5_d2/n1000_c0.25_e0"
    # gt_fname = f"{tpath}/gt.pkl"
    # ca_fname = f"{tpath}/phi.pkl"
    # dat_fname = f"{tpath}/data.pkl"
    main(args)