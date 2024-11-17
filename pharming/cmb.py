

# Created by: L.L. Weber
# Created on: 2024-05-06 09:20:11



import pandas as pd 
import numpy as np 
import networkx as nx 
import itertools
import argparse 
from collections import defaultdict



def vaf_validation(ct, ca, dat, min_cells=10):
    """Compute the observed vaf for each snv in the tree and compare it to the latent vaf"""
    vaf_dict = defaultdict(list)
    for n in ct.tree:
        cells = ca.get_cells(n)
        if len(cells) ==0:
            continue
        
        genos =ct.genotypes[n]
        for j, g in genos.items():
            vaf = (g[2] + g[3])/(g[0] + g[1])
            vaf_dict[vaf, j] += list(cells )
    results= []
    for vaf,j in vaf_dict:
        cells = list(set(vaf_dict[vaf,j]))
        if len(cells) > min_cells:
            var = dat.var[cells,:][:,j].sum(axis=0)
            total = dat.total[cells,:][:,j].sum(axis=0)
            obs_vaf = var/total
            results.append([j,vaf, obs_vaf, var, total])
    df = pd.DataFrame(results, columns=["snv", "latent_vaf", "obs_vaf", "var", "total"])
    print(df.head())
    return df 

    


def loss_clade_cmb(n, ca, ct, dat, min_cells=10):
    if len(ct.mut_loss_mapping[n]) > 0:
        snvs = ct.mut_loss_mapping[n]
        clade_nodes =  nx.descendants(ct.tree, n) | {n}

        within_clade_cells = list(itertools.chain.from_iterable([ca.get_cells(u) for u in clade_nodes]))

        if len(within_clade_cells) >= min_cells:
            within_cmb = dat.compute_cmb(within_clade_cells, snvs)
            df = pd.DataFrame({'cell': within_clade_cells, 'cmb': within_cmb})
            df["within_clade"] = 1

        else:
            df = pd.DataFrame(columns=['cell', 'cmb', 'within_clade'])
        df["clade_snvs"] = len(snvs)
        df["clade"] = n

        return df
    else:
        df = pd.DataFrame(columns=['cell', 'cmb', 'within_clade'])
        df["clade_snvs"] = 0
        df["clade"] = n

   



def clade_cmb(n, ca, ct, dat, min_cells=10):
   
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

    outside_df = None 
    within_df = None
    #returns a vector of length outside_clade_cells
    if len(outside_clade_cells) >= min_cells:
        outside_cmb = dat.compute_cmb(outside_clade_cells, not_lost_snvs)
        outside_df = pd.DataFrame({'cell': outside_clade_cells, 'cmb': outside_cmb})
        outside_df["within_clade"] = 0

    if len(within_clade_cells) >= min_cells:
        within_cmb = dat.compute_cmb(within_clade_cells, not_lost_snvs)
        within_df = pd.DataFrame({'cell': within_clade_cells, 'cmb': within_cmb})
        within_df["within_clade"] = 1

    if outside_df is not None and within_df is not None:
        df = pd.concat([outside_df, within_df])
    elif outside_df is not None:
        df = outside_df
    elif within_df is not None:
        df = within_df
    else:
         df = pd.DataFrame(columns=['cell', 'cmb', 'within_clade'])
    df["clade_snvs"] = len(not_lost_snvs)
    df["clade"] = n

    return df 

   
def cmb_all(ct, ca, dat, min_cells=10):
    ct.update_mappings()
    cmb_list = []
    for n in ct.mut_mapping:
        if len(ct.mut_mapping[n]) > 0:
            cmb_list.append(clade_cmb(n, ca, ct, dat, min_cells=min_cells))

    cmb_df = pd.concat(cmb_list)
    cmb_df["cell"] = cmb_df["cell"].astype(int)
  
    return cmb_df 

def cmb_loss_all(ct, ca, dat, min_cells=10):
    ct.update_mappings()
    cmb_list = []
    for n in ct.mut_mapping:
            cmb_list.append(loss_clade_cmb(n, ca, ct, dat, min_cells=min_cells))

    cmb_df = pd.concat(cmb_list)
    cmb_df["cell"] = cmb_df["cell"].astype(int)
  
    return cmb_df 

def main(args):
    sol = pd.read_pickle(args.solution)
    ct = sol[0].ct
    ca = sol[0].phi
    dat  = pd.read_pickle(args.data)

    if args.vaf is not None:
        vaf_df = vaf_validation(ct, ca, dat, min_cells = args.min_cells)
        vaf_df.to_csv(args.vaf, index=False)

    if args.out:
        cmb_df = cmb_all(ct, ca, dat, args.min_cells)
        cmb_df.to_csv(args.out, index=False)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-s", "--solution", required=False,
                        help = "pickled solution file")
    parser.add_argument("--min-cells",  required=False, type=int, default=10,
                        help = "minimum number of cells to compute the scores for a clade")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="filename where cmb values will be saved")
    parser.add_argument("-v" ,"--vaf", required=False, type=str,
                        help="filename where vaf values will be saved")


    args = parser.parse_args()
    # pth = "/scratch/data/leah/Pharming/dlp/SA921"
    # pth2 = "pharming/decifer/isegs3_tm4_top5_lamb1000"
    # args = parser.parse_args([
    #     "-d", f"{pth}/input/data.pkl",
    #     "-s", f"dlp/SA921/pharming/decifer_k5/s20_r25/isegs10_tm10_top5_lamb1000/solutions.pkl",
    #     "-o", f"dlp/SA921/pharming/decifer_k5/s20_r25/isegs10_tm10_top5_lamb1000/cmb.csv",  
    # ])

    # args = parser.parse_args()
    # pth = "act/TN3"
    # pth2 = "pharming/decifer_k6/isegs10_tm10_top5_lamb1000"
    # args = parser.parse_args([
    #     "-d", f"{pth}/input/data.pkl",
    #     "-s", f"{pth}/{pth2}/solutions.pkl",
    #     "-o", "test/TN3_k6.csv"
    # ])



    # tpath = "/scratch/data/leah/Pharming/simulation_study/sims/s10_m10000_k25_l5_d2/n1000_c0.25_e0"
    # gt_fname = f"{tpath}/gt.pkl"
    # ca_fname = f"{tpath}/phi.pkl"
    # dat_fname = f"{tpath}/data.pkl"
    main(args)