import argparse
from clonal_tree import ClonalTree, load
from data import Data 
from copy import deepcopy
from cell_mapping import CellAssign
import pandas as pd 
from fit_segment import FitSegmentTree, draw 
import networkx as nx 
from sklearn.metrics.cluster import adjusted_rand_score
import logging

from sti_v2 import STI, pickle_object


import multiprocessing  
import numpy as np
from scipy.stats import binom

def dict_to_df(mydict, colnames):
    df = pd.DataFrame.from_dict(mydict, orient='index', columns=[colnames[1]])

# Reset index to create a column from dictionary keys
    df.reset_index(inplace=True)
    df.columns = colnames
    return df 


def compute_ari(alpha1, alpha2) -> float:
        #  gt_mut = self.get_mut_clusters()
         gt =   pd.Series(alpha1)
         pred= pd.Series(alpha2)
         df = pd.concat([gt, pred], axis=1, keys=['gt', 'pred'])
        #  pred_mut = obj.get_mut_clusters()
 
         return adjusted_rand_score(df["gt"].values, df["pred"].values)



def sti_fit(seg, gt, T_m, phi, delta, lamb=1e3, lamb2=1e5):


        gt_seg = deepcopy(gt)
        gt_psi = gt_seg.get_psi()

        snvs = dat.seg_to_snvs[seg]
        gt_seg.filter_snvs(snvs)
        cost = gt_seg.compute_likelihood(dat, phi, lamb)
        print(f"Segment {seg} ground truth cost: {cost} ")
        # gt_seg.draw(f"test/gt_seg{seg}.png", phi, segments=[seg] )
        gt_psi = gt.get_psi()
        
        # pd.Series(gt_psi).to_csv(f"test/s11_psi.csv")
        cn_states, counts = dat.cn_states_by_seg(seg)
        print(f"{seg}: {cn_states}")
        S = gt.get_cna_tree(seg)

        st  = STI(S,T_m, delta, lamb1=lamb)


        trees = st.fit(dat, seg)
        print(f"Segment {seg} ground truth cost: {cost}  best inferred cost: {trees[0].cost}")
        score_results = []
        for i,res in enumerate(trees):
            tree = res.ct
            ca = res.phi
            scores = gt_seg.score_snvs(tree)
            scores["tree"] = i
            scores["segment"] = seg 
            scores["num_states"] = len(cn_states)
            scores["cell_ari"] =phi.compute_ari(ca)
            scores["cost"] = res.cost
            scores["gt_cost"] = cost
            score_results.append(scores)


        pickle_object(trees[:2], f"test/seg{seg}_trees.pkl")
        return trees, score_results



if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-s" ,"--seed", required=False, type=int,
                        help="random number seed (default: 1026)")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=10,
                        help="lambda value, default=1")
    parser.add_argument("-j" ,"--cores", required=False, type=int,default=1,
                        help="number of processes")
    parser.add_argument("-g" ,"--segment", required=False, type=int,
                        help="segment id of tree to build")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory where output files should be written")
    parser.add_argument("-S", "--scores", type=str,
        help = "scores for segment trees")
    # parser.add_argument("--state-trees", required=False, 
    #                     help= "filename of pregenerated state trees or path of program to generate state trees" )
    
    args = parser.parse_args()
    # parser.add_argument('-g', '--segment', type=int, required=False)
    # parser.add_argument("-d", "--data", type=str)

    instance = "s11_m5000_k25_l7"
    # instance = "s12_m5000_k25_l7"
    folder = "n1000_c0.25_e0" 
    pth = f"simulation_study/input"

    args = parser.parse_args([

        "-d", f"{pth}/{instance}/{folder}/data.pkl",
        "-t", f"{pth}/{instance}/gt.pkl",
        "-c", f"{pth}/{instance}/{folder}/cellAssign.pkl",
        "-j", "4",
        # "-s", "14",
        # "--segment", "0",
        # "--out", f"/Users/leah/Documents/Research/projects/Pharming/test",
        "-S", f"/Users/leah/Documents/Research/projects/Pharming/test/scores_c0.25.csv",

    ])



# S = nx.DiGraph([((1,1), (2,1)), ((1,1), (4,1))])
    from data import Data, load_from_pickle

    lamb = args.lamb
    dat = load_from_pickle(args.data)
    print(dat)
    gt = load_from_pickle(args.tree)
    print(gt)
 
    phi = load(args.cell_assign)
    print(phi)


    
    gt.reindex_snvs(dat.mut_lookup)
    # gt.draw("test/gt_tree.png")

    gt_psi = gt.get_psi()
    psi_df  = dict_to_df(gt_psi, ["snv", "node"])
    # psi_df.to_csv(f"test/psi.csv", index=False)
    phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    root = gt.root

    gt_T_m = gt.mutation_cluster_tree()
    gt_T_m.remove_node(root)


    mapping = {}
    nodes = list(gt_T_m.nodes)
    nodes.sort()
    for i,n in enumerate(nodes):
          mapping[n] =i 
    
    T_m = nx.relabel_nodes(gt_T_m, mapping)

    gt_delta = {mapping[n]: gt_dcfs[n] for n in mapping if n != root}

    all_scores = []
    for seg in  [10,20]:#dat.seg_to_snvs:
        cn_states, counts = dat.cn_states_by_seg(seg)
        if len(cn_states) <= 1:
              continue
        print(f"Inferring segment {seg}....")
        try:
            trees, score_results = sti_fit(seg, gt, deepcopy(T_m), phi, gt_delta.copy())
            all_scores.append(score_results)
        except Exception as Argument: 
            logging.exception(f"Error on segment: {seg} ")
        # if len(trees) <= 3:
        #     num_sol = len(trees)
        # else:
        #     num_sol = 3
        # for i in range(num_sol):
        #         trees[i].png(f"test/seg{seg}_tree{i}.png", segments=[seg])
    
    flattened_list = [item for sublist in all_scores for item in sublist]
    pd.DataFrame(flattened_list).to_csv(args.scores, index=False)
    print("done")


    
#     df = dict_to_df(gt_dcfs, ["node", "dcf"])
# #     df = pd.DataFrame.from_dict(gt_dcfs, orient='index', columns=['dcfs'])

# # # Reset index to create a column from dictionary keys
# #     df.reset_index(inplace=True)
# #     df.columns = ['nodes', 'dcfs']
#     df.to_csv("test/dcfs.csv", index=False)
#     # dcfs, _ = gt.compute_dcfs(phi, 20)
#     # for k,delta in dcfs.items():
#     #       print(f"{k}: {dcfs}")
#     cna_genos = gt.get_cna_genos()
#     cost = gt.compute_costs(dat, phi,lamb )
#     # gt.draw(f"{args.out}/gt.png")

#     score_results = []

#     best_trees = {}
#     # skip_segs = [10,12,17,20]
#     skip_segs= []
#     segments_to_process = list(dat.seg_to_snvs.keys())
#     # segments_to_process = [12]  #12

#     res = sti_fit(20, gt,dat, phi)
#     res[0].ct.draw("test/best_tree.png", res[0].phi)


#     # Create a pool of processes
#     # with multiprocessing.Pool(args.cores) as pool:
#     #         # Use the pool to parallelize the processing of segments
#     #     score_results = pool.starmap(infer_segment_tree, [(seg, gt, dat, phi, lamb) for seg in segments_to_process])

#     merged_list = [item for sublist in score_results for item in sublist]
#     results = pd.DataFrame(merged_list)
#     results.to_csv(f"{args.scores}", index=False)
#     # print(f"Skipped segments: {skip_segs}")
#     print('done')
    # for seg in dat.seg_to_snvs:
    #     # if seg == 20:
    #     # if seg in [1, 20,23]:
    #         cna_genos = gt.get_cna_genos()
    #         gt_seg = deepcopy(gt)
    #         snvs = dat.seg_to_snvs[seg]
    #         gt_seg.filter_snvs(snvs)
    #         gt_cost = gt_seg.compute_costs(dat, phi, lamb)
    #         print(phi)
    #         print(f"GT cost for seg {seg}: {gt_cost}")
    #         gt_seg.draw(f"{args.out}/gt_segtree{seg}.png", segments=[seg])
    #         gt_psi = gt.get_psi()

    #         pd.Series(gt_psi).to_csv(f"{args.out}/psi.csv")
    #         dcfs, dcfs_by_state = gt.compute_dcfs(phi, seg)
    #         print(dcfs)
 

        
        
    #         cn_states, counts = dat.cn_states_by_seg(seg)
    #         print(f"{seg}: {cn_states}")
    #         # continue
 
    #         if len(cn_states) <=3:
            
    #             S = gt.get_cna_tree(seg)
    #             draw(S, f"{args.out}/cna_tree_seg{seg}.png")
    #             gt_psi = gt_seg.get_psi()
    #             fs = FitSegmentTree(S, max_clusters=2, psi=gt_psi)
            
    #             try:
    #                 best_tree =fs.fit(seg, dat, lamb, top_n=5)
            
            
    #                 if len(best_tree) >0:
    #                     for i, tup in enumerate(best_tree):
    #                         tree, ca = tup
    #                     # tree = best_trees[seg]
    #                         # ca, cost, _, _ = tree.assign_cells(dat, lamb)
    #                         cost =tree.compute_costs(dat,ca, lamb)
    #                         tree.draw(f"{args.out}/seg{seg}_tree{i}.png", segments=[seg])

    #                         cell_count = ca.get_cell_count()
                 
    #                         for u,n  in cell_count.items():
    #                             if n ==0 and tree.is_leaf(u):
    #                                 print(f"{u} invalid cell assigment")
    #                         pd.Series(tree.get_psi()).to_csv(f"{args.out}/tree_psi.csv")
    #                         scores = gt_seg.score_snvs(tree)
    #                         scores["tree"] = i
    #                         scores["segment"] = seg 
    #                         scores["num_states"] = len(cn_states)
    #                         scores["cell_ari"] =phi.compute_ari(ca)
    #                         score_results.append(scores)

    #             except:
    #                 skip_segs.append(seg)
    #                 print(f"Error building tree for segment {seg}, continuing... ")
 