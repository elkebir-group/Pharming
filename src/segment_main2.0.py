import argparse
from clonal_tree import ClonalTree, load
from data import Data 
from copy import deepcopy
from cell_mapping import CellAssign
import pandas as pd 
import networkx as nx 
from sklearn.metrics.cluster import adjusted_rand_score
import logging
import utils
import itertools
from cna_merge import CNA_Merge

from sti_v2 import STI



import numpy as np


def score_tree(gt, gt_phi, inf,inf_phi, segments):
      
        scores = gt.score_snvs(inf)
        scores["tree"] = i
        scores["segment"] = ":".join([str(ell) for ell in segments])
        scores["cell_ari"] =gt_phi.compute_ari(inf_phi)
        scores["gt_cost"] = gt.cost 
        scores["inf_cost"] = inf.cost 
        
        return scores 
            

     

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
            scores = score_tree(gt_seg, phi, tree, ca, segments=[seg])
            scores["num_states"] = len(cn_states)
            score_results.append(scores)

       
          


        utils.pickle_object(trees[:2], f"test/seg{seg}_trees.pkl")
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
        "-S", f"/Users/leah/Documents/Research/projects/Pharming/test/cna_scores.csv",

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

    with open("test/input/T_m.txt", "w+") as file:
        for u,v in T_m.edges:
              file.write(f"{u}\t{v}\n")

    utils.pickle_object(T_m, "test/T_m.pkl")
    utils.draw(T_m, "test/T_m.png")

    gt_delta = {mapping[n]: gt_dcfs[n] for n in mapping if n != root}

    with open("test/input/dcfs.txt", "w+") as file:
        for key,val in gt_delta.items():

              file.write(f"{val}\n")

    test_segs = [0,10,20, 24]
    snvs = list(itertools.chain(*[dat.seg_to_snvs[seg] for seg in test_segs]))
    gt.filter_snvs(snvs)
    cost = gt.compute_likelihood(dat, phi, lamb)
    gt.draw("test/output/gt_tree.png", phi, segments=test_segs)

    # all_scores = []
    # tree_sols = {}
    # for seg in  test_segs:#dat.seg_to_snvs:
    #     cn_states, counts = dat.cn_states_by_seg(seg)
    #     if len(cn_states) <= 1:
    #           continue
    #     print(f"Inferring segment {seg}....")
    #     try:
    #         trees, score_results = sti_fit(seg, gt, deepcopy(T_m), phi, gt_delta.copy())

    #         all_scores.append(score_results)
    #         tree_sols[seg] = trees 
    #     except Exception as Argument: 
    #         logging.exception(f"Error on segment: {seg} ")
    #     # if len(trees) <= 3:
    #     #     num_sol = len(trees)
    #     # else:
    #     #     num_sol = 3
    #     # for i in range(num_sol):
    #     #         trees[i].png(f"test/seg{seg}_tree{i}.png", segments=[seg])
    
    # flattened_list = [item for sublist in all_scores for item in sublist]
    # pd.DataFrame(flattened_list).to_csv(args.scores, index=False)
    # utils.pickle_object(tree_sols, "test/tree_sols.pkl")
    tree_sols = utils.load_pickled_object("test/tree_sols.pkl")
    merge_list = []
    root = [n for n in T_m if T_m.in_degree[n]==0][0]
    T_m.add_edge(max(T_m)+2, root) 
    bad_trees = []
    for l1, l2 in  itertools.combinations(tree_sols.keys(), 2):
         for sol1 in tree_sols[l1]:
              CT1 = sol1.ct
              for sol2 in tree_sols[l2]:
                CT2 = sol2.ct 
                try:
                    all_sols = CNA_Merge(CT1, CT2, T_m).fit(dat, lamb=1e5)
                    merge_list += all_sols
                except:
                    bad_trees.append((CT1, CT2))
                    all_sols = CNA_Merge(CT1, CT2, T_m, verbose=True).fit(dat, lamb=lamb)

    # utils.pickle_object(bad_trees, "test/bad_trees.png")            
    # print(len(bad_trees))
    merge_list = sorted(merge_list, key=lambda x: x[0])
    for i in range(min(len(merge_list),5)):
        obj, ca, ct = merge_list[i]
        ct.draw(f"best_tree{i}.png", ca, segments=test_segs)

   
    score_results= [score_tree(gt, phi, ct, ca, segments=test_segs) for obj,ca, ct in merge_list]

    pd.DataFrame(score_results).to_csv(args.scores, index=False)
  
                     


              


    print("done")

