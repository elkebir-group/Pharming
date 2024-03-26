import argparse
from clonal_tree import ClonalTree, load
from cell_mapping import CellAssign
import pandas as pd 
from sklearn.metrics.cluster import adjusted_rand_score
from utils import load_pickled_object, draw
import itertools
import numpy as np 

from data import Data
from copy import deepcopy



def score_tree(gt, gt_phi, inf,inf_phi, segments=None):
        if segments is None:
              segments = inf.get_segments()
        gt.filter_snvs(inf.get_all_muts())
        _ = gt.compute_likelihood(dat, gt_phi, lamb)
        scores = gt.score_snvs(inf)
        # scores["tree"] = i
        scores["segment"] = ":".join([str(ell) for ell in segments])
        scores["cell_ari"] =gt_phi.compute_ari(inf_phi)
        scores["gt_cost"] = gt.cost 
        scores["gt_snv_cost"] =gt.snv_cost
        scores["gt_cna_cost"] = gt.cna_cost
        scores["inf_cost"] = inf.cost 
        scores["inf_snv_cost"] =inf.snv_cost
        scores["inf_cna_cost"] = inf.cna_cost
        scores['cna_mad'] = gt.cna_genotype_similarity(gt_phi, inf, inf_phi)

        cna_trees_correct = [compare_CNA_trees(gt,inf, ell ) for ell in segments]
        scores["perc_cna_trees_correct"] = sum(cna_trees_correct)/len(segments)
              
        
        return scores 
            

def compare_CNA_trees(gt, inf, segment):
        S_gt = gt.get_cna_tree(segment)
        # draw(S_gt,f"test/S_{segment}_gt.png")
        S_inf= inf.get_cna_tree(segment)
        # draw(S_inf,f"test/S_{segment}_inf.png")
        return set(S_gt.edges) == set(S_inf.edges)


def eval_segtree(gt, gt_phi, sol, segment, dat, lamb=1e3):
    snvs = sol.ct.seg_to_muts[segment]
    # sol.ct.filter_snvs(snvs)
    # sol.ct.compute_likelihood(dat, sol.phi, lamb)
    gt.filter_snvs(snvs)

    # gt.prune(gt_phi)
    _ = gt.compute_likelihood(dat, gt_phi, lamb)
    gt.draw(f"test/seg{segment}_gt.png", gt_phi, segments=[segment], include_dcfs=True)




    scores  = score_tree(gt, gt_phi, sol.ct, sol.phi, segments=[segment])
    scores["nsnvs"] = len(snvs)

    return scores 



def compute_ari(alpha1, alpha2) -> float:
        #  gt_mut = self.get_mut_clusters()
         gt =   pd.Series(alpha1)
         pred= pd.Series(alpha2)
         df = pd.concat([gt, pred], axis=1, keys=['gt', 'pred'])
        #  pred_mut = obj.get_mut_clusters()
 
         return adjusted_rand_score(df["gt"].values, df["pred"].values)




def save_psi(gt, inf, fname):
      gt_psi = pd.Series(gt.get_psi())
      inf_psi = pd.Series(inf.get_psi())
      df = pd.concat([gt_psi, inf_psi],axis=1)
      df.reset_index(inplace=True)
      df.columns = ["snv", "gt_psi", "inf_psi"]
      print(df.head())
      df.to_csv(fname, index=False)

def seg_score_wrapper(segtrees, ell):
    seg_scores = []
    segtrees[0].png(f"test/seg_inf{ell}.png")
    for sol in segtrees:

          seg_scores.append(eval_segtree(deepcopy(gt), phi, sol, ell, dat, lamb))
    return seg_scores

if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-S", "--solutions", type=str,
                        help = "scores for segment trees")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=10,
                        help="lambda value, default=1")
    parser.add_argument("-g" ,"--segment", required=False, type=int,
                        help="segment id of tree to build")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory where score files should be written")


    
    args = parser.parse_args()


    instance = "s11_m5000_k25_l5"
    # instance = "s12_m5000_k25_l7"
    folder = "n1000_c0.05_e0" 
    pth = f"simulation_study/input"

    args = parser.parse_args([

        "-d", f"{pth}/{instance}/{folder}/data.pkl",
        "-t", f"{pth}/{instance}/{folder}/gt.pkl",
        "-c", f"{pth}/{instance}/{folder}/phi.pkl",
        # "-s", "14",
        # "--segment", "0",
        "-o", "test/scores.csv",
        "-S", f"test/solutions.pkl",

    ])







    lamb = args.lamb
    dat = load_pickled_object(args.data)



    # if args.segment is None:
    #     eval_segs = [ell for ell in dat.segments if dat.num_cn_states(ell) > 1]
    # else:
    #     eval_segs = args.segment 

    # snvs = list(itertools.chain(*[dat.seg_to_snvs[seg] for seg in eval_segs]))

    gt = load_pickled_object(args.tree)
    phi = load(args.cell_assign)

    # gt.reindex_snvs(dat.mut_lookup)

    phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    root = gt.root

    # gt.filter_snvs(snvs)
    Tm = gt.mutation_cluster_tree()
    
    draw(Tm, "test/Tm.png")
    cost = gt.compute_likelihood(dat, phi, lamb)
    seg2 = deepcopy(gt)
    seg2.filter_segments([2])
    cost2 = seg2.compute_likelihood(dat, phi, lamb)
    seg2.draw("test/gt_seg2.png", phi, segments=[2], include_dcfs=True)

    sol_list =   load_pickled_object(args.solutions)
    for i,sol in enumerate(sol_list):
          sol.png(f"test/ct{i}_seg2_0.png")
    score_results= [score_tree(deepcopy(gt), phi, sol.ct, sol.phi) for sol in sol_list]
    pd.DataFrame(score_results).to_csv(args.out, index=False)

    # mysegs = [0,1,2,11,13,14,18,24]
    # all_scores = []
    # for ell in mysegs:
    #       segtrees = load_pickled_object(f"test/segtrees{ell}.pkl")
    #       all_scores += seg_score_wrapper(segtrees, ell)

    # pd.DataFrame(all_scores).to_csv("test/segtree_scores.csv", index=False)    


    # seg16map = {8:8, 4:3, 3:3, 0:0, 1:1, 2:2, 6:5, 7:6, 5:4 }
    # from copy import deepcopy
    # phi16 = deepcopy(phi)
    # phi16.relabel_clones(seg16map)
    # from utils import pickle_object
    # pickle_object(phi16, "test/phi16.pkl")


    


    # sol_list  = load_pickled_object(args.solutions)
    # results= []
    # sol_dict =  load_pickled_object("test/segtree_solutions.pkl")
    # mysegs = [16, 19, 12, 10, 8, 17]
    # for ell in mysegs:
    #     for i, sol in enumerate(sol_dict[ell]):
    #         sol.png(f"test/seg{ell}_{i}.png")

    # for ell, sol_list in sol_dict.items():
    #       num_cn_states = dat.num_cn_states(ell)
    #       gt_seg = deepcopy(gt)
    #       for sol in sol_list:
    #         scores = eval_segtree(gt_seg, phi, sol,ell, dat, lamb )
    #         scores["num_cn_states"] = num_cn_states 
    #         results.append(scores)
    
    # pd.DataFrame(results).to_csv("test/segtree_scores.csv", index=False)      


    # ell = 16
    # num_cn_states = dat.num_cn_states(ell)
    # results= []
    # gt_seg = deepcopy(gt)
    # seg_scores = []
    # gt.filter_snvs(dat.seg_to_snvs[ell])
    # save_psi(gt, sol_list[0].ct, "test/psi_comp.csv")
    # for sol in sol_list:
    #     scores = eval_segtree(gt_seg, phi, sol,ell, dat, lamb )
    #     scores["num_cn_states"] = num_cn_states 
    #     results.append(scores)
    # pd.DataFrame(results).to_csv("test/segtree16_scores.csv", index=False)  

    # segtrees  = load_pickled_object("test/segrees_ilp.pkl")
    # seg_scores = []
    # for seglist in segtrees:
    #       sol = seglist[0]
    #       seg  = list(sol.segments)[0]
    #       seg_scores.append(eval_segtree(deepcopy(gt), phi, sol, seg, dat, lamb))
    # pd.DataFrame(seg_scores).to_csv("test/seg_scores_ilp.csv", index=False)
    # print("done")

    # for i,best_sol in enumerate(sol_list):
 

    #         best_sol.png(f"test/inf{i}.png")

    #     # for ell in best_sol.ct.get_segments():
    #         scores.append(eval_segtree(deepcopy(gt), phi, deepcopy(best_sol), ell, dat, lamb))
    # pd.DataFrame(seg_scores).to_csv("test/seg_scores.csv", index=False)

          



    print("done")
    # sol_list[0].png("test/opt_tree.png")


    # score_results = []
    # clonal_trees  = load_pickled_object("test/clonal_trees.pkl")
    # for i,sol_list in enumerate(clonal_trees):
    #     sol_list = sorted(sol_list, key= lambda x: x.cost )
    #     sol_list[0].png(f"test/opt_tree_tm{i}.png")
    #     score_results += [score_tree(gt, phi, sol.ct, sol.phi, segments=eval_segs) for sol in sol_list]
    # pd.DataFrame(score_results).to_csv(args.out, index=False)



    