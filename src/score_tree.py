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



def score_tree(gt, gt_phi, inf,inf_phi, segments):
      
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
    sol.ct.filter_snvs(snvs)
    sol.ct.compute_likelihood(dat, sol.phi, lamb)
    gt.filter_snvs(snvs)
    gt.prune(gt_phi)
    cost = gt.compute_likelihood(dat, gt_phi, lamb)
    # gt.draw(f"test/seg{segment}_gt.png", gt_phi, segments=[segment], include_dcfs=True)




    scores  = score_tree(gt, gt_phi, sol.ct, sol.phi, segments=[segment])

    return scores 



def compute_ari(alpha1, alpha2) -> float:
        #  gt_mut = self.get_mut_clusters()
         gt =   pd.Series(alpha1)
         pred= pd.Series(alpha2)
         df = pd.concat([gt, pred], axis=1, keys=['gt', 'pred'])
        #  pred_mut = obj.get_mut_clusters()
 
         return adjusted_rand_score(df["gt"].values, df["pred"].values)







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


    instance = "s11_m5000_k25_l7"
    # instance = "s12_m5000_k25_l7"
    folder = "n1000_c0.05_e0" 
    pth = f"simulation_study/input"

    args = parser.parse_args([

        "-d", f"{pth}/{instance}/{folder}/data.pkl",
        "-t", f"{pth}/{instance}/gt.pkl",
        "-c", f"{pth}/{instance}/{folder}/cellAssign.pkl",
        # "-s", "14",
        # "--segment", "0",
        "-o", f"/Users/leah/Documents/Research/projects/Pharming/test/tree_scores.csv",
        "-S", f"test/solution.pkl",

    ])







    lamb = args.lamb
    dat = load_pickled_object(args.data)


    if args.segment is None:
        eval_segs = [ell for ell in dat.segments if dat.num_cn_states(ell) > 1]
    else:
        eval_segs = args.segment 

    snvs = list(itertools.chain(*[dat.seg_to_snvs[seg] for seg in eval_segs]))

    gt = load_pickled_object(args.tree)
    phi = load(args.cell_assign)

    gt.reindex_snvs(dat.mut_lookup)

    phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    root = gt.root

    gt.filter_snvs(snvs)

    cost = gt.compute_likelihood(dat, phi, lamb)


    sol_list  = load_pickled_object(args.solutions)

    # segtrees  = load_pickled_object("test/segrees_ilp.pkl")
    # seg_scores = []
    # for seglist in segtrees:
    #       sol = seglist[0]
    #       seg  = list(sol.segments)[0]
    #       seg_scores.append(eval_segtree(deepcopy(gt), phi, sol, seg, dat, lamb))
    # pd.DataFrame(seg_scores).to_csv("test/seg_scores_ilp.csv", index=False)
    # print("done")

    seg_scores = []
    for i,best_sol in enumerate(sol_list):
            ell=2

            best_sol.png(f"test/seg{ell}_inf{i}.png")
        # for ell in best_sol.ct.get_segments():
            seg_scores.append(eval_segtree(deepcopy(gt), phi, deepcopy(best_sol), ell, dat, lamb))
    pd.DataFrame(seg_scores).to_csv("test/seg_scores.csv", index=False)

          

    score_results= [score_tree(gt, phi, sol.ct, sol.phi, segments=eval_segs) for sol in sol_list]
    pd.DataFrame(score_results).to_csv(args.out, index=False)

    print("done")
    # sol_list[0].png("test/opt_tree.png")


    # score_results = []
    # clonal_trees  = load_pickled_object("test/clonal_trees.pkl")
    # for i,sol_list in enumerate(clonal_trees):
    #     sol_list = sorted(sol_list, key= lambda x: x.cost )
    #     sol_list[0].png(f"test/opt_tree_tm{i}.png")
    #     score_results += [score_tree(gt, phi, sol.ct, sol.phi, segments=eval_segs) for sol in sol_list]
    # pd.DataFrame(score_results).to_csv(args.out, index=False)



    