import argparse
from clonal_tree import ClonalTree, load
from cell_mapping import CellAssign
import pandas as pd 
from sklearn.metrics.cluster import adjusted_rand_score
import utils
import itertools
import numpy as np 

from data import Data, load_from_pickle



def score_tree(gt, gt_phi, inf,inf_phi, segments):
      
        scores = gt.score_snvs(inf)
        # scores["tree"] = i
        scores["segment"] = ":".join([str(ell) for ell in segments])
        scores["cell_ari"] =gt_phi.compute_ari(inf_phi)
        scores["gt_cost"] = gt.cost 
        scores["inf_cost"] = inf.cost 
        scores['cna_mad'] = gt.cna_genotype_similarity(gt_phi, inf, inf_phi)
        
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
    dat = load_from_pickle(args.data)


    if args.segment is None:
        eval_segs = [ell for ell in dat.segments if dat.num_cn_states(ell) > 1]
    else:
        eval_segs = args.segment 

    snvs = list(itertools.chain(*[dat.seg_to_snvs[seg] for seg in eval_segs]))

    gt = load_from_pickle(args.tree)
    phi = load(args.cell_assign)

    gt.reindex_snvs(dat.mut_lookup)

    phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    root = gt.root

    gt.filter_snvs(snvs)

    cost = gt.compute_likelihood(dat, phi, lamb)


    sol_list  = load_from_pickle(args.solutions)

    # score_results= [score_tree(gt, phi, sol.ct, sol.phi, segments=eval_segs) for sol in sol_list]
    # pd.DataFrame(score_results).to_csv(args.out, index=False)
    sol_list[0].png("test/opt_tree.png")


    score_results = []
    clonal_trees  = load_from_pickle("test/clonal_trees.pkl")
    for i,sol_list in enumerate(clonal_trees):
        sol_list[0].png(f"test/opt_tree_tm{i}.png")
        score_results += [score_tree(gt, phi, sol.ct, sol.phi, segments=eval_segs) for sol in sol_list]
    pd.DataFrame(score_results).to_csv(args.out, index=False)



    