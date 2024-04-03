
from utils import load_pickled_object
from clonal_tree import ClonalTree
import networkx as nx 

from cell_mapping import CellAssign
import score_tree
import pandas as pd 
import argparse



def make_clonal_tree(phert):

    tree = phert.tree
    genotypes = {u: {} for u in tree}
    snvs = set(phert.get_all_muts())
    snvs = set(int(j) for j in snvs)
    for u in phert.tree:
        present = nx.ancestors(tree, u) | {u}
        pres_snvs =set(int(j) for v in present for j in phert.mut_mapping[v] )
        not_present = snvs - pres_snvs
        for j in pres_snvs:
            genotypes[u][j] =(1,1,1,0)
        for j in not_present:
            genotypes[u][j] =(1,1,0,0)
    root= [n for n in tree if tree.in_degree[n]==0][0]
    new_root = max(tree)+1
    tree.add_edge(new_root, root)
    genotypes[new_root] = {}
    for j in snvs:

        genotypes[new_root][j]= (1,1,0,0)
    return ClonalTree(tree, genotypes, {0: list(snvs)})
 


def make_phi(phert):
    phi = {}
    cell_map = phert.cell_mapping
    clones = set(phert.tree.nodes)
    for u, cells in cell_map.items():
        for i in cells[0]:
            phi[i] =u
    return CellAssign(phi, clones)





if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--gt-tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-T", "--phert-tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=0,
                        help="lambda value, default=1")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory where score files should be written")

    # parser.add_argument("--gt_phi", required=False, type=str,
    #                     help="directory where ground truth phi saved")
    # parser.add_argument("--gt_tree", required=False, type=str,
    #                     help="directory where ground truth psi saved")

    
    args = parser.parse_args()

    # seed = 10
    # cov = 0.25
    # instance = f"s{seed}_m5000_k25_l5"
    # folder = f"n1000_c{0.25}_e0" 
    # ppth = f"simulation_study/test"


    # args = parser.parse_args([

    #     "-d", f"{ppth}/{instance}/{folder}/data.pkl",
    #     "-t", f"{ppth}/{instance}/{folder}/gt.pkl",
    #     "-c", f"{ppth}/{instance}/{folder}/phi.pkl",
    #     "-T", f"simulation_study/phertilizer/test/{instance}/{folder}/tree.pkl",
    #     "-o", f"test/scores.csv"


    # ])
    dat = load_pickled_object(args.data)


    phert = load_pickled_object(args.phert_tree)
    snvs = phert.get_all_muts()
    snvs = set(int(i) for i in snvs)
    # for u, snvs in phert.mut_mapping.items():
    #     print(f"{u}: {len(snvs)}")

    gt = load_pickled_object(args.gt_tree)
    gt_snvs = set(gt.get_all_muts())

    diff_snvs = snvs.symmetric_difference(gt_snvs)
    print(f"Number of SNVs in symmetric difference: {diff_snvs}")
 
    phi = load_pickled_object(args.cell_assign)

    # if args.gt_psi is not None:
    #     gt.write_psi(args.gt_psi)
    
    # if args.gt_phi is not None:
    #     phi.write_phi(args.gt_phi)

    ct  = make_clonal_tree(phert)
    # print(len(ct.get_all_muts()))

    ca = make_phi(phert)

      
    # ct.draw("test/phert.png", ca)
    # gt.draw("test/gt.png", phi)

    ct.compute_likelihood(dat,ca, args.lamb )



    score = score_tree.score_tree(gt, phi, ct, ca, dat, args.lamb, segments=[0])
    index = ['Row1']
    pd.DataFrame(score, index=index).to_csv(args.out, index=False)

