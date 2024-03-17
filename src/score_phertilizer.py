
from utils import load_pickled_object
from clonal_tree import ClonalTree
import networkx as nx 
from genotype import genotype
from cell_mapping import CellAssign
import score_tree
import pandas as pd 
import argparse



def make_clonal_tree(phert):

    tree = phert.tree
    genotypes = {u: {} for u in tree}
    snvs = set(phert.get_all_muts())
    for u in phert.tree:
        present = nx.ancestors(tree, u) | {u}
        pres_snvs =set(j for v in present for j in phert.mut_mapping[v] )
        not_present = snvs - pres_snvs
        for j in pres_snvs:
            genotypes[u][j] = genotype(1,1,1,0)
        for j in not_present:
            genotypes[u][j] = genotype(1,1,0,0)
    root= [n for n in tree if tree.in_degree[n]==0][0]
    new_root = max(tree)+1
    tree.add_edge(new_root, root)
    genotypes[new_root] = {}
    for j in snvs:

        genotypes[new_root][j]= genotype(1,1,0,0)
    return ClonalTree(tree, genotypes, {0: snvs})
 


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


    
    args = parser.parse_args()


    # instance = "s11_m5000_k25_l7"
    # folder = "n1000_c0.05_e0" 
    # ppth = f"simulation_study/input"


    # args = parser.parse_args([

    #     "-d", f"simulation_study/input/{instance}/{folder}/data.pkl",
    #     "-t", f"simulation_study/input/{instance}/{folder}/gt.pkl",
    #     "-c", f"simulation_study/input/{instance}/{folder}/phi.pkl",
    #     "-T", f"simulation_study/phertilizer/{instance}/{folder}/tree.pkl",
    #     "-o", f"test/scores.csv"


    # ])
    dat = load_pickled_object(args.data)


    phert = load_pickled_object(args.phert_tree)
    gt = load_pickled_object(args.gt_tree)
    phi = load_pickled_object(args.cell_assign)

    ct  = make_clonal_tree(phert)
    ca = make_phi(phert)
    ct.compute_likelihood(dat,ca, args.lamb )


  
    # ct.draw("test/phert.png", ca)
    # gt.draw("test/gt.png", phi)
    score = score_tree.score_tree(gt, phi, ct, ca, segments=[0])
    index = ['Row1']
    pd.DataFrame(score, index=index).to_csv(args.out, index=False)

