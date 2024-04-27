
from utils import load_pickled_object, pickle_object
from clonal_tree import ClonalTree
import networkx as nx 
from solution import Solution
from cell_mapping import CellAssign
import argparse
import pandas as pd 
import utils 

def parse_tree(fname):
    nedges =0
    edge_list = []
    with open(fname, "r+") as file:
        for i, line in enumerate(file):
            line = line.strip().split(" ")
            if i==0:
                nedges = int(line[0])
            elif i <= nedges:
                edge_list.append((int(line[0]), int(line[1])))
            else:
              break
    return nx.DiGraph(edge_list)



def make_genotypes(tree, psi):

    mut_mapping = utils.inverse_dict(psi)
    snvs = set(psi.keys())
    genotypes = {u: {} for u in tree}
    for v in tree:
        if v not in mut_mapping:
            mut_mapping[v] =[]
    for u in tree:
        present = nx.ancestors(tree, u) | {u}
        pres_snvs =set(int(j) for v in present for j in mut_mapping[v])
        not_present = snvs - pres_snvs
        for j in pres_snvs:
            genotypes[u][j] =(1,1,1,0)
        for j in not_present:
            genotypes[u][j] =(1,1,0,0)
    return genotypes
    # root= [n for n in tree if tree.in_degree[n]==0][0]
    # new_root = max(tree)+1
    # tree.add_edge(new_root, root)
    # genotypes[new_root] = {}
    # for j in snvs:

    #     genotypes[new_root][j]= (1,1,0,0)
 
 

# def convert_to_clonal_tree(phert):

#     tree = phert.tree
#     genotypes = {u: {} for u in tree}
#     snvs = set(phert.get_all_muts())
#     snvs = set(int(j) for j in snvs)
#     for u in phert.tree:
#         present = nx.ancestors(tree, u) | {u}
#         pres_snvs =set(int(j) for v in present for j in phert.mut_mapping[v] )
#         not_present = snvs - pres_snvs
#         for j in pres_snvs:
#             genotypes[u][j] =(1,1,1,0)
#         for j in not_present:
#             genotypes[u][j] =(1,1,0,0)
#     # root= [n for n in tree if tree.in_degree[n]==0][0]
#     # new_root = max(tree)+1
#     # tree.add_edge(new_root, root)
#     # genotypes[new_root] = {}
#     # for j in snvs:

#     #     genotypes[new_root][j]= (1,1,0,0)
#     return tree, genotypes
 


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
    parser.add_argument("-t", "--phert-tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-n", "--phi", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-m", "--psi", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=0,
                        help="lambda value, default=1")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory the pickled solution wil be saved")


    
    args = parser.parse_args()

    # seed = 10
    # cov = 0.25
    # instance = f"s{seed}_m5000_k25_l5_d2"
    # folder = f"n1000_c{0.25}_e0" 
    # ppth = f"simulation_study/phertilizer/sims"


    # args = parser.parse_args([

    #     "-d", f"simulation_study/sims/{instance}/{folder}/data.pkl",
    #     "-t", f"{ppth}/{instance}/{folder}/tree.txt",
    #     "-n", f"{ppth}/{instance}/{folder}/phi.csv",
    #     "-m", f"{ppth}/{instance}/{folder}/psi.csv",
    #     "-o", f"test/sol.pkl"


    # ])
    dat = load_pickled_object(args.data)
    T= parse_tree(args.phert_tree)
    cell_df =  pd.read_csv(args.phi)
    phi_dict =dict(zip(cell_df["cell"], cell_df["cluster"]))
    phi = CellAssign(phi_dict, set(T.nodes()))
    df =  pd.read_csv(args.psi)
    df[['segment', 'snv']] = df['mutation'].str.split('_', expand=True)
    df["segment"] = df["segment"].astype(int)
    df["snv"] = df["snv"].astype(int)
    psi = dict(zip(df["snv"], df["cluster"]))
    snv_to_seg = dict(zip(df["snv"], df["segment"]))
    seg_to_snvs = utils.inverse_dict(snv_to_seg)

    genotypes = make_genotypes(T, psi)

    ct = ClonalTree(T, genotypes, seg_to_snvs)
  
    cost =ct.compute_likelihood(dat,phi, args.lamb )
    sol = Solution(cost, ct, phi)
    if args.out is not None:
        pickle_object(sol, args.out)
    
    
    


   

    # score = score_tree.score_tree(gt, phi, ct, ca, dat, args.lamb, segments=[0])
    # index = ['Row1']
    # pd.DataFrame(score, index=index).to_csv(args.out, index=False)

