import argparse
import pandas as pd 
import numpy as np
import networkx as nx
from clonal_tree import ClonalTree
from small_parsimony import SmallParsimony
from Bio import Phylo



def create_mappings(fname):
    clusters = pd.read_csv(fname)
    clust_assign = dict(zip(clusters["index"], clusters["cluster"]))
    cluster_mapping = {}
    for key,val in clust_assign.items():
        if val in cluster_mapping:
            cluster_mapping[val].append( key)
        else:
            cluster_mapping[val] = [key]
    return clust_assign, cluster_mapping


        


def make_genotypes(tree, dat, node_mapping, x_profiles, y_profiles):
    snvs = set(dat.snv_to_seg.keys())
    root = [n for n in tree if tree.in_degree[n]==0][0]

    genotypes = {u: {} for u in tree}

      
    for u in tree:
        u_str = node_mapping[u]
        x_prof = x_profiles[u_str]
        y_prof = y_profiles[u_str]

        for j in snvs:
            ell = dat.snv_to_seg[j]
            cn = (x_prof[ell], y_prof[ell])
            if u != root:
                genotypes[u][j] =(*cn,1,0)
            else:
                genotypes[u][j] =(*cn,0,0)
   
    return genotypes


def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    # new_net_tree = net_tree.copy()
    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'
    
    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree




if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--newick", required=True,
                        help="input file of ground truth tree")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="path where the pickled solution wil be saved")



    
    args = parser.parse_args()
    # instance = "sims4/s14_m5000_k25_l5_d2/n1000_c0.05_e0"
    # args = parser.parse_args(["-d", f"simulation_study/{instance}/data.pkl",
    #                            "-t", f"simulation_study/lazac/{instance}/output_tree.rooted.newick", 
    #                            "-l", "1000", "-o", f"simulation_study/lazac/{instance}/solution.pkl", 
    #                            "--png", f"simulation_study/lazac/{instance}/clonal_tree0.png"])  

    dat = pd.read_pickle(args.data)
    tree = from_newick_get_nx_tree(args.newick)
    leaves =[n for n in tree if tree.out_degree[n]==0]
    max_cn =np.max([dat.copy_y.max(), dat.copy_x.max()])
    min_cn = np.min([dat.copy_y.min(), dat.copy_x.min()])
    root = [n for n in tree if tree.in_degree[n]==0][0]
    sp =SmallParsimony(tree, root, alphabet=list(range(min_cn, max_cn+1)))

    x_profiles = {l: list(dat.copy_x[int(l),]) for l in leaves}
    y_profiles = {l: dat.copy_y[int(l),] for l in leaves}

    
    x_profiles[root] = [1 for _ in range(len(dat.segments))]
    y_profiles[root] = [1 for _ in range(len(dat.segments))]
    score_x, x_labels =sp.sankoff(x_profiles)
    score_y, y_labels = sp.sankoff(y_profiles)

    node_mapping = {n: i for i,n in enumerate(tree.nodes)}
    rev_node_mapping = {i: n for n,i in node_mapping.items()}
    T = nx.relabel_nodes(tree, node_mapping, copy=tree)


    ct = ClonalTree(T, make_genotypes(T, dat, rev_node_mapping, x_labels, y_labels), dat.seg_to_snvs)
    if args.out is not None:
        pd.to_pickle(ct, args.out)

    # ct.draw("test/lazac.png", segments=dat.segments)









   
   




    # tree = nx.nx_agraph.from_agraph(G)


    # T = nx.DiGraph()
    # mut_mapping = {}
    # phi = {}
    # for u,v, _ in tree.edges:
    #     print(f"{u} -> {v}")
    #     if not re.search(r"s[0-9]+", v):
    #         T.add_edge(name_to_cluster[u], name_to_cluster[v])
    #         mut_mapping[name_to_cluster[v]] = mut_cluster_mapping[name_to_cluster[v]]

    #     else:
    #         cell_cluster = int(v.replace("s", ""))
    #         if name_to_cluster[u] in phi:
    #             phi[name_to_cluster[u]].extend(cell_cluster_mapping[cell_cluster])
    #         else:
    #             phi[name_to_cluster[u]] = cell_cluster_mapping[cell_cluster]


    # if args.collapse:
    #     collapse_linear_chains(T, name_to_cluster["Root"], mut_mapping, phi)


    # genotypes = make_genotypes(T, mut_mapping, phi, dat)
    # ct = ClonalTree(T, genotypes, dat.seg_to_snvs)
    # phi_reverse = {v: k for k, vals in phi.items() for v in vals}
    # cell_assign = CellAssign(phi_reverse, set(T.nodes))

    # likelihood = ct.compute_likelihood(dat, cell_assign, lamb=args.lamb)
    # sol = Solution(likelihood, ct, cell_assign)


      
