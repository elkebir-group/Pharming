import re
import argparse
import pandas as pd 
import pygraphviz as pgv
import networkx as nx
from clonal_tree import ClonalTree
from cell_mapping import CellAssign
from solution import Solution




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



def get_consenus_profile(phi, n, tree, dat, root):
    copy_number_node = None
    if n not in phi:
        reverse_root_to_node_path =nx.shortest_path(tree, root, n)[::-1]
        for u in reverse_root_to_node_path:
            if u in phi:
                copy_number_node = u
                break
        
    else: 
        copy_number_node = n
    if copy_number_node is None:
        return {ell: (1,1) for ell in dat.segments}
    
    cells = phi[copy_number_node]

    return {ell: dat.consensus_profile(cells, ell) for ell in dat.segments}

        


def make_genotypes(tree, mut_mapping, phi, dat):
    root = [n for n in tree if tree.in_degree[n]==0][0]
    snvs = set(dat.snv_to_seg.keys())


    genotypes = {u: {} for u in tree}
    for v in tree:
        if v not in mut_mapping:
            mut_mapping[v] =[]
    for u in tree:
        consensus_profiles = get_consenus_profile(phi, u, tree, dat, root)
        present = nx.ancestors(tree, u) | {u}
        pres_snvs =set(int(j) for v in present for j in mut_mapping[v])
        not_present = snvs - pres_snvs
        for j in pres_snvs:
            ell = dat.snv_to_seg[j]
            cn = consensus_profiles[ell]
            genotypes[u][j] =(*cn,1,0)
        for j in not_present:
            ell = dat.snv_to_seg[j]
            cn = consensus_profiles[ell]
            genotypes[u][j] =(*cn,0,0)
    return genotypes


def collapse_linear_chains(T, root, mut_mapping, phi):
    to_remove  =[]
    for u in nx.dfs_postorder_nodes(T, source = root):
        if T.out_degree[u] ==1 and u not in phi and u != root:
            to_remove.append(u)
            for v in nx.dfs_preorder_nodes(T, source = u):
                if v in phi or T.out_degree[v] >1:
                    mut_mapping[v].extend(mut_mapping[u]) 
                    del mut_mapping[u]
                    break
    for n in to_remove:
        parent =list(T.predecessors(n))[0]
        children = list(T.successors(n))
        T.remove_node(n)
        for c in children:
            T.add_edge(parent, c)


if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--scite-tree", required=True,
                        help="input file of ground truth tree")
    parser.add_argument("-g", "--scite-genes", required=True,
                        help="input file of ground truth tree")
    parser.add_argument("-m", "--mut-clusters", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-clusters", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=1000,
                        help="lambda value, default=1000")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="path where the pickled solution wil be saved")
    parser.add_argument("--png", required=False, type=str,
                        help="path where the pickled solution wil be saved")


    
    args = parser.parse_args()
    # instance = "sims4/s14_m5000_k25_l5_d2/n1000_c0.05_e0"
    # args = parser.parse_args(["-d", f"simulation_study/{instance}/data.pkl",
    #                            "-t", f"simulation_study/baseline/{instance}/scite_ml0.gv", 
    #                            "-g", f"simulation_study/baseline/{instance}/scite.genes", 
    #                            "-m", f"simulation_study/baseline/{instance}/mut_cluster.csv", 
    #                            "-c", f"simulation_study/baseline/{instance}/cell_clusters.csv", 
    #                            "-l", "1000", "-o", f"simulation_study/baseline/{instance}/solution.pkl", 
    #                            "--png", f"simulation_study/baseline/{instance}/clonal_tree0.png"])  

    dat = pd.read_pickle(args.data)


    G=pgv.AGraph(args.scite_tree)
    

    mut_clust_assign, mut_cluster_mapping = create_mappings(args.mut_clusters)
    cell_clust_assign, cell_cluster_mapping = create_mappings(args.cell_clusters)


    genes = []
    name_to_cluster ={}
    with open(args.scite_genes, "r+") as f:
        for i,line in enumerate(f):
            genes.append(line.strip())
            name_to_cluster[line.strip()] = i

    genes.append("Root")
    name_to_cluster["Root"] = len(genes)-1
            




    tree = nx.nx_agraph.from_agraph(G)


    T = nx.DiGraph()
    mut_mapping = {}
    phi = {}
    for u,v, _ in tree.edges:
        print(f"{u} -> {v}")
        if not re.search(r"s[0-9]+", v):
            T.add_edge(name_to_cluster[u], name_to_cluster[v])
            mut_mapping[name_to_cluster[v]] = mut_cluster_mapping[name_to_cluster[v]]

        else:
            cell_cluster = int(v.replace("s", ""))
            if name_to_cluster[u] in phi:
                phi[name_to_cluster[u]].extend(cell_cluster_mapping[cell_cluster])
            else:
                phi[name_to_cluster[u]] = cell_cluster_mapping[cell_cluster]

    collapse_linear_chains(T, name_to_cluster["Root"], mut_mapping, phi)


    genotypes = make_genotypes(T, mut_mapping, phi, dat)
    ct = ClonalTree(T, genotypes, dat.seg_to_snvs)
    phi_reverse = {v: k for k, vals in phi.items() for v in vals}
    cell_assign = CellAssign(phi_reverse, set(T.nodes))

    likelihood = ct.compute_likelihood(dat, cell_assign, lamb=args.lamb)
    sol = Solution(likelihood, ct, cell_assign)
    if args.out is not None:
        pd.to_pickle(sol, args.out)
    if args.png is not None:
        sol.png(args.png, segments=dat.segments)
      
