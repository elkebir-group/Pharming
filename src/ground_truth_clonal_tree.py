# Created by: L.L. Weber
# Created on: 2024-03-14 13:01:31


from clonal_tree import ClonalTree
from data import Data
from cell_mapping import CellAssign
from utils import load_pickled_object
import argparse 
import pandas as pd 
import networkx as nx
from genotype import genotype



def genotypes_prep(genotypes_fname, genotypes):
    firstline = True
    mut_to_seg = {}

    with open(genotypes_fname, "r+") as file:
        for line in file:
            if firstline:
                firstline= False
                continue
            else:
                row = line.strip().split("\t")
                row = [int(r) for r in row]
                m = row[4]
                mut_to_seg[m] = row[1]
    
                g= genotype(row[2], row[3], row[5], row[6])
                genotypes[row[0]][m] = g
    return mut_to_seg
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", type=str, help="filename of input tree")
    parser.add_argument("-g" ,"--genotypes", type=str, help="ground truth genotypes")
    parser.add_argument("-d" ,"--data", type=str, 
                        help="filename of data object, if given, the gt genotypes will be filtered/reindexed")
    parser.add_argument("-p", "--phi", required=True,
                        help="input file for cell assignment")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=1e3)
    parser.add_argument("-D" ,"--dcfs", required=False, type=str)
    parser.add_argument("-P", "--assign",
                        help="input file for cell assignment")
    parser.add_argument( "--mut-cluster-tree", required=False, type=str)
    parser.add_argument("-T", "--clonal-tree",  type=str, 
        help="output filename of pickled clonal tree object")
    parser.add_argument( "--draw",  type=str, 
        help="png or pdf of output tree")


    args = parser.parse_args()


    # instance = "s11_m5000_k25_l7"
    # # instance = "s12_m5000_k25_l7"
    # folder = "n1000_c0.05_e0" 
    # pth = f"simulation_study/input"
    # opth = "test"

    # args = parser.parse_args([

    #     # "-d", f"{pth}/{instance}/{folder}/data.pkl",
    #     "-d", f"{opth}/data.pkl",
    #     "-g", f"{pth}/{instance}/node.tsv",
    #     "-t", f"{pth}/{instance}/tree.tsv",
    #     "-p", f"{pth}/{instance}/{folder}/cellAssignments.p0",
    #     "-l", "1e3",
    #     "-D", f"{opth}/dcfs.txt",
    #     "-T", f"{opth}/gt.pkl",
    #     "-P", f"{opth}/phi.pkl",
    #     "--mut-cluster-tree", f"{opth}/Tm.txt",
    #     "--draw", f"{opth}/gt.png"
    # ])


    
    #construct the tree from the edges 
    tree = nx.DiGraph()
    with open(args.tree, "r+") as file:
        for line in file:
               vals = line.strip().split("\t")
               vals = [int(v) for v in vals]
               tree.add_edge(vals[0], vals[1])

    #create the genotypes for the tree
    genotypes = {v:  {} for v in tree}

    mut_to_seg = genotypes_prep(args.genotypes, genotypes)
    seg_to_snvs = {}
    for key, value in mut_to_seg.items():
        if value not in seg_to_snvs:
            seg_to_snvs[value] = [key]
        else:
            seg_to_snvs[value].append(key)
    
    #make ground truth clonal tree
    gt = ClonalTree(tree, genotypes, seg_to_muts= seg_to_snvs )

    #reindex the SNVs to make the indices in the data object
    #filter to include only SNVs in data object
    if args.data is not None:
        dat = load_pickled_object(args.data)
        gt.reindex_snvs(dat.mut_lookup)
        gt.filter_snvs(dat.muts)


    #make ground truth cell assignment phi
    cell_assignment = pd.read_csv(args.phi)
    phi = {}

    for index, row in cell_assignment.iterrows():

        i = row['Cell']
        v = row['Cluster']
        phi[index] = v


    ca = CellAssign(phi, gt.clones())

    ca.relabel(dat.cell_lookup)
    if args.assign is not None:
        ca.save(args.assign)
   
 
    #get additional useful files like DCFs and Tm
    gt_dcfs = gt.compute_dcfs(ca)
    root = gt.root


    gt_T_m = gt.mutation_cluster_tree()
    gt_T_m.remove_node(root)




    mapping = {}
    nodes = list(gt_T_m.nodes)
    nodes.sort()
    for i,n in enumerate(nodes):
          mapping[n] =i 
    
    T_m = nx.relabel_nodes(gt_T_m, mapping)
    if args.mut_cluster_tree is not None:
        with open(args.mut_cluster_tree, "w+") as file:
            for u,v in T_m.edges:
                file.write(f"{u}\t{v}\n")


    gt_delta = {mapping[n]: gt_dcfs[n] for n in mapping if n != root}

    if args.dcfs is not None:
        with open(args.dcfs, "w+") as file:
            for key,val in gt_delta.items():

                file.write(f"{val}\n")

    gt.compute_likelihood(dat, ca, args.lamb )

    if args.draw is not None:
        gt.draw(args.draw, ca, segments=dat.segments, include_dcfs=True)

    if args.clonal_tree is not None:
        gt.save(args.clonal_tree)

        










