from clonal_tree import ClonalTree
from data import Data, load_from_files
import argparse 
import pandas as pd 
import numpy as np 
import networkx as nx
from genotype import genotype



def genotypes_prep(genotypes_fname, genotypes, mut_lookup):
    firstline = True

    with open(genotypes_fname, "r+") as file:
        for line in file:
            if firstline:
                firstline= False
                continue
            else:
                row = line.strip().split("\t")
                row = [int(r) for r in row]
                if row[4] in mut_lookup.values:
                    m = mut_lookup[mut_lookup == row[4]].index[0]
       
                    g= genotype(row[2], row[3], row[5], row[6])
                    genotypes[row[0]][m] = g
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-t", "--tree", type=str,
        help="filename of input tree")
    parser.add_argument("--phi", type=str, help="ground truth cell assignments")
    parser.add_argument("-g" ,"--genotypes", type=str, help="ground truth genotypes")

    parser.add_argument("-D", "--data",  type=str, 
        help="output filename of pickled data object")
    parser.add_argument("-T", "--clonal-tree",  type=str, 
        help="output filename of pickled clonal tree object")
    parser.add_argument( "--draw",  type=str, 
        help="png or pdf of output tree")


    # args = parser.parse_args()
    instance = "s12_n5000_m5000_k25_c0.1_l7"
    pth = f"/scratch/data/leah/pharming/simulation_study/input/{instance}"
    args = parser.parse_args(
        ["-f", f"{pth}/sparse.p0",
         "-c" ,f"{pth}/cells.p0",
         "-g", f"{pth}/genotypes.txt",
         "--phi", f"{pth}/cellAssignments.p0",
         "-t", f"{pth}/tree.txt",
         "-T", f"test/ground_truth.pickle",
         "-D", f"{pth}/data.pickle",
         "--draw", f"test/tree.png"
        ]
    )

    data = load_from_files(args.file, args.profiles)
    print(data)
    if args.data is not None:
        data.save(args.data)
    

    # for s in data.segments:
    #     print(f"{s}: {data.cn_states_by_seg(s)}")
    #     cell_state_map = data.cells_by_cn(s)
    
    tree = nx.DiGraph()
    with open(args.tree, "r+") as file:
        firstline =True
   
        for line in file:
            if firstline:
                firstline = False
                continue
                
            else:
               vals = line.strip().split("\t")
               vals = [int(v) for v in vals]
               tree.add_edge(vals[0], vals[1])
    
    genotypes = {v:  {} for v in tree}

    genotypes_prep(args.genotypes, genotypes, data.mut_lookup)
    
    cell_mapping = {v: [] for v in tree}
    
    if args.phi is not None:
        cell_assignment = pd.read_csv(args.phi)
        # print(cell_assignment.head())
        for index, row in cell_assignment.iterrows():
            i = row['Cell']
            v = row['Cluster']
            cell_mapping[v].append(i)


    ct = ClonalTree(tree, genotypes, cell_mapping=cell_mapping )
    missing_muts = set(data.muts) - set(ct.get_all_muts())

    if args.phi is not None:
        cost2 = ct.compute_pooled_costs(data, lamb=100)
        costs = ct.compute_costs(data, lamb=100)
        print(f"cost 1: {costs} cost 2: {cost2}")
        



    if args.draw is not None:
        ct.draw(args.draw)
    if args.clonal_tree is not None:
        ct.save(args.clonal_tree)








