from clonal_tree import ClonalTree
from data import Data 
import pickle 
import argparse 
import pandas as pd 
import numpy as np 
import networkx as nx
from genotype import genotype
import pickle 



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
    parser.add_argument("-t", "--tree", type=str,
        help="filename of input tree")
    parser.add_argument("--phi", type=str, help="ground truth cell assignments")
    parser.add_argument("-g" ,"--genotypes", type=str, help="ground truth genotypes")
    parser.add_argument("-d", "--data",  type=str, 
        help="input filename of pickled data object")
    parser.add_argument("-T", "--clonal-tree",  type=str, 
        help="output filename of pickled clonal tree object")
    parser.add_argument( "--draw",  type=str, 
        help="png or pdf of output tree")


    args = parser.parse_args()
    # instance = "s12_m5000_k25_l7"
    # dpath = "n500_c0.1"
    # pth = f"simulation_study/input/{instance}"
    # args = parser.parse_args([

    #      "-g", f"{pth}/node.tsv",
    #      "--phi", f"{pth}/{dpath}/cellAssignments.p0",
    #      "-t", f"{pth}/tree.tsv",
    #      "-T", f"test/ground_truth.pickle",
    #      "-d", f"{pth}/{dpath}/data.pickle",
    #      "--draw", f"test/tree.png"
    #     ]
    # )



    if args.data is not None:
        with open(args.data, "rb") as file:
            dat = pickle.load(file)

    # for key,val in data.cells_by_cn(1).items():
    #     print(f"State:{key} # of cells: {len(val)}")
    

    # for s in data.segments:
    #     print(f"{s}: {data.cn_states_by_seg(s)}")
    #     cell_state_map = data.cells_by_cn(s)
    
    tree = nx.DiGraph()
    with open(args.tree, "r+") as file:
        for line in file:
               vals = line.strip().split("\t")
               vals = [int(v) for v in vals]
               tree.add_edge(vals[0], vals[1])

    genotypes = {v:  {} for v in tree}

    genotypes_prep(args.genotypes, genotypes, dat.mut_lookup)
    
    cell_mapping = {v: [] for v in tree}
    
    if args.phi is not None:
        cell_assignment = pd.read_csv(args.phi)
        series = dat.cell_lookup
        # print(cell_assignment.head())
        for index, row in cell_assignment.iterrows():

            i = row['Cell']
            v = row['Cluster']
         
            filtered_series = series[series == i]

            # Get the corresponding indices
            indices = filtered_series.index.tolist()[0]

            cell_mapping[v].append(i)


    ct = ClonalTree(tree, genotypes, seg_to_muts=dat.seg_to_snvs, cell_mapping=cell_mapping )
    ct.trim()
    missing_muts = set(dat.muts) - set(ct.get_all_muts())

    if args.phi is not None:
        # cost2 = ct.compute_pooled_costs(data, lamb=100)
        costs = ct.compute_costs(dat, lamb=100)
        # print(f"cost 1: {costs} cost 2: {cost2}")
        



    if args.draw is not None:
        ct.draw(args.draw)
    if args.clonal_tree is not None:
        ct.save(args.clonal_tree)








