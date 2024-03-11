from clonal_tree import ClonalTree
from data import Data 
import pickle 
import argparse 
import pandas as pd 
import numpy as np 
import networkx as nx
from genotype import genotype
import pickle 




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
                # if row[4] in mut_lookup.values:
                    # m = mut_lookup[mut_lookup == row[4]].index[0]
                mut_to_seg[m] = row[1]
    
                g= genotype(row[2], row[3], row[5], row[6])
                genotypes[row[0]][m] = g
    return mut_to_seg
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", type=str,
        help="filename of input tree")
    # parser.add_argument("--phi", type=str, help="ground truth cell assignments")
    parser.add_argument("-g" ,"--genotypes", type=str, help="ground truth genotypes")
    # parser.add_argument("-m", "--mut_lookup",  type=str, 
    #     help="input filename of mapping of internal data index to mutation label ")
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



    # if args.mut_lookup is not None:
    #     mut_lookup = pd.read_csv(args.mut_lookup)
    #     mut_lookup = mut_lookup.set_index("index")
 

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

    mut_to_seg = genotypes_prep(args.genotypes, genotypes)
    seg_to_snvs = {}
    for key, value in mut_to_seg.items():
        if value not in seg_to_snvs:
            seg_to_snvs[value] = [key]
        else:
            seg_to_snvs[value].append(key)
    
    # cell_mapping = {v: [] for v in tree}
    ct = ClonalTree(tree, genotypes, seg_to_muts= seg_to_snvs )
    # if args.phi is not None:
    #     cell_assignment = pd.read_csv(args.phi)
    #     phi = {}
    #     series = dat.cell_lookup
    #     # print(cell_assignment.head())
    #     for index, row in cell_assignment.iterrows():

    #         i = row['Cell']
    #         v = row['Cluster']

         
    #         filtered_series = series[series == i]

    #         # Get the corresponding indices
    #         index = filtered_series.index.tolist()[0]
    #         phi[index] = v

    #         # cell_mapping[v].append(i)

    #     ca = CellAssign(phi, ct.clones())
    #     if args.cell_assign is not None:
    #         ca.save(args.cell_assign)


    #     ct.trim(ca)

    #     costs = ct.compute_costs(dat, ca, lamb=0.25)
    #     if args.draw is not None:
    #         ct.draw(args.draw, ca)
    # else:
    if args.draw is not None:
            ct.draw(args.draw)
    # missing_muts = set(dat.muts) - set(ct.get_all_muts())
    if args.clonal_tree is not None:
        ct.save(args.clonal_tree)
    # if args.phi is not None:
        # cost2 = ct.compute_pooled_costs(data, lamb=100)
 
        # print(f"cost 1: {costs} cost 2: {cost2}")
        













