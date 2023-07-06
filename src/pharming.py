from genotype_tree import GenotypeTree
import networkx as nx 
import numpy as np
from fit_segment_tree import BuildSegmentTree
import pandas as pd 
from data import Data

        
def convert_tree_string(edge_string,id):
    tree = nx.DiGraph()
    genotypes = {}
    index = 0
    node_id = -1
    node_mapping = {}
  
    while index < len(edge_string):
        u = edge_string[index].split(",")
        v = edge_string[index+1].split(",")
        geno_u = tuple([int(i) for i in u])
        geno_v = tuple([int(i) for i in v])
        if geno_u not in node_mapping:
            node_id += 1
            node_mapping[geno_u] = node_id
            genotypes[node_id] = geno_u

        if geno_v not in node_mapping:
            node_id +=1
            node_mapping[geno_v] = node_id
            genotypes[node_id] = geno_v
        u_node_id = node_mapping[geno_u]
        v_node_id = node_mapping[geno_v]
        tree.add_edge(u_node_id, v_node_id)
        index += 2
    
    return GenotypeTree(tree, genotypes, node_mapping, id=id)
    




    

        
        
            
def read_genotype_trees(fname):
    genotype_trees = []
    id =0
    with open(fname, "r+") as file:
        state_tree_section = False
        for line in file:
            

            line = line.strip()
            if "state_trees" in line:
                state_tree_section = True
                continue

            if state_tree_section:
                tree_string = line.split()
                
                genotype_trees.append(convert_tree_string(tree_string, id))
                id += 1
                
    return genotype_trees 
       
        


if __name__ == "__main__":

    fname = "/scratch/data/leah/phertilizer2.0/src/test_state_trees.txt"
    geno_trees = read_genotype_trees(fname)


    T_CNA = geno_trees[0].generate_CNA_tree()

    T_SNVs = {0: geno_trees[1], 1: geno_trees[2], 2: geno_trees[3], 3: geno_trees[0]}
    clusters = np.array([i for i in range(4)], dtype=int)
    DCF = np.array([[0.855, 0.175, 0, 0.145 ], [1,1,1,0]])
    cn_profile = np.array([2,2,5,2]).reshape(-1,1)
    genotypes = np.array([[1,0,0,0], [1,1,0,0],[1,1,1,0], [0,0,0,1]])

    var =np.full_like(genotypes, 0)
    var[genotypes==1] = 5
    total = np.full_like(genotypes, 10)
    snvs = np.array([i  for i in range(4)], dtype=int)
    cell_lookup = pd.Series([i for i in range(4)])
    mut_lookup =pd.Series([i for i in range(4)])
    snv_to_seg = np.full_like(snvs, 1)

    data = Data(var, total, cn_profile, snv_to_seg, cell_lookup, mut_lookup)
    cells_by_sample= {0: [0,1,3], 1: [2]}
    cell_assign = BuildSegmentTree(0,T_CNA, T_SNVs, DCF, clusters, cells_by_sample, data ).fit()

    


    
