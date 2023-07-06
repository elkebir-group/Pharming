from genotype_tree import GenotypeTree
import networkx as nx 
import numpy as np
from fit_segment_tree import BuildSegmentTree
import pandas as pd 
from data import Data
import argparse
        
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
       
        
class Pharming:
    def __init__(self, max_dcf_clusters=5, nrestarts=10) -> None:
        self.max_clusters = max_dcf_clusters
        self.nrestarts = nrestarts

    def fit_segment(self, s):
        seg_snvs = self.data.seg_to_snvs[s]
        cells_by_cn = self.data.cells_by_cn(s)
        vaf_by_cn = {}
        for cn, cells in cells_by_cn.items():
            print(f"segment: {s} copy number {cn}: #cells {len(cells)}")
            vaf_by_cn[cn] = self.data.compute_vafs(cells = cells, snvs=seg_snvs)
            
            # print(vaf_by_cn[cn])



    def fit(self, data):
        self.data = data 
        self.segments = data.segments
        for s in self.segments:
            self.fit_segment(s)
        print("done")



if __name__ == "__main__":

    # fname = "/scratch/data/leah/pharming/src/test_state_trees.txt"
    # geno_trees = read_genotype_trees(fname)


    # T_CNA = geno_trees[0].generate_CNA_tree()

    # T_SNVs = {0: geno_trees[1], 1: geno_trees[2], 2: geno_trees[3], 3: geno_trees[0]}
    # clusters = np.array([i for i in range(4)], dtype=int)
    # DCF = np.array([[0.855, 0.175, 0, 0.145 ], [1,1,1,0]])
    # cn_profile = np.array([2,2,5,2]).reshape(-1,1)
    # genotypes = np.array([[1,0,0,0], [1,1,0,0],[1,1,1,0], [0,0,0,1]])

    # var =np.full_like(genotypes, 0)
    # var[genotypes==1] = 5
    # total = np.full_like(genotypes, 10)
    # snvs = np.array([i  for i in range(4)], dtype=int)
    # cell_lookup = pd.Series([i for i in range(4)])
    # mut_lookup =pd.Series([i for i in range(4)])
    # snv_to_seg = np.full_like(snvs, 1)

    # data = Data(var, total, cn_profile, snv_to_seg) #cell_lookup, mut_lookup)
    # cells_by_sample= {0: [0,1,3], 1: [2]}
    # cell_assign = BuildSegmentTree(0,T_CNA, T_SNVs, DCF, clusters, cells_by_sample, data ).fit()


    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c" ,"--copy_numbers", required=True,
                        help="input files of copy numbers by segment with unlabeled columns [segment cell totalCN]")
    tpath = "/scratch/data/leah/pharming/test"
    args = parser.parse_args([
        "-f", f"{tpath}/input/read_counts.tsv",
        "-c", f"{tpath}/input/copy_numbers.tsv"
    ])

    print("\nWelcome to the Pharm! Let's start pharming.....\n")
    col_names = ['chr', 'segment', 'mutation_label', 'cell_label','var', 'total']


    
    read_counts = pd.read_table(
        args.file, sep="\t", header=None, names=col_names, skiprows=[0])
    
    read_counts['chr_mutation'] = read_counts['chr'].astype('str') + "_" + \
             read_counts['mutation_label'].astype(str)
    

    cell_labels = np.sort(read_counts['cell_label'].unique())
    mut_labels = np.sort(read_counts['chr_mutation'].unique())



    #create indexed series of mapping of cell index to label
    cell_lookup = pd.Series(data=cell_labels, name="cell_label").rename_axis("cell")     
    mut_lookup = pd.Series(data=mut_labels, name="chr_mutation").rename_axis("mut")



    read_counts = pd.merge(read_counts, cell_lookup.reset_index(), on='cell_label', how='left')
    read_counts = pd.merge(read_counts, mut_lookup.reset_index(), on='chr_mutation', how='left')
   

    copy_numbers = pd.read_table(args.copy_numbers, names=["segment", "cell_label", "cn"])
    segs = copy_numbers.loc[:, ["segment"]].drop_duplicates()
    seg_labels = np.sort(segs['segment'].unique())
    seg_lookup = pd.Series(data=seg_labels, name="segment").rename_axis("seg_id")

    copy_numbers = pd.merge(copy_numbers, cell_lookup.reset_index(), on='cell_label', how='left').drop("cell_label", axis=1)
    copy_numbers= pd.merge(copy_numbers, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)


    
    read_counts = pd.merge(read_counts, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)
    seg_to_mut_mapping = read_counts.loc[:, ["seg_id", "mut"]].drop_duplicates()
    snv_to_seg = seg_to_mut_mapping.set_index("mut")["seg_id"].to_dict()
    seg_to_snvs =  {value: [k for k, v in snv_to_seg.items() if v == value] for value in set(snv_to_seg.values())}


    #check that the copy number data frame contains the exact same cells as read counts
    # cn_cells_in_read_counts = copy_numbers['cell'].isin(read_counts['cell']).all()
    # if cn_cells_in_read_counts.all():
        
    #     print(f"checked that both input datasets contains the same {cell_lookup.shape[0]} cells")
      
    # else:
    #     print("cells appear in copy number input that are not in read counts input these cells will be dropped")
  
    

    

    read_counts= read_counts.set_index(["cell", "mut"])
    var = read_counts["var"].unstack(level="mut", fill_value=0).to_numpy()
    total = read_counts["total"].unstack(level="mut", fill_value=0).to_numpy()
    copy_numbers = copy_numbers.set_index(["seg_id", "cell"])
    copy_numbers= copy_numbers.unstack(level="seg_id", fill_value=0).to_numpy()

    dat = Data(var, total, copy_numbers, snv_to_seg, seg_to_snvs)
    Pharming().fit(dat)



 


 
    
    #create cell indices for 

            #wrangle input data into the correct format
       
        

       
  


        
        # variant_count_data= variant_count_data.set_index(["cell", "mutation"])


 
          
 
        # bin_count_data['cell_index'] = cell_series[bin_count_data['cell']].values
        # bin_count_data = bin_count_data.sort_values(by=['cell_index'])

        # bin_count_data.drop(['cell', 'cell_index'], inplace=True, axis=1)
        # bin_count_data = bin_count_data.to_numpy()
    


    
