from genotype_tree import GenotypeTree
import networkx as nx 
import numpy as np
from fit_segment_tree import BuildSegmentTree
from dcf_clustering import DCF_Clustering
import pandas as pd 
from data import Data
import argparse
from clonal_tree import ClonalTree
        

       
        
class Pharming:
    def __init__(self, T_CNAs, T_SNVs, max_dcf_clusters=5, nrestarts=25, rng=None) -> None:
        self.T_CNAs = T_CNAs 
        self.T_SNVs = T_SNVs
        self.max_clusters = max_dcf_clusters
        self.nrestarts = nrestarts
        if rng is None:
            self.rng = np.random.default_rng(102)
        else:
            self.rng = rng


    def fit_segment(self, segment):

        bst = BuildSegmentTree(segment, self.T_CNAs, self.T_SNVs)
        T_Seg, mut_mapping = bst.fit(self.data, segment)
        SegTree = ClonalTree(segment, T_Seg, mut_mapping )
        return SegTree
       

        # for s,cc in zip(snvs,cell_counts_by_snv):
        #     print(f"{mut_lookup[s]}: {cc}")


   
        # for s, id in T_SNV_dict.items():
        #     print(f"{mut_lookup[s]}: {id} ")
       
        #get genotypes tree for these snvs
        
    
        # T_SNVs = {s: self.T_SNVs[s] for s in snvs}
        #get all possible cna trees 
        #get T_SNVs
      
        # clust= DCF_Clustering(self.T_CNAs, self.T_SNVs, T_SNV_Clusters, clusters=4, nrestarts=self.nrestarts, rng=self.rng)
        # #clust results is a 
        # clust_results =clust.fit(tree_id_to_indices, snvs, alt, total)
        # print(f"Segment {segment} obj per snvs: {clust_results.likelihood/len(snvs)}\nDCFs:")
        # print(clust_results.DCFs)
        # clusters = clust_results.get_clusters()
        # T_SNV_dict = clust_results.snv_to_tree(snvs, self.T_SNVs)
            



   







            # print(vaf_by_cn[cn])



    def fit(self, data):
        self.data = data 
        self.segments = data.segments
        # for s in self.segments:
        test_seg= 20
        segTree = self.fit_segment(test_seg)
        print(segTree)
        print(f"Segment {test_seg} complete!")

def convert_tree_string(edge_string,id):
    tree = nx.DiGraph()

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
   

        if geno_v not in node_mapping:
            node_id +=1
            node_mapping[geno_v] = node_id

  
        u_node_id = node_mapping[geno_u]
        v_node_id = node_mapping[geno_v]
        tree.add_edge(u_node_id, v_node_id)
        tree.nodes[v_node_id]["genotype"]= geno_v
        tree.nodes[u_node_id]["genotype"]= geno_u
        index += 2
    
    return GenotypeTree(tree, node_mapping, id=id)
                   
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
    parser.add_argument("-s" ,"--seed", required=False, type=int,
                        help="random number seed (default: 1026)")
    parser.add_argument("--state-trees", required=True, 
                        help= "filename of state trees" )
    tpath = "/scratch/data/leah/pharming/test"
    args = parser.parse_args([
        "-f", f"{tpath}/input/read_counts.tsv",
        "-c", f"{tpath}/input/copy_numbers.tsv",
        "-s", "11",
        "--state-trees", "/scratch/data/leah/pharming/src/test_state_trees.txt"
    ])

    # df = pd.read_csv(f"{tpath}/decifer_out.seg20.csv")
    # print(df.head())
    


    print("\nWelcome to the Pharm! Let's start pharming.....\n")

    snv_trees = read_genotype_trees(args.state_trees)
    # for g in snv_trees:
    #     print(g.id)
    #     print(g)
    T_CNAs = [snv_trees[0].generate_CNA_tree()]
    # for g in snv_trees:
    #     print(g.is_refinement(T_CNA))

    # test_tree = snv_trees[2]
    # for d, cn in zip([0.175, 1],[2,5]):
    #     v = posterior_dcf(d, cn ,5, 10, test_tree)
    #     print(f"cn_sample: {cn} d: {d} v: {v}")

    # for index, row in df.iterrows():
    # # Extract the necessary column values
    #     tree_index = int(row['tree_index'])
    #     cn =5

    #     vaf= row['VAR_1']/row['TOT_1']
    #     dcf = row['point_estimate_DCF1']
    #     T_SNV = snv_trees[tree_index]
    #     v= T_SNV.dcf_to_v(dcf,cn)
    #     d = T_SNV.v_to_dcf(vaf,cn)
    #     print(f"m:{int(row['mut_index'])} cn: {cn}\nobs_vaf: {vaf} est_vaf: {v}\nobs_dcf: {dcf} est_dcf: {d}\n")
    #     if np.abs(d - dcf) > 0.2:
    #         print(T_SNV)
    #         # v_check = T_SNV.dcf_to_v(dcf,cn)
    #         d_check =  T_SNV.v_to_dcf(v,cn)
    #         d_check2 =  T_SNV.v_to_dcf(vaf,cn)
    #         print(f"d_check: {d_check} check2: {d_check2}")
    #         print("check")
    # print("done")
    # Pass the values to your function
  


    col_names = ['chr', 'segment', 'mutation_label', 'cell_label','var', 'total']

    rng = np.random.default_rng(args.seed)
    
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
    Pharming(T_CNAs, snv_trees, rng=rng).fit(dat)



 


 
    
    #create cell indices for 

            #wrangle input data into the correct format
       
        

       
  


        
        # variant_count_data= variant_count_data.set_index(["cell", "mutation"])


 
          
 
        # bin_count_data['cell_index'] = cell_series[bin_count_data['cell']].values
        # bin_count_data = bin_count_data.sort_values(by=['cell_index'])

        # bin_count_data.drop(['cell', 'cell_index'], inplace=True, axis=1)
        # bin_count_data = bin_count_data.to_numpy()
    


    
