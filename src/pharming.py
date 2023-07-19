from genotype_tree import CNATree
import networkx as nx 
import numpy as np
from fit_segment_tree import BuildSegmentTree
from data import load, Data
import argparse
from clonal_tree import ClonalTree
import os 
import itertools

       
        
class Pharming:
    def __init__(self,seed=1026, max_dcf_clusters=4, nrestarts=25, start_state=(1,1)) -> None:

        self.rng = np.random.default_rng(seed)
        self.max_clusters = max_dcf_clusters
        self.nrestarts = nrestarts

    
        self.start_state = start_state

    @staticmethod
    def enumerate_cna_trees(graph):
        spanning_trees = []
        nodes = graph.nodes()

        if len(nodes) ==1:
            G = nx.DiGraph()
            G.add_nodes_from(graph.nodes)

            return [G]

        # Iterate over all possible subsets of edges
        for edges in itertools.combinations(graph.edges(), len(nodes) - 1):
            # Create a new graph with the selected edges
            spanning_tree = nx.Graph(list(edges))

            # Check if the new graph is a valid spanning tree
            if nx.is_tree(spanning_tree):
                spanning_tree.to_directed()
                # if spanning_tree.in_deg
                digraph = nx.DiGraph()
                for e in spanning_tree.edges:
                    digraph.add_edge(*e)
                spanning_trees.append(digraph)

        return spanning_trees
    
    def get_cna_trees(self, cn_states):
   

        if self.start_state not in cn_states:
                cn_states.append(self.start_state)
        
        node_mapping = {self.start_state: 0}

        node_id =0
        for cn in cn_states:
            if cn != self.start_state:
                node_id += 1
                node_mapping[cn] = node_id

        G = nx.complete_graph(len(cn_states))
        cna_trees  =self.enumerate_cna_trees(G)
        T_CNAs = [CNATree(t, node_mapping, id=i) for i,t in enumerate(cna_trees)]
        return T_CNAs
   
    def fit_segment(self,  g):
        BestSegTree = None
        cn_states =  self.data.cn_states_by_seg(g)
            # cn_states = [(1,1), (3,1), (4,1)]
        T_CNAs = self.get_cna_trees(cn_states)
        best_like = np.NINF
        for T_CNA in T_CNAs:
 
            SegTree =BuildSegmentTree(T_CNA).fit(self.data, g)

            if SegTree.loglikelihood > best_like:
                best_like = SegTree.loglikelihood
                BestSegTree= SegTree
        
        return BestSegTree
        
      

    @staticmethod
    def enumerate_T_CNA(T_SNVs):
        id = 0    
        T_CNAs = []
        TCNA_to_TSNVs = {}
        for t_snv in T_SNVs:
            
            cand = t_snv.generate_CNA_tree()
            inlist = False 
            for ct in T_CNAs:
                if ct.is_identical(cand):
                
                    TCNA_to_TSNVs[ct.id].append(t_snv.id)
                    inlist = True
                    break
                    
            if not inlist:
                cand.set_id(id)
                T_CNAs.append(cand)
                TCNA_to_TSNVs[id] = [t_snv.id]
                id += 1
        return T_CNAs, TCNA_to_TSNVs
            

    def combine_trees(self):
        pass 

    def fit(self, data, segments= None):
        
        self.data = data
          
        if segments is None:
            segments = data.segments
      
        self.SegTrees = {}
    
        for g in segments:
            # try:
                self.SegTrees[g] = self.fit_segment(g)
                print(f"Segment {g} complete!")
            # except:
           
                # print(f"Warning: Segment {g} failed!")
                # self.SegTrees[g]= None
           
        
        self.combine_trees()
        return self.SegTrees






       










     

  

if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c" ,"--copy-numbers", required=True,
                        help="input files of copy numbers by segment with unlabeled columns [segment cell totalCN]")
    parser.add_argument("-s" ,"--seed", required=False, type=int,
                        help="random number seed (default: 1026)")

    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory where output files should be written")
    # parser.add_argument("--state-trees", required=False, 
    #                     help= "filename of pregenerated state trees or path of program to generate state trees" )
    
    parser.add_argument('-g', '--segment', type=int, required=False)
    tpath = "/scratch/data/leah/pharming/test"
    args = parser.parse_args([
        "-f", f"{tpath}/input/read_counts.tsv",
        "-c", f"{tpath}/input/copy_numbers.tsv",
        "-s", "11",
        "--segment", "5",
        "--out", "test/Seg5"
        # "--state-trees", "/scratch/data/leah/pharming/src/test_state_trees.txt"
        # "--state-trees", "/scratch/data/leah/pharming/decifer/build/generatestatetrees"
    ])


    print("\nWelcome to the Pharm! Let's start pharming.....\n")
    

    dat = load(args.file, args.copy_numbers )
    if args.segment is None:
        segments = dat.segments
    
    else:
        segments = [args.segment]
    SegTrees = Pharming(args.seed).fit(dat, segments)
    if not os.path.exists(args.out):
        os.makedirs(args.out)
        print(f"Directory '{args.out}' created successfully.")
    else:
        print(f"Directory '{args.out}' already exists.")
    for g, T_Seg in SegTrees.items():
        if T_Seg is not None:
            pred_cell, pred_mut = T_Seg.generate_results(dat.cell_lookup, dat.mut_lookup)

            print(pred_cell.head())
            print(pred_mut.head())
            pred_cell.to_csv(f"{args.out}/pred_cell.csv", index=False)
            pred_mut.to_csv(f"{args.out}/pred_mut.csv",index=False)
            T_Seg.draw(f"{args.out}/tree.png")
            T_Seg.save(f"{args.out}/tree.pickle")
            T_Seg.save_text(f"{args.out}/tree.txt")






    print("\nPharming complete...let's go sell our trees at the local Pharmer's Market!")
 


 
    
    #create cell indices for 

            #wrangle input data into the correct format
       
        

       
  


        
        # variant_count_data= variant_count_data.set_index(["cell", "mutation"])


 
          
 
        # bin_count_data['cell_index'] = cell_series[bin_count_data['cell']].values
        # bin_count_data = bin_count_data.sort_values(by=['cell_index'])

        # bin_count_data.drop(['cell', 'cell_index'], inplace=True, axis=1)
        # bin_count_data = bin_count_data.to_numpy()
    


    
