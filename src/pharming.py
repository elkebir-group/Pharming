from genotype_tree import CNATree
import networkx as nx 
import numpy as np
from fit_segment_tree import BuildSegmentTree
from data import load_from_files, load_from_pickle
import argparse
from clonal_tree import ClonalTree
import os 
import itertools
import multiprocessing
import cProfile


       
        
class Pharming:
    def __init__(self,seed=1026, max_dcf_clusters=3, start_state=(1,1)) -> None:

        self.rng = np.random.default_rng(seed)
        self.max_clusters = max_dcf_clusters
       

    
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
            seed = self.rng.integers(1e8, size=1)[0]
            try:
                SegTree =BuildSegmentTree(T_CNA,seed, max_clusters=self.max_clusters).fit(self.data, g)

                if SegTree.loglikelihood > best_like:
                    best_like = SegTree.loglikelihood
                    BestSegTree= SegTree
            except:
                print(f"Warning: Segment {g} failed")
               
        
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
           
            #     print(f"Warning: Segment {g} failed!")
            #     self.SegTrees[g]= None
           
        
        self.combine_trees()
        return self.SegTrees
    
    def fit_parallel(self, data, segments= None, num_cores=4):
        
        self.data = data
          
        if segments is None:
            segments = data.segments
      
     
    
        pool = multiprocessing.Pool(processes=num_cores)

        # Map the build_segment_tree function to the segments using the pool
        results = pool.map(self.fit_segment, segments)

        # Close the pool and wait for all processes to complete
        pool.close()
        pool.join()

        self.SegTrees = {g: T_Seg for g,T_Seg in zip(segments,results)}
        
        return self.SegTrees






       










# def my_function():
#      _ = Pharming(args.seed).fit(dat, [0])
#      return 

  

if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")

    parser.add_argument("-f", "--file", required=False,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c" ,"--copy-numbers", required=False,
                        help="input files of copy numbers by segment with unlabeled columns [segment cell totalCN]")
    parser.add_argument("-s" ,"--seed", required=False, type=int,
                        help="random number seed (default: 1026)")
    parser.add_argument("-j" ,"--num-cores", required=False, type=int,
                        help="Max number of cores to use for inferring segment trees")

    parser.add_argument("-g" ,"--segment", required=False, type=int,
                        help="segment id of tree to build")
    parser.add_argument("-o" ,"--out", required=False, type=str,
                        help="directory where output files should be written")
    parser.add_argument("-L", "--likelihoods", type=str,
        help = "filename of marginal likelihood")
    # parser.add_argument("--state-trees", required=False, 
    #                     help= "filename of pregenerated state trees or path of program to generate state trees" )
    
    args = parser.parse_args()
    # parser.add_argument('-g', '--segment', type=int, required=False)
    # parser.add_argument("-d", "--data", type=str)

    # instance = "s13_n1000_m15000_c5_p0.1_l0"
    # tpath = f"/scratch/data/leah/pharming/sim_study/pharming/{instance}"

    # args = parser.parse_args([
    #     # "-f", f"{tpath}/input/read_counts.tsv",
    #     # "-c", f"{tpath}/input/copy_numbers.tsv",
    #     "-d", f"{tpath}/data.pickle",
    #     "-s", "13",
    #     "--segment", "16",
    #     "--out", f"{tpath}/SegTrees",
    #     # "-L", f"{tpath}/like.csv"
    #     # "--state-trees", "/scratch/data/leah/pharming/src/test_state_trees.txt"
    #     # "--state-trees", "/scratch/data/leah/pharming/decifer/build/generatestatetrees"
    # ])


    print("\nWelcome to the Pharm! Let's start pharming.....\n")
    
    if args.data is not None:
        dat = load_from_pickle(args.data)
    elif args.file and args.copy_numbers is not None:
        dat = load_from_files(args.file, args.copy_numbers )
    else:
        IOError("Either both read counts and copy number files must be specified or alternatively, the path to the preprocessed pharming data object!")
    
    if args.segment is None:
        segments = dat.segments
    
    else:
        segments = [args.segment]
    # segments = segments[:4]
    SegTrees = Pharming(args.seed).fit_parallel(dat, segments,num_cores=args.num_cores)
    # SegTrees = Pharming(args.seed).fit(dat, segments)
    # SegTrees ={}
    # cProfile.run("my_function()")

    if not os.path.exists(args.out):
        os.makedirs(args.out)
        print(f"Directory '{args.out}' created successfully.")
    else:
        print(f"Directory '{args.out}' already exists.")
    likelihoods = {}
    for g, T_Seg in SegTrees.items():
        if T_Seg is not None:
            pred_cell, pred_mut = T_Seg.generate_results(dat.cell_lookup, dat.mut_lookup)
            likelihoods[g]= T_Seg.compute_likelihood(dat,g)
            pred_cell.to_csv(f"{args.out}/pred_cell_g{g}.csv", index=False)
            pred_mut.to_csv(f"{args.out}/pred_mut_g{g}.csv",index=False)
            T_Seg.draw(f"{args.out}/tree_g{g}.png")
            T_Seg.save(f"{args.out}/tree_g{g}.pickle")
            T_Seg.save_text(f"{args.out}/tree_g{g}.txt")


    if args.likelihoods is not None:
            print("Saving segment likelihoods...")
            with open(args.likelihoods, "w+") as file:
                file.write("segment,likelihood\n")
                for g, like in likelihoods.items():
                    file.write(f"{g},{like}\n")
     



    print("\nPharming complete...let's go sell our trees at the local Pharmer's Market!")
 


 
    

       
        

       
  


        
 


    
