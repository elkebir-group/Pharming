
# Created by: L.L. Weber
# Created on: 2024-02-29 18:40:59



import networkx as nx 
import numpy as np
from sti_v2 import STI
import itertools
import clonelib
import multiprocessing
import cProfile
from tree_merging import ClonalTreeMerging
from utils import get_top_n, pickle_object, load_pickled_object, draw



        
class Pharming:
    def __init__(self, dcfs=None, cnatrees=None, k=3, T_m=None,
                  start_state=(1,1), seed=102,verbose=False, top_n=3, ilp=False) -> None:
        self.verbose =verbose 
        self.rng = np.random.default_rng(seed)

        if dcfs is not None:
            self.k = len(dcfs)
            self.delta  = dcfs 
        else:
            self.k = k 
            self.delta = {q: self.rng.random(size=1) for q in range(self.k)}

        if T_m is None:
            self.scriptTm = self.enumerate_mutcluster_trees()

        else:
            
            if not isinstance(T_m,list):
                self.scriptTm = [T_m]
                all_trees =  self.enumerate_mutcluster_trees()
                choices = self.rng.choice(len(all_trees), size=2)
                sample_trees = [all_trees[i] for i in choices]
                for i,t in enumerate(sample_trees):
               
                    if len(set(T_m.edges).difference(t.edges))==0:
                        print(f"sampled tree {i} is ground truth tree")
                        raise Exception("sampled tree is ground truth tree")
                self.scriptTm += sample_trees

            
            for T_m in self.scriptTm:
                for n in T_m:
                    if n not in self.delta:
                        raise ValueError(f" Node {n} does not match a cluster id. \
                                         Each node label in the mutation cluster \
                                         tree must map to a unique value in [k] ")
            
            # for i, T_m in enumerate(self.scriptTm):
            #     draw(T_m, f"test/Tm_{i}.png")
                
        if cnatrees is None:
            self.cnatrees = {} 
        else:
            self.cnatrees = cnatrees

    
        self.start_state = start_state
        self.top_n = top_n
        self.ilp = ilp

    def check_dcfs(self, T):
        for u in T:
            desc_dcf = sum( self.delta[u] for u in  sorted(T.successors(u)))
            if self.delta[u] < desc_dcf or desc_dcf > 1:
                return False 
        return True 

    def enumerate_mutcluster_trees(self, ntrees=None):
        '''
        Enumerate the set of mutation cluster trees that respect 
        the sum condition of the DCFs
    
        '''

        G = nx.DiGraph()
        G.add_nodes_from([q for q in range(self.k)])
        for u,v in itertools.combinations(range(self.k),2):
            if self.delta[u] > self.delta[v]:
                G.add_edge(u,v, weight=1)
            elif self.delta[u] < self.delta[v]:
                G.add_edge(v,u, weight=1)
            else:
                G.add_edge(u,v, weight=1)
                G.add_edge(v,u, weight=1)
            
        trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)

        if ntrees is None:
            return [tree for tree in trees if self.check_dcfs(tree) ]
        
        scriptTm = []
        for tree in trees:
            if self.check_dcfs(tree):
                scriptTm.append(tree)
            
            if len(scriptTm) == ntrees:
                return scriptTm


        return scriptTm

        
    def enumerate_cna_trees(self, cn_states):
   

        trees = clonelib.get_cna_trees(cn_states, *self.start_state )
    
        def convert_to_CNAtree(tree):
    
            S = nx.DiGraph()
            if len(tree) == 0:

                S.add_node(self.start_state)
            else:
                S.add_edges_from(tree)
       
            return S

        T_CNAS = [convert_to_CNAtree(tree) for tree in trees]
        
        if self.verbose:
            print(f"Enumerated CNA trees for {cn_states}")
            for T in T_CNAS:
                print(T)
        
        return T_CNAS


   

    

        
    def fit_segment(self, ell,T_m):
        """
        Fits a segment tree for the data in segment ell given mutation cluster tree T_m.
        
        This function takes a segment ID `ell` and a mutation cluster tree `T_m`, and 
        fits a segment tree for the data in the specified segment. 

        :param int ell: Segment ID.
        :param nx.DiGraph T_m: Mutation cluster tree with nodes labeled by cluster q in [k].
        :return: A list of the top n solutions.
        :rtype: list of Solutions
        """
        cn_states, counts = self.data.cn_states_by_seg(ell)

        if len(self.cnatrees) == 0 or ell not in self.cnatrees:
            cnatrees = self.enumerate_cna_trees(cn_states)
        else:
            cnatrees = [self.cnatrees[ell]]
        all_trees = []
        Tm_edges = list(T_m.edges)
        for S in cnatrees:
            st  = STI(S, Tm_edges, self.delta, lamb1=self.lamb, ilp=self.ilp)
            trees = st.fit(self.data, ell)
            all_trees.append(trees)
        print(f"Segment {ell} complete!")
        return get_top_n(all_trees, self.top_n)


    def segment_trees_inference(self, Tm, segments):
        """
        nx:DiGraph Tm: a mutation cluster tree
        iterable segments: an iterable of segments for which segtrees should be inferred, given Tm

        return list of lists of segment trees, one list per each segment
        """

         #infer a list of lists of segtrees
        if self.cores  <= 1:
            segtrees = []
            for ell in segments:
                # try:
                    if self.data.num_cn_states(ell,self.start_state) > 1:
                        segtrees.append(self.fit_segment(ell, Tm))
                        
        else:
            arguments = [(ell, Tm) for ell in segments 
                            if self.data.num_cn_states(ell, self.start_state) > 1]
            pool = multiprocessing.Pool(processes=self.cores)
            segtrees = pool.starmap(self.fit_segment, arguments)
        
        return segtrees
   

    def integrate(self, Tm_edges, segtrees):
            print(f"All segment trees constructed, integrating trees...")
            ctm = ClonalTreeMerging(self.rng, n_orderings=1, top_n=self.top_n)
            
            Tm = nx.DiGraph(Tm_edges)
            #add the normal clone as the root
            root = [n for n in Tm if Tm.in_degree[n]==0][0]
            Tm.add_edge(max(Tm)+2, root) 
            
            top_trees = ctm.fit(segtrees, Tm, self.data, self.lamb, cores=self.cores)
            return top_trees

    def fit(self, data, lamb=1e3, segments= None, cores=1, ninit_segs=3, ninit_Tm=1):
        '''
        @params Data data: the input data (C,A,D) to fit
        @params float lamb (float): a regularization parameter the cost function
        @params list segments: the list of segment ids to fit
        @params int cores: the number of processors to use 

        Fits a clonal tree T (with assosciated genotypes) and a mapping phi of cells to clones with minimum cost
        for the input data and specified segments. 

        returns a list of the top_n Solutions to the Clonal Tree Inference with Copy Number problem (CTICN)
        '''
        
        self.data = data
        self.lamb = lamb 
        self.cores = cores 
        
        clonal_trees = []
        best_costs = []

        if segments is None:
            segments = data.segments
        

        
        if ninit_segs < len(segments):
            init_segs = self.data.get_largest_segments( ninit_segs, min_cn_states=3)
        else:
            init_segs = segments
      
      

        #find the most promising mutation cluster trees 
        costs = []
        init_trees = []
        for i,Tm in enumerate(self.scriptTm):
            print(f"Starting inference for Tm_{i}")
            segtrees = self.segment_trees_inference(Tm, segments=init_segs)
            # pickle_object(segtrees, "test/segrees_ilp.pkl")
            top_trees = self.integrate(Tm, segtrees)
            best_cost = top_trees[0].cost
            print(f"Integration complete for tree Tm_{i}: {best_cost}")
            costs.append(best_cost)
            init_trees.append(top_trees)

      
        #identify the mutation cluster trees that yield minimum cost over the initial segments
        sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])

        if len(sorted_indices) > ninit_Tm:
            smallest_indices = sorted_indices[:ninit_Tm]
        else:
            smallest_indices  =sorted_indices
        
        
         #fit the remaining segment trees for each top Tm and integrate
        remaining_segs = set(segments) - set(init_segs)
        remaining_segs = [ell for ell in remaining_segs if self.data.num_cn_states(ell) > 1]
        self.clonal_trees = []
        best_costs = []
        
        if len(remaining_segs) > 1:
    
            # segtrees = []
            for i in smallest_indices:
                Tm = self.scriptTm[i]
                segtrees = [list(init_trees[i])]
                segtrees += self.segment_trees_inference(Tm, segments=remaining_segs)
                top_trees = self.integrate(Tm.edges, segtrees)

                best_cost = min(tr.cost for tr in top_trees)
                best_costs.append(best_cost)
                self.clonal_trees.append(top_trees)
                print(f"Integration complete for tree Tm_{i}: {best_cost}")
                # pickle_object(clonal_trees, "test/clonal_trees.pkl")
        
        else:
            self.clonal_trees = init_trees
            best_costs = []
        
        best_trees =  get_top_n(self.clonal_trees, self.top_n)

            
        self.post_process(best_trees)
        
        return  best_trees, best_costs




         
            # #infer a list of lists of segtrees
            # if cores  <= 1:
            #     segtrees = []
            #     for ell in segments:
            #         # try:
            #             if self.data.num_cn_states(ell,self.start_state) > 1:
            #                 segtrees.append(self.fit_segment(ell, T_m))
                         
            # else:
            #     arguments = [(ell, T_m) for ell in segments 
            #                  if self.data.num_cn_states(ell, self.start_state) > 1]
            #     pool = multiprocessing.Pool(processes=cores)
            #     segtrees = pool.starmap(self.fit_segment, arguments)
            # pickle_object(segtrees, "test/segtrees.pkl")
            
            # segtrees = load_pickled_object("test/segtrees.pkl")

        #     best_cost = min(tr.cost for tr in top_trees)
        #     best_costs.append(best_cost)
        #     clonal_trees.append(top_trees)
        #     print(f"Integration complete for tree Tm_{i}: {best_cost}")
        #     # pickle_object(clonal_trees, "test/clonal_trees.pkl")
        # clonal_trees =  get_top_n(clonal_trees, self.top_n)   

    


  

    def post_process(self, sol_list):
        """
        Removes linear chain from each tree with no cells assigned and
        maps SNVs in segments with only 1 CN states 
        
        This function takes in a Solution list and post-processes the clonal tree of each solution
        such that linear chains with no cell assignments are removed from the tree. Each clonal
        tree object is modified in place. 

        :param list sol_list: a list of Solution objects
 

        """

        for sol in sol_list:
            sol.prune()

            # ct = sol.ct 
            # for ell in self.data.segments:
            #     if 
            # snvs = dat_snvs - inf_snvs

        


        #TODO: assign SNVs to mutation cluster nodes only
        
        