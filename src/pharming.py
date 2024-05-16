
# Created by: L.L. Weber
# Created on: 2024-02-29 18:40:59



import networkx as nx 
import numpy as np
from sti_v2 import STI
import itertools
import clonelib
import multiprocessing
from tree_merging import ClonalTreeMerging
import utils
from dcf_clustering_v2 import DCF_Clustering
from copy import deepcopy 


        
class Pharming:
    def __init__(self, dcfs=None, cnatrees=None, k=3,
                  start_state=(1,1), seed=102,
                  verbose=False, top_n=3, ilp=False, collapse=False,
                  order = "weighted-random",
                  ninit_segs = None,
                  ninit_Tm = None,
                  cell_threshold=5,
                  max_loops = 3,
                  thresh_prop = 0
                  ) -> None:
        self.verbose =verbose
        self.verbose = True 
        self.rng = np.random.default_rng(seed)

        if dcfs is not None:
            self.k = len(dcfs)
            self.delta  = dcfs 
        else:
            self.k = k 
            self.delta = None
        
 

        if cnatrees is None:
            self.cnatrees = {} 
        else:
            self.cnatrees = cnatrees

    
        self.start_state = start_state
        self.top_n = top_n
        self.ilp = ilp
        self.collapse = collapse
        self.cell_threshold = cell_threshold
        self.order = order
        self.ninit_segs= ninit_segs
        self.ninit_Tm = ninit_Tm
        self.max_loops = max_loops
        print(f"Max loops: {self.max_loops}")
        self.thresh_prop = thresh_prop

        # if True:
        #     self.ground_truth_tm = np.Inf
   
            #self.enumerate_mutcluster_trees()
            # print(f"\nTotal number of mutation cluster trees: {len(self.scriptTm)}")
            # for i,T in enumerate(self.scriptTm):
            #     if set(T.edges) == set(T_m.edges):
            #         print(f"Found ground truth tree at index {i}!")
            #         self.ground_truth_tm = i
            # utils.pickle_object(self.scriptTm, "test/scriptTm.pkl")

     
                # all_trees =  self.enumerate_mutcluster_trees()
                # choices = self.rng.choice(len(all_trees), size=10)
                # sample_trees = [all_trees[i] for i in choices]
                # for i,t in enumerate(sample_trees):
               
                #     if len(set(T_m.edges).difference(t.edges))==0:
                #         print(f"sampled tree {i} is ground truth tree")
                #         raise Exception("sampled tree is ground truth tree")
                # self.scriptTm += sample_trees

            

            
            # for i, T_m in enumerate(self.scriptTm):
            #     draw(T_m, f"test/Tm_{i}.png")
                
     

    def check_dcfs(self, T, delta):
        # return True

        for u in T:
            desc_dcf = sum( delta[u] for u in  sorted(T.successors(u)))
            if delta[u] < desc_dcf or desc_dcf > 1:
                return False 
        return True 

    def enumerate_mutcluster_trees(self, delta, ntrees=None):
        '''
        Enumerate the set of mutation cluster trees that respect 
        the sum condition of the DCFs
    
        '''

        G = nx.DiGraph()
        G.add_nodes_from([q for q in range(self.k)])
        for u,v in itertools.combinations(range(self.k),2):
            if delta[u] >= delta[v]:
                G.add_edge(u,v, weight=1)
            elif delta[u] < delta[v]:
                G.add_edge(v,u, weight=1)
  
            
        trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)

        if ntrees is None:
            return [tree for tree in trees if self.check_dcfs(tree,delta)]
        
        scriptTm = []
        for tree in trees:
            if self.check_dcfs(tree, delta):
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
        
        # if self.verbose:
        #     print(f"Enumerated CNA trees for {cn_states}")
        #     for T in T_CNAS:
        #         print(T)
        
        return T_CNAS


   

    

        
    def fit_segment(self, ell,T_m, stis):
        """
        Fits a segment tree for the data in segment ell given mutation cluster tree T_m.
        
        This function takes a segment ID `ell` and a mutation cluster tree `T_m`, and 
        fits a segment tree for the data in the specified segment. 

        :param int ell: Segment ID.
        :param nx.DiGraph T_m: Mutation cluster tree with nodes labeled by cluster q in [k].
        :return: A list of the top n solutions.
        :rtype: list of Solutions
        """
        if self.verbose:
            print(f"Starting segment {ell}...")
        
        segtrees = []
        Tm_edges = list(T_m.edges)
        for st in stis:
            trees = st.fit(Tm_edges, self.data, ell)
            if len(trees) > 0:
                segtrees.append(trees)
            
        if len(segtrees) ==0:
            print(f"Segment {ell} failed for {Tm_edges}!")
        else:
            segtrees = utils.concat_and_sort(segtrees)

        if self.verbose:
            print(f"Segment {ell} complete!")

       
        return segtrees


    def segment_trees_inference(self, Tm, stis):
        """
        nx:DiGraph Tm: a mutation cluster tree
        iterable segments: an iterable of segments for which segtrees should be inferred, given Tm

        return list of lists of segment trees, one list per each segment
        """


        if self.cores  <= 1:
            segtrees = []
            for ell in stis:
                # try:
                    if self.data.num_cn_states(ell,self.thresh_prop) > 1:
                
                        segtrees.append(self.fit_segment(ell, Tm, stis[ell]))
                        
        else:
            arguments = [(ell, Tm, stis[ell]) for ell in stis 
                            if self.data.num_cn_states(ell, self.thresh_prop) > 1]
            with multiprocessing.Pool(processes=self.cores) as pool:                
                segtrees = pool.starmap(self.fit_segment, arguments, chunksize=1)
        
        return segtrees
   
    @utils.timeit_decorator
    def integrate(self, Tm_edges, segtrees):
            
            print(f"Starting integration for {len(segtrees)} segments...")
            ctm = ClonalTreeMerging(self.k, self.rng, top_n=self.top_n, order=self.order,
                                    collapse=self.collapse, cell_threshold=self.cell_threshold)
            
            Tm = nx.DiGraph(Tm_edges)

            #add the normal clone as the root
            root = [n for n in Tm if Tm.in_degree[n]==0][0]
            Tm.add_edge(max(Tm)+2, root) 
            
            top_trees = ctm.fit(segtrees, Tm, self.data, self.lamb, cores=self.cores)


            return top_trees
    

    def partition_segments(self, segments,  min_cn_states=2):
        if self.ninit_segs is None or self.ninit_segs > len(segments):
            init_segs = [ell for ell in segments 
                         if self.data.num_cn_states(ell, self.thresh_prop) >= min_cn_states and self.data.num_snvs(ell) > 0]

        else:
      
            init_segs = sorted([ell for ell in segments if 
                                self.data.num_cn_states(ell, self.thresh_prop) >= min_cn_states],
                                reverse=True, key= lambda x: self.data.num_snvs(x))
            if len(init_segs) > self.ninit_segs:
                init_segs = init_segs[:self.ninit_segs]

        
        init_segs = set(init_segs)
        

        
        remaining_segs = set(segments) -  init_segs
        infer_segs = set([ell for ell in remaining_segs if self.data.num_cn_states(ell, self.thresh_prop) > 1 and self.data.num_snvs(ell)> 0])
        place_segs = set([ell for ell in segments if self.data.num_cn_states(ell, self.thresh_prop)==1 and self.data.num_snvs(ell)> 0])

        no_snvs_segs = set([ell for ell in segments if  self.data.num_snvs(ell) ==0])

        
        return init_segs, infer_segs, place_segs, no_snvs_segs
      

    def infer(self, Tm_list, stis, init_trees=None):
        """
        infer a clonal tree for a list of mutation cluster trees
        and a subset of segments.  

        list Tm_list: a list of networkx DiGraphs of mutation cluster trees
        iterable seg_list: a set/list of segments to infer clonal trees for each segment
        list init_trees: a list of lists of clonal trees previoyly inferred on disjoint segments from seg_list
        """
        costs = []
        all_trees = []
        for i,Tm in enumerate(Tm_list):

            segtrees = self.segment_trees_inference(Tm, stis)
            if init_trees is None:
                tree_list = segtrees
            else:
                tree_list = [init_trees[i]] + segtrees
            top_trees = self.integrate(Tm, tree_list )
            if len(top_trees) > 0:
                best_cost = top_trees[0].cost
                costs.append(best_cost)
                all_trees.append(top_trees)
            else:
                print(f"Integration failed for Tm_{i}")
                all_trees.append([])
                costs.append(np.Inf)

        return all_trees, costs 

    def place_snvs(self, solutions, segments):
        """
        Place SNVs that occur in segments with only a single copy number state
        in the clonal tree of each solution.

        list solutions: a list of solutions
        iterable segments: an iterable of segments that consist of only 1 copy number state
        """
        rho = {ell: {} for ell in segments}
        seg_to_snvs = {}
        states = {}
        for ell in segments:
     
            rho[ell] = {}
            cn_states, _ = self.data.cn_states_by_seg(ell)
            state = cn_states.pop()
            states[ell] = state
            T_SNV = nx.DiGraph()
     
      
            T_SNV.add_edge((*state, 0,0), (*state, 1,0))
       
            rho[ell] = {q: [T_SNV] for q in range(self.k)}
            seg_to_snvs[ell] = self.data.seg_to_snvs[ell]

        for sol in solutions:
            sol.ct.assign_genotypes(self.data, sol.phi, rho, seg_to_snvs, states)
            sol.ct.add_rho(rho)

    def place_cnas(self, solutions, segments):


        
        for ell in segments:
           states, counts = self.data.cn_states_by_seg(ell)
           setS = self.enumerate_cna_trees(states)

    # @utils.timeit_decorator
    def preprocess_helper(self, ell, delta):

        stis =[]
        cn_prop = self.data.thresholded_cn_prop(ell, self.thresh_prop, self.start_state)
        cn_states = set(cn_prop.keys())
        if len(self.cnatrees) == 0 or ell not in self.cnatrees:
            cnatrees = self.enumerate_cna_trees(cn_states)
        else:
            cnatrees = [self.cnatrees[ell]]

            # Tm_edges = list(T_m.edges)
        for S in cnatrees:
            st =STI(ell, S, delta, lamb=self.lamb, ilp=self.ilp)
            st.precompute_costs(self.data)
            stis.append(st)
        return ell, stis
    
    def infer_dcfs(self):
       dcf_clust = DCF_Clustering(rng= self.rng, nrestarts=18, cna_restriction=1)
    #    like, dcfs , _, _, _ = dcf_clust.decifer(self.data, np.array([0.179, 0.241, 0.32, 0.424, 0.985]) )
    #    print(like)
       like, dcfs , _, _, _= dcf_clust.run(self.data, k_vals=[self.k], cores=6)
       self.delta = {i: dcfs[i] for i in range(len(dcfs))}

    def preprocess(self, seg_list, delta):

        if self.cores <= 1:
            stis = {}
            dictlist = []
            for ell in seg_list:
                _, stis[ell]  = self.preprocess_helper(ell, delta)
           
        else:
            args = [(ell, delta) for ell in seg_list]
            with multiprocessing.Pool(self.cores) as pool:
                dictlist = pool.starmap(self.preprocess_helper, args)
            
            stis = {ell: vals for ell, vals in dictlist}
        return stis
       
  

    @utils.timeit_decorator
    def fit(self, data, lamb=1e3, segments= None, cores=1, Tm=None):
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
        


        if segments is None:
            segments = data.segments

        init_segs, infer_segs, place_segs, no_snvs_segs = self.partition_segments(segments, min_cn_states=2)
        print(f"\nSegment partition:\ninitial segments: {len(init_segs)}\ninference segments: {len(infer_segs)}\nplace segments: {len(place_segs)}\n")
        
        print("Plowing the field.... ")
        if self.delta is None:
            self.delta = self.infer_dcfs()
        # else:
        #     self.delta  = self.check_for_duplicates()

        if Tm is not None:
            scriptTm = [Tm]
            for T_m in scriptTm:
                for n in T_m:
                    if n not in self.delta:
                        raise ValueError(f" Node {n} does not match a cluster id. \
                                         Each node label in the mutation cluster \
                                         tree must map to a unique value in [k] ")
        else:
            scriptTm = self.enumerate_mutcluster_trees(self.delta)
      
        loop = 0

                    
    
  
       
        
        all_best_trees = []
        allscriptTm = []
        delta = self.delta.copy()
 
        while loop < self.max_loops and len(scriptTm) > 0:
            print(f"Starting mutation cluster trees iteration {loop} with {len(scriptTm)} trees...")
            print(f"DCFs delta: {delta}")
            stis_init = self.preprocess(init_segs, delta)

            allscriptTm += scriptTm
        
    
    
            print("Planting the seeds.... ")
            init_trees, costs = self.infer(scriptTm, stis_init)
            self.clonal_trees = init_trees

            # return utils.concat_and_sort(init_trees)

            #identify the mutation cluster trees that yield minimum cost over the initial segments
            sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])
            if self.ninit_Tm is None or len(sorted_indices) <= self.ninit_Tm:
                smallest_indices = sorted_indices
            else:
                smallest_indices = sorted_indices[:self.ninit_Tm]
        
            
            print(f"Best mutation cluster trees:")
            for i in smallest_indices:
                print(f"{i}: {list(scriptTm[i].edges)}")
                # if i == self.ground_truth_tm:
                #     print("Including the ground truth mutation cluster tree!")
            
            print("\nWatering the fields.... ")
            if len(infer_segs) > 0:
                init_Tm = [scriptTm[i] for i in smallest_indices]
                stis_infer = self.preprocess(infer_segs, delta)
                self.clonal_trees, costs = self.infer(init_Tm, stis_infer, [init_trees[i] for i in smallest_indices] )
            
            best_trees =  utils.get_top_n(self.clonal_trees, self.top_n)

                
            self.post_process(best_trees)
    
            self.place_snvs(best_trees, place_segs)

            

        # self.place_cnas(best_trees, no_snvs_segs)
        
        
            print("\nHarvesting....")
            for sol in best_trees:
                sol.optimize(self.data, self.lamb)
            
            all_best_trees.append(best_trees)
            best_trees = sorted(best_trees, key=lambda x: x.cost)
            #check to see if there are any new mutation cluster trees and repeat
            if len(best_trees) > 0:
                best = best_trees[0]
                delta = best.ct.compute_dcfs(best.phi)
                if len(delta) < self.k:
                    for q in range(self.k):
                        if q not in delta:
                            delta[q] = self.rng.uniform()
           

                newTms = self.enumerate_mutcluster_trees(delta)
                scriptTm = []
                for T in newTms:
                    present = False 
                    for T_prime in allscriptTm:
                        if set(T.edges) == set(T_prime.edges):
                            present = True 
                            break 
                    if not present:
                        scriptTm.append(T)
            loop += 1
                            
            
        
        best_trees = utils.get_top_n(all_best_trees, self.top_n)


        return  best_trees
        


  

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
            sol.optimize(self.data, self.lamb)
            # sol.prune()
       


            # ct = sol.ct 
            # for ell in self.data.segments:
            #     if 
            # snvs = dat_snvs - inf_snvs

        


        #TODO: assign SNVs to mutation cluster nodes only
        
        