
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
import concurrent.futures
import time
import stopit



        
class Pharming:
    def __init__(self, dcfs=None, cnatrees=None, k=3,
                  start_state=(1,1), seed=102,
                  verbose=False, top_n=3, ilp=False, collapse=False,
                  order = "weighted-random",
                  ninit_segs = None,
                  ninit_Tm = None,
                  cell_threshold=5,
                  max_loops = 3,
                  thresh_prop = 0,
                  sum_condition = True,
                  timeout = 20
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
        print(f"Start state: {self.start_state}")
        self.ninit_segs= ninit_segs
        print(f"Max # of initial segments: {self.ninit_segs}")
        self.ninit_Tm = ninit_Tm
        print(f"Max # of full inference mutation cluster trees {self.ninit_Tm}")
        self.max_loops = max_loops
        print(f"Max loops: {self.max_loops}")
        self.top_n = top_n
        print(f"Top n: {self.top_n}")
        self.ilp = ilp
        self.collapse = collapse
        self.cell_threshold = cell_threshold

        print(f"Collapse CNA internal nodes: {self.collapse} with cell threshold: {self.cell_threshold}")
        self.order = order
        print(f"Integration ordering: {self.order}")
      
        self.thresh_prop = thresh_prop
        print(f"CN states threshold proportion of cells per segment: {self.thresh_prop}")
        self.sum_condition = sum_condition
        print(f"Filtering mutation cluster trees using the sum condition: {self.sum_condition}")


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
                
     
    @staticmethod
    def check_dcfs(T, delta):
   

        for u in T:
            desc_dcf = sum( delta[u] for u in  sorted(T.successors(u)))
            if delta[u] < desc_dcf or desc_dcf > 1:
                return False 
        return True 

    def enumerate_mutcluster_trees(self, delta):
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

        if self.sum_condition:
            return [tree for tree in trees if self.check_dcfs(tree,delta)]
        else:
            return trees
        
   
    def enumerate_cna_trees_python(self, cn_states):
        cn_states = [cn for cn in cn_states if cn != self.start_state]


        G = nx.DiGraph()
        G.add_nodes_from(cn_states)
        for u in G.nodes:
            for v in G.nodes:
                if u == v:
                    continue
                if (u[0] == 0 and not v[0] > 0) or (u[1] == 0 and not v[1] > 0):
                    continue
                G.add_edge(u,v, weight=1)

        G.add_node(self.start_state)
        for u in cn_states:
            G.add_edge(self.start_state, u,weight=1)

        cnatrees = nx.algorithms.tree.branchings.ArborescenceIterator(G)

        #check the CNA  tree to make sure we don't get a ressurrection of an allele
        def check_tree(T):
            for u,v in T.edges:
                if u[0] ==0 and v[0] > 0:
                    return False 
                if u[1] ==0 and v[1] > 1:
                    return False
            return True

        scriptS = [S for S in cnatrees if check_tree(S)]


     
        return scriptS
       


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

            print(f" segment | cost | snv | cna")
        
            cost,snv, cna = segtrees[0].compute_likelihood(self.data, self.lamb)
            print(f"|{ell} | {cost} | {snv} | {cna} |")

        if self.verbose:
            print(f"Segment {ell} complete!")

        # best_tree = segtrees[0]
        # best_tree.png(f"test/s11/inf{ell}.png")
        return segtrees


    def segment_trees_inference(self, Tm, stis):
        """
        nx:DiGraph Tm: a mutation cluster tree
        iterable segments: an iterable of segments for which segtrees should be inferred, given Tm

        return list of lists of segment trees, one list per each segment
        """


        if self.cores  <= 1:
        # if True:
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
    def integrate(self, Tm_edges, segtrees, restarts=1):
            
            print(f"Starting integration for {len(segtrees)} segments...")
            # ctm = ClonalTreeMerging(self.k, self.rng, top_n=self.top_n, order=self.order,
            #                         collapse=self.collapse, cell_threshold=self.cell_threshold)
            
            all_trees = []
            for i in range(restarts):
                Tm = nx.DiGraph(Tm_edges)

                #add the normal clone as the root
                root = [n for n in Tm if Tm.in_degree[n]==0][0]
                Tm.add_edge(max(Tm)+2, root) 

                
                ctm = ClonalTreeMerging(self.k, self.rng, top_n=self.top_n, order='in-place',
                                        collapse=self.collapse, cell_threshold=self.cell_threshold)
                
                top_trees = ctm.fit(segtrees, Tm, self.data, self.lamb, cores=self.cores)
                all_trees.extend(top_trees)
            
            return top_trees

    def order_segments(self, segs:set):
            
            segs = [ell for ell in segs]
            weights = [self.data.num_snvs(ell) for ell in segs]
            weights = weights/np.sum(weights)
           
            ordering = self.rng.choice(len(segs), size=len(segs), replace=False, p=weights)   
            seg_order = [segs[i] for i in ordering]
            return  seg_order

    def partition_segments(self, segments,  min_cn_states=2):
        if self.ninit_segs is None or self.ninit_segs > len(segments):
            init_segs = [ell for ell in segments 
                         if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) >= min_cn_states and self.data.num_snvs(ell) > 0]

        else:
      
            init_segs = sorted([ell for ell in segments if 
                                self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) >= min_cn_states],
                                reverse=True, key= lambda x: self.data.num_snvs(x))
            if len(init_segs) > self.ninit_segs:
                init_segs = init_segs[:self.ninit_segs]

        
        init_segs = set(init_segs)
        

        
        remaining_segs = set(segments) -  init_segs
        infer_segs = set([ell for ell in remaining_segs if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False) > 1 and self.data.num_snvs(ell)> 0])
        place_segs = set([ell for ell in segments if self.data.num_cn_states(ell, self.thresh_prop, include_start_state=False)==1 and self.data.num_snvs(ell)> 0])

        no_snvs_segs = set([ell for ell in segments if  self.data.num_snvs(ell) ==0])

        print("Init Segs:")
        print(init_segs)

        print("Infer Segs:")
        print(infer_segs)

        print("Place Segs:")
        print(place_segs)

        print("No Snvs Segs:")
        print(no_snvs_segs)
        return init_segs, infer_segs, place_segs, no_snvs_segs
      

    def infer(self, Tm_list, stis, init_trees=None, init_order=None):
        """
        infer a clonal tree for a list of mutation cluster trees
        and a subset of segments.  

        list Tm_list: a list of networkx DiGraphs of mutation cluster trees
        iterable seg_list: a set/list of segments to infer clonal trees for each segment
        list init_trees: a list of lists of clonal trees previoyly inferred on disjoint segments from seg_list
        """
        costs = []
        all_trees = []

        #if given a list of initial trees, then initialize the integrated tree with the first list of tres in the list
        init_on_first = init_trees is not None

        if init_order is not None:
            order_dict = {ell: i for i,ell in enumerate(init_order)}

        for i,Tm in enumerate(Tm_list):
            print(f"Staring Tm {i}: {list(Tm.edges)}")
            # if set([(2, 0), (2, 1), (3, 2)]) != set(Tm.edges):
            #     continue

            segtrees = self.segment_trees_inference(Tm, stis)
      
            
            if init_order is not None:
                segtrees = sorted(segtrees, key=lambda tl: order_dict[list(tl[0].ct.get_segments())[0]])
        
            if init_trees is None:
                tree_list = segtrees
            else:
                tree_list =  [init_trees[i]] + segtrees

            top_trees = self.integrate(Tm, tree_list)
          
            if len(top_trees) > 0:
                best_cost = top_trees[0].cost
                costs.append(best_cost)
                all_trees.append(top_trees)
            else:
                print(f"Integration failed for Tm_{i}")
                all_trees.append([])
                costs.append(np.Inf)
            print(f"Tm {i} complete!")

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
            cn_prop = self.data.thresholded_cn_prop(ell, self.thresh_prop,
                                                          self.start_state, include_start_state=False)
            assert len(cn_prop) == 1
            for k in cn_prop:
                state= k 
                
            states[ell] = state
            T_SNV = nx.DiGraph()
     

            T_SNV.add_edge((*state, 0,0), (*state, 1,0))
            if state != self.start_state:
                T_SNV.add_edge((*self.start_state, 0,0), (*state, 0,0))
       
            rho[ell] = {q: [T_SNV] for q in range(self.k)}
            seg_to_snvs[ell] = self.data.seg_to_snvs[ell]

        for sol in solutions:
            sol.ct.assign_genotypes(self.data, sol.phi, rho, seg_to_snvs, states, start_state=self.start_state)
            sol.ct.add_rho(rho)
            sol.update_segments()

    def place_cnas(self, solutions, segments):


        
        for ell in segments:
           states, counts = self.data.cn_states_by_seg(ell)
           setS = self.enumerate_cna_trees(states)



    def preprocess_helper(self, ell, delta):
        print(f"Segment {ell}: starting preprocessing..")
  

        stis = []
        cn_prop = self.data.thresholded_cn_prop(ell, self.thresh_prop, self.start_state)
        cn_states = set(cn_prop.keys())

        # try:
        
        if len(self.cnatrees) == 0 or ell not in self.cnatrees:
            cnatrees =  self.enumerate_cna_trees_python(cn_states)
            #  cnatrees =  self.enumerate_cna_trees(cn_states)
        else:
            cnatrees = [self.cnatrees[ell]]
        for S in cnatrees:
            st = STI(ell, S, delta, lamb=self.lamb, prop_thresh=self.thresh_prop)
            st.precompute_costs(self.data)
            stis.append(st)
    

        
        print(f"Segment {ell}: preprocessing complete")
        # except Exception as e:
            # print(f"Segment {ell} timed out enumerating CNA trees, skipping..")

 
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
            print(f"DCFs delta: {delta}")
            print(f"Starting mutation cluster trees iteration {loop} with {len(scriptTm)} trees...")
           
            stis_init = self.preprocess(init_segs, delta)

            allscriptTm += scriptTm
        
    
    
            print("Planting the seeds.... ")
            init_order = self.order_segments(init_segs)
            print("Segment integration order:")
            print(init_order)
            init_trees, costs = self.infer(scriptTm, stis_init, init_order=init_order)

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
            
            best_tree_int = utils.get_top_n(self.clonal_trees, self.top_n)
            print("Best trees after initial integration ")
            for i,sol in enumerate(best_tree_int):
               cost, snv, cna = sol.compute_likelihood(self.data, self.lamb)
              
            #    sol.png(f"test/s11/intSol{i}.png")
            #    print(i)
            #    print("cna trees 7 and 8")
            #    print(sol.ct.get_cna_tree(7).edges)
            #    print(sol.ct.get_cna_tree(8).edges)
               if self.verbose:
                print(f"{i}: {cost}, {snv}, {cna}")
            # return best_trees
            print("\nWatering the fields.... ")
            if len(infer_segs) > 0:
                init_Tm = [scriptTm[i] for i in smallest_indices]
                init_order = self.order_segments(infer_segs)
                stis_infer = self.preprocess(infer_segs, delta)
                self.clonal_trees, costs = self.infer(init_Tm, stis_infer, 
                                                      [init_trees[i] for i in smallest_indices], 
                                                      init_order = init_order )
            
            best_trees =  utils.get_top_n(self.clonal_trees, self.top_n)

            print(f" tree | cost | snv | cna")
            for i,b in enumerate(best_trees):
                cost,snv, cna = b.compute_likelihood(self.data, self.lamb)
                print(f"|{i} | {cost} | {snv} | {cna} |")
                
            self.post_process(best_trees)
    
            self.place_snvs(best_trees, place_segs)

            

        # self.place_cnas(best_trees, no_snvs_segs)
        
        
            print("\nHarvesting....")
            for sol in best_trees:
                sol.optimize(self.data, self.lamb)
            
            all_best_trees.append(best_trees)
            best_trees = sorted(best_trees, key=lambda x: x.cost)
            print(f" tree | cost | snv | cna")
            for i,b in enumerate(best_trees):
                cost,snv, cna = b.compute_likelihood(self.data, self.lamb)
                print(f"|{i} | {cost} | {snv} | {cna} |")
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
        for sol in best_trees:
            sol.prune_leaves(self.k)



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
            
       


            # ct = sol.ct 
            # for ell in self.data.segments:
            #     if 
            # snvs = dat_snvs - inf_snvs

        


        #TODO: assign SNVs to mutation cluster nodes only
        
        