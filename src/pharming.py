
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




        
class Pharming:
    def __init__(self, dcfs=None, cnatrees=None, k=3, T_m=None,
                  start_state=(1,1), seed=102,verbose=False, top_n=3, ilp=False, collapse=False,
                  cell_threshold=10) -> None:
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
            # utils.pickle_object(self.scriptTm, "test/scriptTm.pkl")

        else:
            
            if not isinstance(T_m,list):
                self.scriptTm = [T_m]
                # all_trees =  self.enumerate_mutcluster_trees()
                # choices = self.rng.choice(len(all_trees), size=10)
                # sample_trees = [all_trees[i] for i in choices]
                # for i,t in enumerate(sample_trees):
               
                #     if len(set(T_m.edges).difference(t.edges))==0:
                #         print(f"sampled tree {i} is ground truth tree")
                #         raise Exception("sampled tree is ground truth tree")
                # self.scriptTm += sample_trees

            
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
        self.collapse = collapse
        self.cell_threshold = cell_threshold

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
        # cn_states, counts = self.data.cn_states_by_seg(ell)

        # if len(self.cnatrees) == 0 or ell not in self.cnatrees:
        #     cnatrees = self.enumerate_cna_trees(cn_states)
        # else:
        #     cnatrees = [self.cnatrees[ell]]
        print(f"Starting segment {ell}...")
        
        all_trees = []
        Tm_edges = list(T_m.edges)
        for st in stis:
            trees = st.fit(Tm_edges, self.data, ell)
            if len(trees) > 0:
                all_trees.append(trees)
            # if len(trees) ==0:
            #     print(f"Segment {ell} failed")
            #     utils.pickle_object(st, f"test/st{ell}.pkl")
            #     utils.pickle_object(Tm_edges, f"test/sti_Tm{Tm_edges}.pkl")
        # for S in cnatrees:
        #     st  = STI(S, Tm_edges, self.delta, lamb1=self.lamb, ilp=self.ilp)
            # trees = st.fit(self.data, ell)
            # else:
          
            
        if len(all_trees) ==0:
            print(f"Segment {ell} failed for {Tm_edges}!")
        else:
            segtrees = utils.concat_and_sort(all_trees)

        print(f"Segment {ell} complete!")

       
        return segtrees

    @utils.timeit_decorator
    def segment_trees_inference(self, Tm, stis):
        """
        nx:DiGraph Tm: a mutation cluster tree
        iterable segments: an iterable of segments for which segtrees should be inferred, given Tm

        return list of lists of segment trees, one list per each segment
        """

         #infer a list of lists of segtrees
        if self.cores  <= 1:
            segtrees = []
            for ell in stis:
                # try:
                    if self.data.num_cn_states(ell,self.start_state) > 1:
                
                        segtrees.append(self.fit_segment(ell, Tm, stis[ell]))
                        
        else:
            arguments = [(ell, Tm, stis[ell]) for ell in stis 
                            if self.data.num_cn_states(ell, self.start_state) > 1]
            with multiprocessing.Pool(processes=self.cores) as pool:                
                segtrees = pool.starmap(self.fit_segment, arguments, chunksize=1)
        
        return segtrees
   
    @utils.timeit_decorator
    def integrate(self, Tm_edges, segtrees):
            print(f"All segment trees constructed, integrating trees...")
            ctm = ClonalTreeMerging(self.k, self.rng, n_orderings=1, top_n=self.top_n, 
                                    collapse=self.collapse, cell_threshold=self.cell_threshold)
            
            Tm = nx.DiGraph(Tm_edges)
            #add the normal clone as the root
            root = [n for n in Tm if Tm.in_degree[n]==0][0]
            Tm.add_edge(max(Tm)+2, root) 
            
            top_trees = ctm.fit(segtrees, Tm, self.data, self.lamb, cores=self.cores)

            # for sol in top_trees:
            #     sol.optimize(self.data, self.lamb)

            return top_trees
    

    def partition_segments(self, segments, ninit_segs=2, min_cn_states=2):
        if ninit_segs < len(segments):
      
            init_segs = sorted([ell for ell in segments if self.data.num_cn_states(ell) >= min_cn_states],
                                reverse=True, key= lambda x: self.data.num_snvs(x))
            if len(init_segs) > ninit_segs:
                init_segs = init_segs[:ninit_segs]
        else:
            init_segs = segments
        
        init_segs = set(init_segs)
        

        
        remaining_segs = set(segments) -  init_segs
        infer_segs = set([ell for ell in remaining_segs if self.data.num_cn_states(ell) > 1])
        place_segs = set([ell for ell in segments if self.data.num_cn_states(ell)==1])
        return init_segs, infer_segs, place_segs
      
    @utils.timeit_decorator
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
        
    @utils.timeit_decorator
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

    # @utils.timeit_decorator
    def preprocess_helper(self, ell):
        #TODO parallelize
        # stis = {ell: [] for ell in seg_list}
        stis =[]
        cn_states, counts = self.data.cn_states_by_seg(ell)

        if len(self.cnatrees) == 0 or ell not in self.cnatrees:
            cnatrees = self.enumerate_cna_trees(cn_states)
        else:
            cnatrees = [self.cnatrees[ell]]

            # Tm_edges = list(T_m.edges)
        for S in cnatrees:
            st =STI(ell, S, self.delta, lamb=self.lamb, ilp=self.ilp)
            st.precompute_costs(self.data)
            stis.append(st)
        return ell, stis
    
    def preprocess(self, seg_list):

        if self.cores <= 1:
            stis = {}
            dictlist = []
            for ell in seg_list:
                _, stis[ell]  = self.preprocess_helper(ell)
           
        else:
            
            with multiprocessing.Pool(self.cores) as pool:
                dictlist = pool.map(self.preprocess_helper, seg_list)
            
            stis = {ell: vals for ell, vals in dictlist}
        return stis
       
        # stis = {k: v for d in dictlist for k, v in d.items()}
        # return stis
        # #TODO parallelize
        # stis = {ell: [] for ell in seg_list}
        # for ell in seg_list:
        #     cn_states, counts = self.data.cn_states_by_seg(ell)

        #     if len(self.cnatrees) == 0 or ell not in self.cnatrees:
        #         cnatrees = self.enumerate_cna_trees(cn_states)
        #     else:
        #         cnatrees = [self.cnatrees[ell]]

        #     # Tm_edges = list(T_m.edges)
        #     for S in cnatrees:
        #         st =STI(ell, S, self.delta, lamb1=self.lamb, ilp=self.ilp)
        #         st.precompute_costs(self.data)
        #         stis[ell].append(st)
        # return stis
   

    @utils.timeit_decorator
    def fit(self, data, lamb=1e3, segments= None, cores=1, ninit_segs=3, ninit_Tm=3):
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

        init_segs, infer_segs, place_segs = self.partition_segments(segments, ninit_segs, min_cn_states=2)
        stis = self.preprocess(init_segs)
        init_trees, costs = self.infer(self.scriptTm, stis)
        self.clonal_trees = init_trees

        # return utils.concat_and_sort(init_trees)

        #identify the mutation cluster trees that yield minimum cost over the initial segments
        sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])

        if len(sorted_indices) > ninit_Tm:
            smallest_indices = sorted_indices[:ninit_Tm]
        else:
            smallest_indices  =sorted_indices
        

        print(f"SMALLEST INDICES: {smallest_indices}")
        for i in smallest_indices:
            print(f"{i}: {list(self.scriptTm[i].edges)}")
        
        if len(infer_segs) > 0:
            init_Tm = [self.scriptTm[i] for i in smallest_indices]
            stis = self.preprocess(infer_segs)
            self.clonal_trees, costs = self.infer(init_Tm, stis, [init_trees[i] for i in smallest_indices] )
        
        best_trees =  utils.get_top_n(self.clonal_trees, self.top_n)

            
        self.post_process(best_trees)
        # pickle_object(best_trees, "test/best_trees.pkl")
        # best_trees = load_pickled_object("test/best_trees.pkl")
        # self.clonal_trees = [[]]

        self.place_snvs(best_trees, place_segs)
        for sol in best_trees:
            sol.optimize(self.data, self.lamb)
        
        return  best_trees
        

        #            self.clonal_trees = []
        # best_costs = []
        
        # # if len(remaining_segs) > 1:
        #     init_Tm = [T for T in self.scriptTm[i] for i in smallest_indices]

    
            # # segtrees = []
            # for i in smallest_indices:
            #     Tm = self.scriptTm[i]
            #     segtrees = [list(init_trees[i])]
            #     segtrees += self.segment_trees_inference(Tm, segments=remaining_segs)
            #     top_trees = self.integrate(Tm.edges, segtrees)

            #     best_cost = min(tr.cost for tr in top_trees)
            #     best_costs.append(best_cost)
            #     self.clonal_trees.append(top_trees)
            #     print(f"Integration complete for tree Tm_{i}: {best_cost}")
        
      
        # #find the most promising mutation cluster trees 
        # costs = []
        # init_trees = []
        # for i,Tm in enumerate(self.scriptTm):
        #     print(f"Starting inference for Tm_{i}")
        #     segtrees = self.segment_trees_inference(Tm, segments=init_segs)
        #     # pickle_object(segtrees, "test/segrees_ilp.pkl")
        #     top_trees = self.integrate(Tm, segtrees)
        #     if len(top_trees) > 0:
        #         best_cost = top_trees[0].cost
        #         print(f"Integration complete for tree Tm_{i}: {best_cost}")
        #         costs.append(best_cost)
        #         init_trees.append(top_trees)
        #     else:
        #         print(f"Integration failed for Tm_{i}")
        #         costs.append(np.Inf)

      

        
         #fit the remaining segment trees for each top Tm and integrate
        # remaining_segs = set(segments) - set(init_segs)
        # remaining_segs = [ell for ell in remaining_segs if self.data.num_cn_states(ell) > 1]
        # self.clonal_trees = []
        # best_costs = []
        
        # if len(remaining_segs) > 1:
    
        #     # segtrees = []
        #     for i in smallest_indices:
        #         Tm = self.scriptTm[i]
        #         segtrees = [list(init_trees[i])]
        #         segtrees += self.segment_trees_inference(Tm, segments=remaining_segs)
        #         top_trees = self.integrate(Tm.edges, segtrees)

        #         best_cost = min(tr.cost for tr in top_trees)
        #         best_costs.append(best_cost)
        #         self.clonal_trees.append(top_trees)
        #         print(f"Integration complete for tree Tm_{i}: {best_cost}")
        #         # pickle_object(clonal_trees, "test/clonal_trees.pkl")
        
        # else:
        #     self.clonal_trees = init_trees
        #     best_costs = []
        





         
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
            sol.optimize(self.data, self.lamb)
            # sol.prune()
       


            # ct = sol.ct 
            # for ell in self.data.segments:
            #     if 
            # snvs = dat_snvs - inf_snvs

        


        #TODO: assign SNVs to mutation cluster nodes only
        
        