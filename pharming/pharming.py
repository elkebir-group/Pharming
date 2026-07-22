
# Created by: L.L. Weber
# Created on: 2024-02-29 18:40:59

## Ori's personal commments are made using 2 hashes (## comment)
import itertools
import networkx as nx 
import numpy as np
import multiprocessing

from .sti_v2 import STI
from . import clonelib
from .tree_merging import ClonalTreeMerging
from .utils import concat_and_sort, timeit_decorator, get_top_n
from .dcf_clustering_v2 import DCF_Clustering


class Pharming:
    ## take in our arguments passed in from main.py
    def __init__(self, 
                dcfs=None, 
                k=3,
                start_state=(1,1), 
                seed=102,
                verbose=False, 
                top_n=3, 
                collapse=False,
                order = "weighted-random",
                ninit_segs = None,
                ninit_Tm = None,
                cell_threshold=5,
                thresh_prop = 0,
                sum_condition = True
                  ) -> None:
        
        
        self.verbose =verbose
        self.verbose = True 

        ## introduces an element of randomness
        self.rng = np.random.default_rng(seed)

        if dcfs is not None:
            self.k = len(dcfs)
            self.delta  = dcfs 
        else:
            self.k = k 
            self.delta = None
        
 

      
        self.cnatrees = {} 


        ## number of maternal + paternal alleles starting out with
        self.start_state = start_state
        print(f"Start state: {self.start_state}")
        self.ninit_segs= ninit_segs
        print(f"Max # of initial segments: {self.ninit_segs}")
        self.ninit_Tm = ninit_Tm
        print(f"Max # of full inference mutation cluster trees {self.ninit_Tm}")
        self.top_n = top_n
        print(f"Top n: {self.top_n}")

        ## look at main.py for note about collapse
        self.collapse = collapse
        if cell_threshold is not None:
            self.cell_threshold = cell_threshold
        else:
            cell_threshold = 0

        print(f"Collapse CNA internal nodes: {self.collapse} with cell threshold: {self.cell_threshold}")
        self.order = order
        print(f"Integration ordering: {self.order}")
      
        self.thresh_prop = thresh_prop
        print(f"CN states threshold proportion of cells per segment: {self.thresh_prop}")
        self.sum_condition = sum_condition
        print(f"Filtering mutation cluster trees using the sum condition: {self.sum_condition}")

                
     
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
            return [tree for tree in trees]
        
   
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
            segtrees = concat_and_sort(segtrees)

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
   
    @timeit_decorator
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
      
            print("All specified segment tree inference complete!")
            if init_order is not None:
                segtrees = sorted(segtrees, key=lambda tl: order_dict[list(tl[0].ct.get_segments())[0]])
        
            if init_trees is None:
                tree_list = segtrees
            else:
                tree_list =  [init_trees[i]] + segtrees
            print("Starting integration...")
            top_trees = self.integrate(Tm, tree_list)
            print("Integration complete...")
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
    

    def calculate_cluster_distance(self, profile_a, profile_b):
        """
        Calculates the Euclidean distance between two DCF cluster profiles.
        
        :param np.ndarray or dict profile_a: First DCF profile vector.
        :param np.ndarray or dict profile_b: Second DCF profile vector.
        :return: float representing the distance.
        """
        # If the profiles are passed as dictionaries, convert them to sorted numpy arrays
        if isinstance(profile_a, dict):
            profile_a = np.array([profile_a[i] for i in sorted(profile_a.keys())])
        if isinstance(profile_b, dict):
            profile_b = np.array([profile_b[i] for i in sorted(profile_b.keys())])
            
        return np.linalg.norm(profile_a - profile_b)   
    
    def hierarchical_clustering_dcfs(self, dcfs_runs, max_distance_threshold=0.05):
        k = len(dcfs_runs[0])
        profiles = np.array([[run[i] for i in range(k)] for run in dcfs_runs]) # converts array of dictionary of cluster_num -> dcfs val to array or arrays where the dcfs[run][cluster_num] = dcfs_val
        # Format: (centroid_array, weight_int, list_of_original_indices)
        clusters = [(profiles[i], 1, [i]) for i in range(len(profiles))] 
        
        while len(clusters) > 1:
            min_dist = float('inf')
            to_merge = None
            
            # Find the two closest clusters using our separate function
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    dist = self.calculate_cluster_distance(clusters[i][0], clusters[j][0])
                    if dist < min_dist:
                        min_dist = dist
                        to_merge = (i, j)
                        
            # --- TERMINATION CONDITION ---
            # If the closest clusters are further apart than our threshold, terminate!
            if min_dist > max_distance_threshold:
                print(f"Clustering terminated: Next closest clusters are {min_dist:.4f} apart.")
                break
                
            # Execute the merge if under the threshold
            idx1, idx2 = to_merge
            c1, w1, items1 = clusters[idx1]
            c2, w2, items2 = clusters[idx2]
            
            new_weight = w1 + w2
            new_centroid = (c1 * w1 + c2 * w2) / new_weight
            new_items = items1 + items2
            
            del clusters[idx2]
            del clusters[idx1]
            clusters.append((new_centroid, new_weight, new_items))
            
        # Format output back into DCF dictionary representations
        final_groups = []
        for centroid, weight, items in clusters:
            dcf_dict = {cluster_id: float(val) for cluster_id, val in enumerate(centroid)}
            final_groups.append({"dcf_profile": dcf_dict, "weight": weight, "members": items})
            
        final_groups = sorted(final_groups, key=lambda x: x["weight"], reverse=True)
        return final_groups


    def assemble_tree(self, dcf, normalized_weight,init_segs, init_order, Tm):
        if Tm is not None:
            scriptTm = [Tm]
        else:
            scriptTm = self.enumerate_mutcluster_trees(dcf)
        # delta = self.delta.copy()
        stis_init = self.preprocess(init_segs, dcf)
        init_trees, costs = self.infer(scriptTm, stis_init, init_order=init_order)
        self.clonal_trees = init_trees

        #identify the mutation cluster trees that yield minimum cost over the initial segments
        ## find the lowest cost indices of the clonal trees (find the indices of the best clonal trees)
        sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])
        if self.ninit_Tm is None or len(sorted_indices) <= self.ninit_Tm:
            smallest_indices = sorted_indices
        else:
            smallest_indices = sorted_indices[:self.ninit_Tm]
        
        best_tree_int = get_top_n(self.clonal_trees , self.top_n)

        # calculate how good these trees really are
        print("Best trees after initial integration ")
        best_score_within_dcf = np.inf

        # review this scoring function with prof and chat
        alpha = 1 #placeholder val for weight on cost as opposed to likelihood
        n = .01 # placeholder val as a constant to control the impact of the size of the cluster on the algorithm
        for i,sol in enumerate(best_tree_int):
            likelihood, snv, cna = sol.compute_likelihood(self.data, self.lamb)
            score = alpha * sol.cost + likelihood - normalized_weight * n
            if score < best_score_within_dcf:
                best_score_within_dcf = score
        return {"dcf" : dcf, "score" : best_score_within_dcf, "scriptTm" : scriptTm, "smallest_indices" : smallest_indices, "init_trees" : init_trees}

    def generate_best_delta(self, init_segs, init_order, num_runs=50, max_distance_threshold=0.05, Tm=None):
        dcfs_runs = []
        if self.delta:
            dcfs_runs.append(self.delta.copy())
        else:
            for i in range(num_runs):
                self.delta = self.infer_dcfs()          
                dcfs_runs.append(self.delta.copy())
        dcf_clusterings = self.hierarchical_clustering_dcfs(dcfs_runs, max_distance_threshold)


        
        best_overall_score = np.inf
        # best_dcf = None
        best_tree_features = None
        for dcf_dict in dcf_clusterings:
            dcf = dcf_dict["dcf_profile"]
            weight = dcf_dict["weight"]
            normalized_weight = weight/num_runs
            tree_features = self.assemble_tree(dcf, normalized_weight, init_segs, init_order, Tm=Tm)
            score = tree_features["score"]
            if score < best_overall_score:
                best_overall_score = score
                best_tree_features = tree_features
        return best_tree_features

    @timeit_decorator
    def fit(self, data, lamb=1e3, segments= None, cores=1, Tm=None):
        print(".....NEW FIT.....")
        self.data = data


        self.lamb = lamb 
        self.cores = cores 
        


        if segments is None:
            segments = data.segments

        ## splits our segments into those with at least 1 snv and min_cn_states unique number of (x, y) CNA states (stored in init_segs) and those which have less but still more than 1 (infer_segs). For example, if a segment contained both (1, 1) and (2, 1), that would could have 2 cn_states
        ## place_segs are all those sges that have 1 cn_state but still have snvs
        init_segs, infer_segs, place_segs, no_snvs_segs = self.partition_segments(segments, min_cn_states=2)
        print(f"\nSegment partition:\ninitial segments: {len(init_segs)}\ninference segments: {len(infer_segs)}\nplace segments: {len(place_segs)}\n")
        print(".....NEW FIT.....")
       ## this returns an order for the segments which is random, but is biased towards segments with higher numbers of snvs
        init_order = self.order_segments(init_segs)
        print("Segment integration order:")
        print(init_order)

        ## if we are given the dcfs, then we wil use that, otherwise, we must calculate it ourselves
        print("Plowing the field.... ")
        num_runs = 50
        max_threshold = .05
        features = self.generate_best_delta(
                    init_segs=init_segs, 
                    init_order=init_order, 
                    num_runs=num_runs, 
                    max_distance_threshold=max_threshold, 
                    Tm=Tm
                )
        self.delta = features["dcf"]
        delta = self.delta.copy()
        scriptTm = features["scriptTm"]
        smallest_indices = features["smallest_indices"]
        init_trees = features["init_trees"]

        print("\nWatering the fields.... ")
        if len(infer_segs) > 0:
            # FIX 2: Prevent IndexError by handling Tm being provided (where scriptTm length is 1)
            if Tm is not None:
                init_Tm = scriptTm
                selected_init_trees = init_trees
            else:
                init_Tm = [scriptTm[i] for i in smallest_indices]
                selected_init_trees = [init_trees[i] for i in smallest_indices]
                
            init_order_infer = self.order_segments(infer_segs)
            stis_infer = self.preprocess(infer_segs, delta)
            
            self.clonal_trees, costs = self.infer(
                init_Tm, 
                stis_infer, 
                selected_init_trees, 
                init_order=init_order_infer
            )     

        best_trees =  get_top_n(self.clonal_trees, self.top_n)

        # print(f" tree | cost | snv | cna")
        for i,b in enumerate(best_trees):
            cost, snv, cna = b.compute_likelihood(self.data, self.lamb)
            # print(f"|{i} | {cost} | {snv} | {cna} |")
        
        ## Removes linear chain from each tree with no cells assigned and maps SNVs in segments with only 1 CN states    
        ## This function takes in a Solution list and post-processes the clonal tree of each solution such that linear chains with no cell assignments are removed from the tree. Each clonal tree object is modified in place.
        ## in other words, if we have a node which is not assigned any particular cell and has only 1 child, we remove the obsolete node in the post process
        self.post_process(best_trees)

        ## Place SNVs that occur in segments with only a single copy number state in the clonal tree of each solution.
        # List solutions: a list of solutions iterable segments: an iterable of segments that consist of only 1 copy number state
        self.place_snvs(best_trees, place_segs)

    
    
        print("\nHarvesting....")
        ## optmize every tree by trying to find the best placement of snvs and cells to nodes using power iteration (coordinate descent). The clonal tree structure itself is not changes
        for sol in best_trees:
            sol.optimize(self.data, self.lamb)
        
        # all_best_trees = []
        # all_best_trees.append(best_trees)
        best_trees = sorted(best_trees, key=lambda x: x.cost)
        # print(f" tree | cost | snv | cna")

        ## show the likelihood of every tree
        for i,b in enumerate(best_trees):
            cost,snv, cna = b.compute_likelihood(self.data, self.lamb)
            # print(f"|{i} | {cost} | {snv} | {cna} |")
 
    
        ## save only the top n best trees and prune all of the leaves without any cells assigned to them
        best_trees = get_top_n(best_trees, self.top_n)
        for sol in best_trees:
            sol.prune_leaves(self.k)

        ## returns the optimized version of the best clonal trees

        return  best_trees

    # @timeit_decorator
    # def fit(self, data, lamb=1e3, segments= None, cores=1, Tm=None):
    #     ## data is your CNA file, snv file, and total read counts all wrapped in the Data data structure
    #     ## Tm specified whether we are already given an snv cluster tree to constrain our fitting to or not
    #     '''
    #     @params Data data: the input data (C,A,D) to fit
    #     @params float lamb (float): a regularization parameter the cost function
    #     @params list segments: the list of segment ids to fit
    #     @params int cores: the number of processors to use 

    #     Fits a clonal tree T (with assosciated genotypes) and a mapping phi of cells to clones with minimum cost
    #     for the input data and specified segments. 

    #     returns a list of the top_n Solutions to the Clonal Tree Inference with Copy Number problem (CTICN)
    #     ''' 
    #     ## initilaizations
    #     self.data = data


    #     self.lamb = lamb 
    #     self.cores = cores 
        


    #     if segments is None:
    #         segments = data.segments

    #     ## splits our segments into those with at least 1 snv and min_cn_states unique number of (x, y) CNA states (stored in init_segs) and those which have less but still more than 1 (infer_segs). For example, if a segment contained both (1, 1) and (2, 1), that would could have 2 cn_states
    #     ## place_segs are all those sges that have 1 cn_state but still have snvs
    #     init_segs, infer_segs, place_segs, no_snvs_segs = self.partition_segments(segments, min_cn_states=2)
    #     print(f"\nSegment partition:\ninitial segments: {len(init_segs)}\ninference segments: {len(infer_segs)}\nplace segments: {len(place_segs)}\n")
        
    #    ## this returns an order for the segments which is random, but is biased towards segments with higher numbers of snvs
    #     init_order = self.order_segments(init_segs)
    #     print("Segment integration order:")
    #     print(init_order)

    #     ## if we are given the dcfs, then we wil use that, otherwise, we must calculate it ourselves
    #     print("Plowing the field.... ")
    #     if self.delta is None:
    #         self.delta = self.infer_dcfs()
       
        
    #     ## if we are given an snv cluster tree to structure our clonal tree off of, we make sure every cluster appearing in our tree also appears as a node in our dcfs calculation
    #         ## if we are not given an snv cluster tree to structure our clonal tree off of, we build all possible snv cluster trees from scratch given our dcfs satisfying the sum rule
    #     if Tm is not None:
    #         scriptTm = [Tm]
    #         for T_m in scriptTm:
    #             for n in T_m:
    #                 if n not in self.delta:
    #                     raise ValueError(f" Node {n} does not match a cluster id. \
    #                                      Each node label in the mutation cluster \
    #                                      tree must map to a unique value in [k] ")
    #     else:
    #         scriptTm = self.enumerate_mutcluster_trees(self.delta)


                    
    
  
       
        
    #     all_best_trees = []

    #     delta = self.delta.copy()
 
    #     # while loop < self.max_loops and len(scriptTm) > 0:
    #     print(f"DCFs delta: {delta}")
    #     print(f"Starting mutation cluster trees iteration with {len(scriptTm)} trees...")
        
    #     ## This is the segment tree inference step. For each segment in init_segs, we infer a list of clonal trees for just that segment (look at fig 2C in the paper)
    #     stis_init = self.preprocess(init_segs, delta)
    
    #     print("Planting the seeds.... ")

    #     # I moved this to earlier in the code
    #     # ## this returns an order for the segments which is random, but is biased towards segments with higher numbers of snvs
    #     # init_order = self.order_segments(init_segs)
    #     # print("Segment integration order:")
    #     # print(init_order)


    #     ## This step is equivalent to the merging step when the best fitting clonal trees from each segment are merged 
    #     init_trees, costs = self.infer(scriptTm, stis_init, init_order=init_order)

    #     self.clonal_trees = init_trees
    

    #     # return utils.concat_and_sort(init_trees)

    #     #identify the mutation cluster trees that yield minimum cost over the initial segments
    #     ## find the lowest cost indices of the clonal trees (find the indices of the best clonal trees)
    #     sorted_indices = sorted(range(len(costs)), key=lambda i: costs[i])
    #     if self.ninit_Tm is None or len(sorted_indices) <= self.ninit_Tm:
    #         smallest_indices = sorted_indices
    #     else:
    #         smallest_indices = sorted_indices[:self.ninit_Tm]
    
        
    #     print(f"Best mutation cluster trees:")
    #     ## print out the cluster tree edges (enough for viewer to visualize tree) that match up to the best clonal trees
    #     for i in smallest_indices:
    #         print(f"{i}: {list(scriptTm[i].edges)}")
    #         # if i == self.ground_truth_tm:
    #         #     print("Including the ground truth mutation cluster tree!")
        
    #     ## this is never used aside from here
    #     best_tree_int = get_top_n(self.clonal_trees, self.top_n)
    #     print("Best trees after initial integration ")
    #     for i,sol in enumerate(best_tree_int):
    #         cost, snv, cna = sol.compute_likelihood(self.data, self.lamb)
            
    #     ## UP UNTIL THIS POINT WE FOUND THE BEST CLONAL TREES ON ONLY THE SEGMENTS WITH HIGH AMOUNTS OF COPY NUMBER STATES. NOW WE USE OUR RESULTS AS A HEURISTIC TO BUILD THE CLONAL TREE OF ALL SEGS
    #     ## WE NOW DO THE SAME PROCESS, BUT IN OUR INPUT, WE USE THE SNV CLUSTER TREE LIST, AS WELL AS THE CLONAL TREES FROM OUR HIGH CN STATES SAMPLE AS A BASELINE FOR HOW WE BUILD THE CLONAL TREE 
    #         ## think of this as the tweaking to our skeletal baseline clonal tree
    #     print("\nWatering the fields.... ")
    #     if len(infer_segs) > 0:
    #         init_Tm = [scriptTm[i] for i in smallest_indices]
    #         init_order = self.order_segments(infer_segs)
    #         stis_infer = self.preprocess(infer_segs, delta)
    #         self.clonal_trees, costs = self.infer(init_Tm, stis_infer, 
    #                                                 [init_trees[i] for i in smallest_indices], 
    #                                                 init_order = init_order )
        
    #     best_trees =  get_top_n(self.clonal_trees, self.top_n)

    #     # print(f" tree | cost | snv | cna")
    #     for i,b in enumerate(best_trees):
    #         cost, snv, cna = b.compute_likelihood(self.data, self.lamb)
    #         # print(f"|{i} | {cost} | {snv} | {cna} |")
        
    #     ## Removes linear chain from each tree with no cells assigned and maps SNVs in segments with only 1 CN states    
    #     ## This function takes in a Solution list and post-processes the clonal tree of each solution such that linear chains with no cell assignments are removed from the tree. Each clonal tree object is modified in place.
    #     ## in other words, if we have a node which is not assigned any particular cell and has only 1 child, we remove the obsolete node in the post process
    #     self.post_process(best_trees)

    #     ## Place SNVs that occur in segments with only a single copy number state in the clonal tree of each solution.
    #     # List solutions: a list of solutions iterable segments: an iterable of segments that consist of only 1 copy number state
    #     self.place_snvs(best_trees, place_segs)

    
    
    #     print("\nHarvesting....")
    #     ## optmize every tree by trying to find the best placement of snvs and cells to nodes using power iteration (coordinate descent). The clonal tree structure itself is not changes
    #     for sol in best_trees:
    #         sol.optimize(self.data, self.lamb)
        
    #     all_best_trees.append(best_trees)
    #     best_trees = sorted(best_trees, key=lambda x: x.cost)
    #     # print(f" tree | cost | snv | cna")

    #     ## show the likelihood of every tree
    #     for i,b in enumerate(best_trees):
    #         cost,snv, cna = b.compute_likelihood(self.data, self.lamb)
    #         # print(f"|{i} | {cost} | {snv} | {cna} |")
 
    
    #     ## save only the top n best trees and prune all of the leaves without any cells assigned to them
    #     best_trees = get_top_n(all_best_trees, self.top_n)
    #     for sol in best_trees:
    #         sol.prune_leaves(self.k)

    #     ## returns the optimized version of the best clonal trees

    #     return  best_trees
        


  

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
            
       
        
        