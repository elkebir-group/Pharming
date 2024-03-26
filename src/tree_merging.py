# from superimposition import Superimposition
import numpy as np
from cna_merge import CNA_Merge
import itertools
from utils import get_top_n, pickle_object, concat_and_sort
import multiprocessing

RANDOM = 'random'
NSNVS = 'nsnvs'
INPLACE = "in place"

class ClonalTreeMerging:
    def __init__(self, k, rng=None, seed=1026, order = INPLACE, progressive=True, top_n=1,
        n_orderings=5, collapse = False, cell_threshold=10, inter_opt=False ):
        
        self.k = k
        if rng is not None:
            self.rng = rng 
        if rng is None:
            self.rng = np.random.default_rng(seed)

        self.top_n = top_n

        self.n_orderings = n_orderings
    
    
        
        if order not in [RANDOM, NSNVS, INPLACE]:
            self.order = RANDOM
        else:
            self.order = order 

        if progressive:
            
            self.merge = self.progressive_merge
 
        else:
            self.merge = self.pairwise_merge
        
        self.collapse = collapse
  
        self.cell_threshold = cell_threshold
        self.segment_failures = set()
        self.inter_opt = inter_opt
   
    
    def fit(self, tree_list, T_m, data, lamb, cores=1):
        '''
            @params tree_list: list of list of ClonalTrees on disjoint subsets of segments
            @params data: a Pharming Data object that is used to fit the superimposed trees 
        '''
        self.data = data 
        self.lamb = lamb
        self.T_m = T_m
        self.cores =cores 
        
        cand_merged_lists = []
        if len(tree_list) <= 0:
                raise ValueError("List must contain at least one tree.")
        
        if self.order == RANDOM:
            for _ in range(self.n_orderings):
                permutated_order = self.rng.permutation(len(tree_list))
                ordered_list = [tree_list[i] for i in permutated_order]
        
        else:
                ordered_list = tree_list

            
            #sort the trees according to other criteria, like number of SNVs or normalized costs
            # pass 
        cand_merged_lists.append(self.merge(ordered_list))
        return  get_top_n(cand_merged_lists, self.top_n)


    def merge_parallel(self, tree1, tree2):
        cnm = CNA_Merge(tree1.get_tree(), tree2.get_tree(), self.T_m.edges, verbose=False)

 
        merged_tree_list = cnm.fit(self.data, self.lamb, self.top_n)
        # if len(merged_tree_list) ==0:
        #     pickle_object(tree1, "test/tree1_fail.pkl")
        #     pickle_object(tree2, "test/tree2_fail.pkl")
        #     pickle_object(cnm, "test/cnm_fail.pkl")
        #     pickle_object(self.data, "tset/cnm_data_fail.pkl")
        # try:
  
        # except:
        #     pickle_object(cnm, "test/cnm.pkl")
        #     pickle_object(self.data, "test/data.pkl")
            # assert False 
        return merged_tree_list

    def merge_helper(self, tree_list1, tree_list2):
            '''
            Assume the resolution of tree_list1 is correct and search through tree_list2
            to find solutions
            @params tree_list1 list of ClonalTrees on the same subset of segments
            @params tree_list2 list of ClonalTrees on the same subset of segments
            '''

            if len(tree_list2) ==0:
                return tree_list1
            
            if len(tree_list1) ==0:
                return tree_list2
            
            tree_list1 = sorted(tree_list1, key=lambda x: x.cost )
            tree_list2 =  sorted(tree_list2, key=lambda x: x.cost )
        
            segs1 = tree_list1[0].segments
            segs2 = tree_list2[0].segments
            
            candidates = []
            sol_list1 = tree_list1[:self.top_n]
            while len(candidates) ==0 and len(tree_list2) > 0:
              
                sol_list2 = tree_list2[:self.top_n]
      
                if len(tree_list2) > self.top_n:
                    tree_list2 = tree_list2[self.top_n:]
                else:
                    tree_list2 = []

                if self.collapse:
                    for sol in sol_list1 + sol_list2:
                        sol.collapse(self.k, self.cell_threshold)
                        

                if self.cores <= 1:
            
                    for sol1, sol2 in itertools.product(sol_list1, sol_list2):
                        merged_tree_list = self.merge_parallel(sol1, sol2)
                        candidates.append(merged_tree_list)

                        # cnm = CNA_Merge(sol1.get_tree(), sol2.get_tree(), self.T_m.edges, verbose=False)
                        # merged_tree_list = cnm.fit(self.data, self.lamb, self.top_n)
                    # for sol in merged_tree_list:
                    #     sol.optimize(self.data, self.lamb)
               
                else:
                    arguments = [(tree1, tree2) for tree1, tree2 in itertools.product(sol_list1, sol_list2)]
                    with multiprocessing.Pool(processes=self.cores) as pool:
                        candidates = pool.starmap(self.merge_parallel, arguments, chunksize=1)
            
                candidates = concat_and_sort(candidates)

            if len(candidates) ==0:
                print(f"Warning, integration failed for segments {segs2}, skipping..")
                self.segment_failures.union(segs2)
                candidates = tree_list1
           
            elif self.inter_opt:
                num = min(self.top_n, candidates)
                # if len(candidates) > self.top_n:
                #     num = self.top_n
                # else:
                #     num = len(candidates)
                
                for i in range(num):
                    candidates[i].optimize(self.data, self.lamb)
            return candidates
      


    def progressive_merge(self, tree_list):
      
        if len(tree_list) == 1:
            return tree_list[0]


        # Initialize the merged_tree with the first two trees in the list
        merged_tree_list = self.merge_helper(tree_list[0], tree_list[1])
        
        # Merge every other tree in the list with the current merged_tree
        for i in range(2, len(tree_list)):
            merged_tree_list=  self.merge_helper(merged_tree_list, tree_list[i])


        return merged_tree_list


    def pairwise_merge(self, tree_list):
        if len(tree_list) == 1:
            # Base case: If there's only one list of trees left, return it
            return tree_list[0]
        # raise NotImplementedError
     
        
        # Create pairs of trees for merging
        pairs = [tree_list[i:i + 2] for i in range(0, len(tree_list), 2)]

        # Initialize a new list for merged trees
        new_tree_list = []

        # Merge pairs of trees and add the merged trees to the new list
        for pair in pairs:
            if len(pair) == 2:
                new_tree_list.append( self.merge_helper(pair[0], pair[1]))
            else:
                # If there's an odd number of trees, add the unpaired tree to the new list
                new_tree_list.append(pair[0])

        # Recursively call merge_trees_in_pairs with the new merged list
        return self.pairwise_merge(new_tree_list)


    


