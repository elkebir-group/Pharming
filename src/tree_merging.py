import numpy as np
from cna_merge import CNA_Merge
import itertools
from utils import get_top_n, pickle_object, concat_and_sort
import multiprocessing
import concurrent.futures
RANDOM = 'random'
NSNVS = 'nsnvs'
INPLACE = "in-place"
WRANDOM = "weighted-random"

class ClonalTreeMerging:
    def __init__(self, k, rng=None, seed=1026, order = INPLACE, progressive=True, top_n=1,
         collapse = False, cell_threshold=10, inter_opt=False ):
        
        self.k = k
        if rng is not None:
            self.rng = rng 
        if rng is None:
            self.rng = np.random.default_rng(seed)

        self.top_n = top_n

    
    
        
        if order not in [RANDOM, NSNVS, INPLACE, WRANDOM]:
            print("Warning: specified order param not valid, using random instead.")
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
        #remove any segments that have no inferrred segment trees 
        tree_list =[lst for lst in tree_list if len(lst) > 0]

        if len(tree_list) ==0:
            print("Warning, no trees to integrate!")
            return cand_merged_lists

        if self.order == RANDOM:
            ordering = self.rng.permutation(len(tree_list))
            
        
        elif self.order == WRANDOM:

            weights = [sol[0].m for sol in tree_list]
            weights = weights/np.sum(weights)
            ordering = self.rng.choice(len(tree_list), size=len(tree_list), replace=False, p=weights)

        elif self.order == NSNVS:
            nsnvs = [sol[0].m for sol in tree_list]
            ordering = np.argsort(nsnvs)[::-1]
            
        else:    
            ordering = [i for i in range(len(tree_list))]
        
        ordered_list = [tree_list[i] for i in ordering]
        cand_merged_lists.append(self.progressive_merge(ordered_list, data, lamb))
        
        return  get_top_n(cand_merged_lists, self.top_n)


    def merge_parallel(self, tree1, tree2, data, lamb):
        try:
            cnm = CNA_Merge(tree1.get_tree(), tree2.get_tree(), self.T_m.edges, verbose=False)

        
            merged_tree_list = cnm.fit(data, lamb, self.top_n)
            print(f"Merged list length: {len(merged_tree_list)}")

        except Exception as e:
        # Log the error or perform other actions
            print(f"An error occurred during tree merging: {e}")
            merged_tree_list = [[]]
   
        return merged_tree_list


    def merge_helper(self, tree_list1, tree_list2, data, lamb):
        '''
        Assume the resolution of tree_list1 is correct and search through tree_list2
        to find solutions
        @params tree_list1 list of ClonalTrees on the same subset of segments
        @params tree_list2 list of ClonalTrees on the same subset of segments
        '''

        if len(tree_list2) == 0:
            return tree_list1
        
        if len(tree_list1) == 0:
            return tree_list2
        
        tree_list1 = sorted(tree_list1, key=lambda x: x.cost)
        tree_list2 = sorted(tree_list2, key=lambda x: x.cost)


        candidates = []
        sol_list1 = list(tree_list1[:self.top_n])

        while len(candidates) == 0 and len(tree_list2) > 0:
            
            sol_list2 = list(tree_list2[:self.top_n])
            
            if len(tree_list2) > self.top_n:
                tree_list2 = list(tree_list2[self.top_n:])
            else:
                tree_list2 = []

            if self.collapse:
                for sol in sol_list1 + sol_list2:
                    sol.collapse(self.k, self.cell_threshold)
                    
            # if self.cores <= 1:
            if True:
                for sol1, sol2 in itertools.product(sol_list1, sol_list2):
                    merged_tree_list = self.merge_parallel(sol1, sol2, data, lamb)
                    candidates.append(merged_tree_list)
                    
            else:
                print(f"Sol list size: {len(sol_list1)}: {len(sol_list2)}")
                arguments = [(tree1, tree2, data, lamb) for tree1, tree2 in itertools.product(sol_list1, sol_list2)]
                with concurrent.futures.ProcessPoolExecutor(max_workers=self.cores) as executor:
                    futures = [executor.submit(self.merge_parallel, *args) for args in arguments]
                    for future in concurrent.futures.as_completed(futures):
                        candidates.append(future.result())

            candidates = concat_and_sort(candidates)

        if len(candidates) == 0:
            candidates = tree_list1
        
        elif self.inter_opt:
            num = min(self.top_n, len(candidates))
            for i in range(num):
                candidates[i].optimize(self.data, self.lamb)

        return candidates

    # def merge_helper(self, tree_list1, tree_list2, data, lamb):
    #         '''
    #         Assume the resolution of tree_list1 is correct and search through tree_list2
    #         to find solutions
    #         @params tree_list1 list of ClonalTrees on the same subset of segments
    #         @params tree_list2 list of ClonalTrees on the same subset of segments
    #         '''

    #         if len(tree_list2) ==0:
    #             return tree_list1
            
    #         if len(tree_list1) ==0:
    #             return tree_list2
            
    #         tree_list1 = sorted(tree_list1, key=lambda x: x.cost )
    #         tree_list2 =  sorted(tree_list2, key=lambda x: x.cost )
        

    #         segs2 = tree_list2[0].segments
            
    #         candidates = []
    #         sol_list1 = list(tree_list1[:self.top_n])
    #         while len(candidates) ==0 and len(tree_list2) > 0:
              
    #             sol_list2 = list(tree_list2[:self.top_n])
      
    #             if len(tree_list2) > self.top_n:
    #                 tree_list2 = list(tree_list2[self.top_n:])
    #             else:
    #                 tree_list2 = []

    #             if self.collapse:
    #                 for sol in sol_list1 + sol_list2:
    #                     sol.collapse(self.k, self.cell_threshold)
                        

    #             if self.cores <= 1:
            
    #                 for sol1, sol2 in itertools.product(sol_list1, sol_list2):
    #                     merged_tree_list = self.merge_parallel(sol1, sol2, data, lamb)
    #                     candidates.append(merged_tree_list)

    #                     # cnm = CNA_Merge(sol1.get_tree(), sol2.get_tree(), self.T_m.edges, verbose=False)
    #                     # merged_tree_list = cnm.fit(self.data, self.lamb, self.top_n)
    #                 # for sol in merged_tree_list:
    #                 #     sol.optimize(self.data, self.lamb)
               
    #             else:
    #                 print(f"Sol list size: {len(sol_list1)}: {len(sol_list2)}")
    #                 arguments = [(tree1, tree2, data, lamb) for tree1, tree2 in itertools.product(sol_list1, sol_list2)]
    #                 with multiprocessing.Pool(processes=self.cores) as pool:
    #                     candidates = pool.starmap(self.merge_parallel, arguments, chunksize=1)
            
    #             candidates = concat_and_sort(candidates)

    #         if len(candidates) ==0:
    #             # print(f"Warning, integration failed for segments {segs2}, skipping..")
    #             # self.segment_failures.union(segs2)
    #             candidates = tree_list1
           
    #         elif self.inter_opt:
    #             num = min(self.top_n, candidates)
    #             # if len(candidates) > self.top_n:
    #             #     num = self.top_n
    #             # else:
    #             #     num = len(candidates)
                
    #             for i in range(num):
    #                 candidates[i].optimize(self.data, self.lamb)
    #         return candidates
      


    def progressive_merge(self, tree_list, data, lamb):
      
        if len(tree_list) == 1:
            return tree_list[0]


        # Initialize the merged_tree with the first two trees in the list
        merged_tree_list = self.merge_helper(tree_list[0], tree_list[1], data, lamb)
        
        # Merge every other tree in the list with the current merged_tree
        for i in range(2, len(tree_list)):
            if i % 5 ==0:
                print(f"Integrating tree list {i} of {len(tree_list)}")
            merged_tree_list=  self.merge_helper(merged_tree_list, tree_list[i],data, lamb)


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


    


