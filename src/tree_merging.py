from superimposition import Superimposition
import numpy as np

class ClonalTreeMerging:
    def __init__(self, rng=None, seed=1026, order = 'random', pairwise=False, top_n=1,
        lamb=0, threshold=10, threads=1, timelimit=100, n_orderings=5 ):
        
        if rng is not None:
            self.rng = rng 
        if rng is None:
            self.rng = np.random.default_rng(seed)

        self.top_n = 1

        #ILP superimposition parameters to use for all problem instances
        self.threshold = threshold 
        self.threads = threads 
        self.timelimit = timelimit 
    
        RANDOM = 'random'
        NSNVS = 'N_SNVS'
        
        if order not in ['random', 'N_SNVS']:
            self.order = RANDOM

        if pairwise:
            self.merge = self.pairwise_merge
        else:
            self.merge = self.progressive_merge
    
    def run(self, tree_list, data):
        '''
            @params tree_list: list of list of ClonalTrees on disjoint subsets of segments
            @params data: a Pharming Data object that is used to fit the superimposed trees 
        '''
        cand_merged_lists = []
        if len(tree_list) <= 0:
                raise ValueError("List must contain at least one tree.")
        if self.order == RANDOM:
            for _ in range(self.n_orderings):
                permutated_list = self.rng.permutation(tree_list)
        
                cand_merged_lists.append(self.merge(permutated_list))
        else:
            #sort the trees according to other criteria, like number of SNVs or normalized costs
            pass 

        flattened_candidates = list(itertools.chain.from_iterable(cand_merged_lists))
        top_n_list = sorted(flattened_candidates, key=lambda x: x.get_cost())[:self.top_n]
        return top_n_list


    def merge_helper(self, tree_list1, tree_list2):
            '''
            @params tree_list1 list of ClonalTrees on the same subset of segments
            @params tree_list2 list of ClonalTrees on the same subset of segments
            '''
            merged_tree_cand = []
            costs =[]
            for tree1, tree2 in itertools.product(tree_list1, tree_list2):

                sp = Superimposition(tree1, tree_list2)
                merged_tree = sp.solve(self.data, lamb= self.lamb, threads=self.threads, timelimit=self.timelimit) 
                merge_tree_cand.append(merge_tree)
        
            merged_tree_list = sorted(object_list, key=lambda x: x.get_cost())[:self.top_n]
            return merged_tree_list

    def progressive_merge(self, tree_list):
        self.data = data 
        if len(tree_list) == 1:
            return tree_list[0]


        # Initialize the merged_tree with the first two trees in the list
        merged_tree_list = self.merge_helper(tree_list[0], tree_list[1])
        
        # Merge every other tree in the list with the current merged_tree
        for i in range(2, len(tree_list)):
            merged_tree_list=  self.merge_helper(merged_tree_listm tree_list[i])


        return merged_tree_list


   def pairwise_merge(self, tree_list):
        
        if len(tree_list) == 1:
            # Base case: If there's only one list of trees left, return it
            return tree_list[0]
        
        # Create pairs of trees for merging
        pairs = [tree_list[i:i + 2] for i in range(0, len(tree_list), 2)]

        # Initialize a new list for merged trees
        new_tree_list = []

        # Merge pairs of trees and add the merged trees to the new list
        for pair in pairs:
            if len(pair) == 2:
                mew_tree_list.append( self.merge_helper(pair[0], pair[1]))
            else:
                # If there's an odd number of trees, add the unpaired tree to the new list
                new_tree_list.append(pair[0])

        # Recursively call merge_trees_in_pairs with the new merged list
        return self.pairwise_merge(new_tree_list)


    


