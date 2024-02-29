import networkx as nx 

from itertools import product
from utils import powerset, draw, merge_lists, load_pickled_object 

import numpy as np
from copy import deepcopy
from clonal_tree import ClonalTree
from genotype import genotype



class CNA_Merge:
    def __init__(self, T1, T2, T_m, verbose=False):
        self.CT1 = T1 
        self.CT2 = T2
        self.T1 = T1.get_tree()
        self.T_m = T_m 
        self.T2 = T2.get_tree()
        self.genotypes1 = T1.get_genotypes()
        self.genotypes2 = T2.get_genotypes()
        self.snvs1 = T1.get_all_muts()
        self.snvs2 = T2.get_all_muts()
        self.old_to_new, self.new_to_old = {}, {}
        self.seg_to_snvs = T1.get_seg_to_muts() | T2.get_seg_to_muts()
        self.k = max(self.T_m)
        for i,u in enumerate(self.T2):
            if u not in self.T_m:
                self.old_to_new[u] = self.k + i + 1
                self.new_to_old[self.k+i +1]  =u

        self.T2 = nx.relabel_nodes(self.T2, self.old_to_new)
        self.all_nodes  = set([n for tree in [self.T1, self.T2] for n in tree])
        self.verbose = verbose 
        if self.verbose:
            draw(self.T1, "test/T1.png")
            draw(self.T2, "test/T2.png")


     
    def construct_clonal_tree(self, T, data, lamb):
        if self.verbose:
            draw(T, "test/current_tree.png")
        
        assert all(n in T for n in self.all_nodes)
  
        def get_predecessor(u, T_orig):
            '''
            get the predessor of node u in T that is also in T_orig
            '''
            parent = list(T.predecessors(u))
            if len(parent) ==0:
                return u  #it must the be root 
            else:
                parent = parent[0]
            if parent in T_orig:
                return parent 
            else:
                return get_predecessor(parent, T_orig)
            
        genos = {u: {} for u in T}
        root = [u for u in T if T.in_degree[u]==0][0] 

        for n in nx.dfs_preorder_nodes(T, root):
            
            #it's a mutation cluster node so concatentae genotypes
            if n in self.T1 and n in self.T2:
            
                genos[n] = deepcopy(self.genotypes1[n]) | deepcopy(self.genotypes2[n])
          
            elif n in self.T1:
                genos[n] = deepcopy(self.genotypes1[n])
                parent =get_predecessor(n, self.T2)
                if parent in self.T_m:
                    v = parent 
                else:
                    v = self.new_to_old[parent]
                for j in self.snvs2:
                    
                    geno = self.genotypes2[v][j].to_tuple()
                    genos[n][j] = genotype(*geno)
          
            
            else:  #n is in T2
                parent = get_predecessor(n, self.T1)
                genos[n] = deepcopy(self.genotypes2[self.new_to_old[n]])
                for j in self.snvs1:
                    geno = self.genotypes1[parent][j].to_tuple()
                    genos[n][j] = genotype(*geno)
 
        
        ct = ClonalTree(T, genos, self.seg_to_snvs.copy())
        # missing_snvs = ct.mut_mapping[ct.root]
        # if len(missing_snvs) > 0:
        #     self.CT1.draw("test/CT1.png", segments=[10])
        #     self.CT2.draw("test/CT2.png", segments=[20])
        #     draw(T, "test/tree.png")
        #     for j in missing_snvs:
        #         if j in self.snvs1:
        #             mut_mapping = self.CT1.mut_mapping
        #             index =1
        #         else:
        #             mut_mapping = self.CT2.mut_mapping 
        #             index =2
        #         for n, snvs in mut_mapping.items():
        #             if j in snvs:
        #                 gained =n
        #                 break 
        #         print(f"gained: {gained}")
        #         if index==1:
        #             for n in self.T1:
        #                 print(self.genotypes1[n][j])
        #         else:
        #             for n in self.T2:
        #                 print(self.genotypes2[n][j])
                
                
            
        cost , ca = ct.assign_cells_by_likelihood(data, lamb=lamb)

        return cost, ca, ct 



    def fit(self, data, lamb, top_n=3):
        all_trees = self.enumerate_trees()
        all_results = []
        for tree in all_trees:
            cost, ct, ca = self.construct_clonal_tree(tree, data, lamb)
            all_results.append((cost,ct,ca))
        
        all_results = sorted(all_results, key= lambda x: x[0])
        if top_n >= len(all_results):
            return all_results[:top_n]
        else:
            return all_results

        




    def enumerate_trees(self):
        desc_nodes   = self.get_desc_dict()
        # for key,val in desc_nodes.items():
        #     print(f"{key}: {val}")

        trees_by_node = {u: self.construct_subgraphs(u, desc_nodes)  for u in desc_nodes }
        
        # mynodes  = [8, 6,5,0] 
        # for n in mynodes:
        #     for tree in trees_by_node[n]:
        #         draw(tree, f"test/node_tree{n}.png")
        
        all_tree_lists = [trees_by_node[u] for u in trees_by_node]
        all_trees = []
        for trees in product(*all_tree_lists):
            tree = nx.compose_all(list(trees))
            all_trees.append(tree)
            # draw(tree, "test/tree.png")
        return all_trees
        
    
  
    def get_desc(self, u,D, T):
        ancestors_of_D = set()
        for node in D:
            ancestors_of_D.update(nx.ancestors(T, node))

        descendants_of_u = nx.descendants(T, u)
        descendants_of_u = descendants_of_u.intersection(ancestors_of_D)
        nodes_to_return = set()
        for node in descendants_of_u:
            if all(nx.has_path(T,node, d) for d in D):
                if node not in self.T_m:
                    nodes_to_return.add(node)
        return nodes_to_return

    def construct_subgraphs(self, u, dict):
        if self.verbose:
            if u ==6:
                print(u)
                print(dict[u])
        all_subgraphs  = []
        prev_key = None
        # if u in [6]:
        #     print("here")
        #if the node has no CN states 
        if all(len(cnset)==0 for key, cn_list in dict[u].items() for cnset in cn_list):
            tree = nx.DiGraph()
            tree.add_node(u)
            for v in self.T_m.successors(u):
                tree.add_edge(u,v) 
            return [tree]
        all_keys = set(dict[u].keys())
        for key, cn_lists in dict[u].items():
            all_keys.remove(key)
            tree_list = self.enumerate_partial_trees(u, key, cn_lists)
            if prev_key is None and tree_list is not None:
                all_subgraphs = []
                for tree, par_dict in tree_list:
                    if len(all_keys) > 0:
                        for v in key:
                            if all(w != v for k in all_keys for w in k if k !=key):
                                tree.add_edge(par_dict[v], v)
                    all_subgraphs.append((tree,par_dict))


                # all_subgraphs = [ tree for tree in tree_list]
            elif prev_key is None and tree_list is None:
                t = nx.DiGraph()
                t.add_node(u)
                par_dict ={v: u for v in key}
                if len(all_keys) > 0:
                    for v in key:
                        if all(w != v for k in all_keys for w in k if k !=key):
                            t.add_edge(par_dict[v], v)

                all_subgraphs.append((t, par_dict))
            elif tree_list is not None:
                new_subgraphs =[]
                for tree, par_dict in all_subgraphs:
                    # draw(tree, "test/start.png")
                    for T, par_dict2 in tree_list:
                        # draw(T, "test/subtree.png")
                
                    
                        if key == (-1,-1):
                            tree = nx.compose(tree, T)
                            # draw(tree, "test/partial_tree.png")
                        else:
                            pt = T.copy()
                        
                            new_root= list(pt.successors(u))[0]
                            pt.remove_node(u)
                        
                            tree  = nx.compose(tree, pt)
                            for v in key:
                                tree.add_edge(par_dict[v], new_root)
                            
                                #if the child cn state is not the remaining keys
                                if all(w != v for k in all_keys for w in k if k !=key):
                                    tree.add_edge(par_dict2[v], v)
                                else:
                                    par_dict[v] = par_dict2[v]
                    
                        new_subgraphs.append((tree, par_dict))
                all_subgraphs = new_subgraphs
                    # par_dict = par_dict | par_dict2
            else:
                if key != (-1,-1):

                    new_subgraphs =[]
                    for tree, par_dict in all_subgraphs:
                        # draw(tree, "test/start.png")
                        for v in key:
                            #if the child cn state is not the remaining keys
                            if all(w != v for k in all_keys for w in k if k !=key):
                                tree.add_edge(par_dict[v], v)
                        # draw(tree, "test/partial_tree.png")
                        new_subgraphs.append((tree,par_dict))
                    all_subgraphs = new_subgraphs
            
            if len(all_subgraphs) > 0:
                prev_key = key

        return [tree for tree, _ in all_subgraphs]

    def get_desc_dict(self):

        cna_nodes_T1 = set([n for n in self.T1 if n not in self.T_m])
        cna_nodes_T2 = set([n for n in self.T2 if n not in self.T_m]) 
        desc_nodes = {}  
   
        for u in self.T_m:
            desc_nodes[u] = {}
            children =  self.T_m.successors(u)
            desc_sets = powerset(children)
            desc_sets = sorted(desc_sets, reverse=True, key=lambda x: len(x))
            cna_nodes = {}

            for d in desc_sets:
                if len(d) > 0:
                    cna_nodes1 = self.get_desc(u, d, self.T1)
                    cna_nodes1 = set(cna_nodes1).intersection(cna_nodes_T1)
                    for n in cna_nodes1:
    
                        cna_nodes_T1.remove(n)
                    cna_nodes2 = self.get_desc(u, d, self.T2)  
                    cna_nodes2 = set(cna_nodes2).intersection(cna_nodes_T2)
                    for n in cna_nodes2:
                        
                        cna_nodes_T2.remove(n) 
                
                    desc_nodes[u][d] = [cna_nodes1, cna_nodes2]
                    # if len(cna_nodes1) > 0 and len(cna_nodes2) > 0:
                    #     mult_keys.append((u, (-1,-1)))
        for u in self.T_m:
            cna_nodes1 = set()
            for v in self.T1.successors(u):
                if v not in self.T_m and v in cna_nodes_T1:
                    cna_nodes1.add(v)
                    cna_nodes_T1.remove(v)
            cna_nodes2 = set()
            for v in self.T2.successors(u):
                if v not in self.T_m and v in cna_nodes_T2:
                    cna_nodes2.add(v)
                    cna_nodes_T2.remove(v)
            desc_nodes[u][(-1,-1)] = [cna_nodes1, cna_nodes2]
        
        return desc_nodes


    @staticmethod
    def get_linear_order(nodeset, tree):
        if len(nodeset) <= 1:
            return [n for n in nodeset]
        root = [n for n in tree if tree.in_degree[n]==0][0]
        return [n for n in nx.dfs_preorder_nodes(tree, root) if n in nodeset]

    def enumerate_partial_trees(self, u, key, cn_node_lists):
        if all(len(cn_set) ==0 for cn_set in cn_node_lists):
                return None 
        tree_list = []
        if  key != (-1,-1):
            order_lists = []
            for nodeset, T in zip(cn_node_lists, [self.T1, self.T2]):
                order_lists.append(self.get_linear_order(nodeset, T))

            if all(len(order_list) ==0 for order_list in order_lists):
                return []
            all_permutations = merge_lists(order_lists[0], order_lists[1])
            
            for perm in all_permutations:
                parent = u
                tree = nx.DiGraph()
                tree.add_node(parent)
                for n in perm:
                    tree.add_edge(parent, n)
                    parent =n 
                parent_dict = {v: parent for v in key if v is not None}

                tree_list.append((tree, parent_dict))
        else:
            
            G = nx.DiGraph()
            G.add_node(u)
            for cn_set, tr in zip(cn_node_lists, [self.T1,self.T2]):
            
            
                for n in cn_set:
                    G.add_edge(u,n)
                    descendants = nx.descendants(tr, n)
            
                    descendants.add(n)

                    # Create a subgraph containing the subtree rooted at node `u`
                    subtree = tr.subgraph(descendants)
                    tree = nx.compose(G, subtree)
        
            for n1 in cn_node_lists[0]:
                for n2 in cn_node_lists[1]:
                    G.add_edge(n1,n2)
                    G.add_edge(n2, n1)
            for n1,n2 in G.edges:
                G[n1][n2]["weight"] = 1
                
            new_trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
            for tr in new_trees:
                tree_list.append((tree, {}))
     


        return tree_list
    


# T_m = nx.DiGraph([(0,1), (0,2)])
# T1 = nx.DiGraph([(0,4),  (4,5), (5,1), (5,2), (2,8)])
# T2 = nx.DiGraph([(0,6), (0,10), (6,7), (7,1), (6,2), (1,9)])

# CT1 = ClonalTree(T1, {}, {})
# CT2 = ClonalTree(T2, {}, {})

# # all_trees = CNA_Merge(CT1, CT2, T_m).enumerate_trees()

# pth = "/Users/leah/Documents/Research/projects/Pharming/test"

# instance = "s11_m5000_k25_l7"
# # instance = "s12_m5000_k25_l7"
# folder = "n1000_c0.25_e0" 
# pth1 = f"simulation_study/input"
# from data import Data, load_from_pickle


# data_file = f"{pth1}/{instance}/{folder}/data.pkl"
# dat = load_from_pickle(data_file)


# ct10 = load_pickled_object(f"{pth}/seg10_trees.pkl")
# SOL1= ct10[0]
# # SOL1.png("test/CT1.png", segments=[10])

# ct20= load_pickled_object(f"{pth}/seg20_trees.pkl")
# SOL2 = ct20[0]
# # SOL2.png("test/CT2.png", segments=[20])

# T_m = load_pickled_object("test/T_m.pkl")
# root = [n for n in T_m if T_m.in_degree[n]==0][0]
# T_m.add_edge(SOL1.ct.root, root) 
# draw(T_m, "test/T_m.png")
# # for i,CT in enumerate([CT1, CT2]):
# #     CT.draw("test/")
# # draw(CT1, "test/CT1.png")
# # draw(CT2, "test/CT2.png")
# all_trees = CNA_Merge(SOL1.ct, SOL2.ct, T_m).fit(dat, lamb=1e5)
# for cost, ca, ct in all_trees:
#     ct.draw("test/ct.png", ca, segments=[10,20])






# class TreeMerging:
#     def __init__(self, rng=None, seed=1026, order = 'random', pairwise=False, top_n=1,
#         lamb=0, threshold=10, threads=1, timelimit=100, n_orderings=5 ):
        
#         if rng is not None:
#             self.rng = rng 
#         if rng is None:
#             self.rng = np.random.default_rng(seed)

#         self.top_n = 1

#         #ILP superimposition parameters to use for all problem instances
#         self.threshold = threshold 
#         self.threads = threads 
#         self.timelimit = timelimit 
    
#         RANDOM = 'random'
#         NSNVS = 'N_SNVS'
        
#         if order not in ['random', 'N_SNVS']:
#             self.order = RANDOM

#         if pairwise:
#             self.merge = self.pairwise_merge
#         else:
#             self.merge = self.progressive_merge
    
#     def run(self, tree_list, data):
#         '''
#             @params tree_list: list of list of ClonalTrees on disjoint subsets of segments
#             @params data: a Pharming Data object that is used to fit the superimposed trees 
#         '''
#         cand_merged_lists = []
#         if len(tree_list) <= 0:
#                 raise ValueError("List must contain at least one tree.")
#         if self.order == RANDOM:
#             for _ in range(self.n_orderings):
#                 permutated_list = self.rng.permutation(tree_list)
        
#                 cand_merged_lists.append(self.merge(permutated_list))
#         else:
#             #sort the trees according to other criteria, like number of SNVs or normalized costs
#             pass 

#         flattened_candidates = list(itertools.chain.from_iterable(cand_merged_lists))
#         top_n_list = sorted(flattened_candidates, key=lambda x: x.get_cost())[:self.top_n]
#         return top_n_list


#     # def merge_helper(self, tree_list1, tree_list2):
#     #         '''
#     #         @params tree_list1 list of ClonalTrees on the same subset of segments
#     #         @params tree_list2 list of ClonalTrees on the same subset of segments
#     #         '''
#     #         merged_tree_cand = []
#     #         costs =[]
#     #         for tree1, tree2 in itertools.product(tree_list1, tree_list2):

#     #             sp = Superimposition(tree1, tree_list2)
#     #             merged_tree = sp.solve(self.data, lamb= self.lamb, threads=self.threads, timelimit=self.timelimit) 
#     #             merge_tree_cand.append(merge_tree)
        
#     #         merged_tree_list = sorted(object_list, key=lambda x: x.get_cost())[:self.top_n]
#     #         return merged_tree_list




  



#     def merge(self, tree1, tree2, T_m ):
#         '''
#         @params: tree1: ClonalTree
#         @params: tree2: ClonalTree
#         @params: T_m: nx.digraph mutational cluster tree
#         '''

#         old_to_new, new_to_old  ={}, {}
           
#         t1_nodes = tree1.clones()
#         t2_nodes = tree2.clones()
#         t1_snvs = tree1.get_all_muts()
#         t2_snvs = tree2.get_all_muts()

#         k = max(t1_nodes) 
#         for i,u in enumerate(tree2.clones()):
#             if u not in T_m:
#                 old_to_new[u] = k + i + 1
#                 new_to_old[k+i +1]  
        
#         T2 = nx.relabel_nodes(tree2.tree, old_to_new)
#         T1 = tree1.tree.copy()




#         def powerset(iterable):
#             s = list(iterable)
#             return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

#         G= nx.DiGraph(self.tree1.tree)
        
#         def find_cn_nodes(segtree):
#             pass 

#         #merge genotypes 
#         for u in tree1.preorder():
#             if u in T_m: 
#                children =  T_m.successors(u)
#                desc_sets = powerset(children)
#                desc_sets = sorted(desc_sets, reverse=True, key=lambda x: len(x))
#                for dset in desc_sets:
#                    cn_nodes1 = find_cn_nodes(tree1, u ,dset)
#                    cn_nodes2 = find_cn_nodes(tree2, u, dset)
#                    for v in cn_nodes1:
#                        for w in cn_nodes2:
#                            G.add_edge(v,w)
#                            G.add_edge(w,v)
        
#         for u,v in G.edges:
#             G[u][v]["weight"] = 1
     

#         def generate_genotypes(tree):
#             genotypes = deepcopy(tree1.genotypes )
 
#             root = [n for n in tree if tree.in_degree[n]==0][0]
#             t1_par = root 
#             t2_par = root 
        
#             for u in nx.dfs_preorder_nodes(tree):
#                 if u in t1_nodes and v in t2_nodes:
#                     genotypes[u] = genotypes[u] | tree2.genotypes[u]

#                 elif u in t1_nodes:
#                     for j in t2_snvs:

#                         tup = tree2.genotypes[t2_par][j].to_tuple()
#                         genotypes[u][j] = genotype(*tup)


#             for u in T_m:
#                 for 
#                 genotypes[u]


#         clonal_trees = []
#         trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
#         for tree in trees:
    
     
        
#             genotypes = generate_genotypes(tree)
#             seg_snv_mapping 
#             clonal_trees.append(ClonalTree(tree,genotypes,{self.ell : self.snvs}))
        
#         return clonal_trees

#     def progressive_merge(self, tree_list):
#         self.data = data 
#         if len(tree_list) == 1:
#             return tree_list[0]


#         # Initialize the merged_tree with the first two trees in the list
#         merged_tree_list = self.merge_helper(tree_list[0], tree_list[1])
        
#         # Merge every other tree in the list with the current merged_tree
#         for i in range(2, len(tree_list)):
#             merged_tree_list=  self.merge_helper(merged_tree_listm tree_list[i])


#         return merged_tree_list


#    def pairwise_merge(self, tree_list):
        
#         if len(tree_list) == 1:
#             # Base case: If there's only one list of trees left, return it
#             return tree_list[0]
        
#         # Create pairs of trees for merging
#         pairs = [tree_list[i:i + 2] for i in range(0, len(tree_list), 2)]

#         # Initialize a new list for merged trees
#         new_tree_list = []

#         # Merge pairs of trees and add the merged trees to the new list
#         for pair in pairs:
#             if len(pair) == 2:
#                 mew_tree_list.append( self.merge_helper(pair[0], pair[1]))
#             else:
#                 # If there's an odd number of trees, add the unpaired tree to the new list
#                 new_tree_list.append(pair[0])

#         # Recursively call merge_trees_in_pairs with the new merged list
#         return self.pairwise_merge(new_tree_list)











    




         
    # all_nodes = [u for nodessets in cn_list for u in nodessets]
    #         if (u, (v,w)) in mult_keys and (v,w) != (-1,-1):
               
                
    #             for perm in permutations(all_nodes):
    #                 T= tree.copy()
    #                 for n in perm:
    #                     T.add_edge(parent,n)
    #                     parent = n
    #                 new_tree_set.append(T)
       

# def construct_subgraphs(u, T1, T2, dict):
#     all_subgraphs  = []
#     prev_key = None
#     all_keys = set(dict[u].keys())
#     for key, cn_lists in dict[u].items():
#         all_keys.remove(key)
#         tree_list = enumerate_partial_trees(u, key, cn_lists, T1, T2)
#         if prev_key is None:
#             all_subgraphs = [ tree for tree in tree_list]
#         elif tree_list is not None:
#             new_subgraphs =[]
#             for tree, par_dict in all_subgraphs:
#                 # draw(tree, "test/start.png")
#                 for T, par_dict2 in tree_list:
#                     # draw(T, "test/subtree.png")
               
                 
#                     if key == (-1,-1):
#                         tree = nx.compose(tree, T)
#                         # draw(tree, "test/partial_tree.png")
#                     else:
#                         pt = T.copy()
                       
#                         new_root= list(pt.successors(u))[0]
#                         pt.remove_node(u)
                       
#                         tree  = nx.compose(tree, pt)
#                         for v in key:
#                             tree.add_edge(par_dict[v], new_root)
                          
#                             #if the child cn state is not the remaining keys
#                             if all(w != v for k in all_keys for w in k if k !=key):
#                                 tree.add_edge(par_dict2[v], v)
#                             else:
#                                 par_dict[v] = par_dict2[v]
                   
#                     new_subgraphs.append((tree, par_dict))
#             all_subgraphs = new_subgraphs
#                 # par_dict = par_dict | par_dict2
#         else:
#             new_subgraphs =[]
#             for tree, par_dict in all_subgraphs:
#                 # draw(tree, "test/start.png")
#                 for v in key:

                          
#                     #if the child cn state is not the remaining keys
#                     if all(w != v for k in all_keys for w in k if k !=key):
#                         tree.add_edge(par_dict[v], v)
#                 # draw(tree, "test/partial_tree.png")
#                 new_subgraphs.append((tree,par_dict))
#             all_subgraphs = new_subgraphs
#         prev_key = key

#     return [tree for tree, _ in all_subgraphs]
      
                        

                      



            # new_tree_set = []
            # cn_list = dict[(v,w)]
            # all_nodes = [u for nodessets in cn_list for u in nodessets]
            # if (u, (v,w)) in mult_keys and (v,w) != (-1,-1):
               
                
            #     for perm in permutations(all_nodes):
            #         T= tree.copy()
            #         for n in perm:
            #             T.add_edge(parent,n)
            #             parent = n
            #         new_tree_set.append(T)
       
            # elif (v,w) == (-1,-1):
            #     for cn_set, tr in zip(dict[(v,w)], [T1,T2]):
            #         G = nx.DiGraph()
            #         G.add_node(u)
            #         for n in cn_set:
            #             G.add_edge(u,n)
            #             descendants = nx.descendants(tr, n)
                
            #             descendants.add(n)

            #             # Create a subgraph containing the subtree rooted at node `u`
            #             subtree = tr.subgraph(descendants)
            #             tree = nx.compose(G, subtree)
            #     cn_lists = dict[(v,w)]
            #     for n1 in cn_lists[0]:
            #         for n2 in cn_lists[1]:
            #             G.add_edge(n1,n2)
            #             G.add_edge(n2, n1)
            #     for n1,n2 in G.edges:
            #         G[n1][n2]["weight"] = 1
                
            #     new_trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
            #     for tr in new_trees:
            #         new_tree_set.append(nx.compose(tr, tree))
            # else:
            #     st = nx.DiGraph()
            #     parent = None
            #     for cn in dict[(v,w)]:
            #         if len(cn) > 0:
            #             for n in cn:
            #                 if parent is not None:
                    
            #                     st.add_edge(parent, n)
            #                 parent = n 
            #             for n in (v,w):
            #                 if n is not None:
            #                     st.add_edge(parent,n )

                            



                
                
                    

                

             



    
    
             
