import networkx as nx 
from itertools import product
from utils import powerset, draw, merge_lists, pickle_object
from clonal_tree import ClonalTree
from solution import Solution
import utils
from copy import deepcopy
from enumerate import Enumerate
import itertools
import numpy as np

class CNA_Merge:
    def __init__(self, T1, T2, Tm_edges, cell_threshold=0,  maxtrees=5000, verbose=False):
        
 
        self.CT1 = deepcopy(T1) 
        self.CT2 =  deepcopy(T2)

        self.CT1_segments = self.CT1.get_segments()

        self.CT2_segments = self.CT2.get_segments()
        # print("Integrating..")
        # print(self.CT1_segments)
        # print(self.CT2_segments)

        # if self.CT2_segments == {8} and self.CT1_segments == {0,2,6,7,11,14}:
        #     self.CT1.draw("test/s11/badCT1.png", segments = self.CT1_segments)
        for ell in self.CT2_segments:
            assert ell not in self.CT1_segments
        
        self.all_segs = self.CT1_segments.union(self.CT2_segments)

        self.init_cna_trees = {}
        for ell in self.CT1_segments:
            self.init_cna_trees[ell] = set(self.CT1.get_cna_tree(ell).edges)
        
        for ell in self.CT2_segments:
            self.init_cna_trees[ell] = set(self.CT2.get_cna_tree(ell).edges)

        self.T1 = self.CT1.get_tree()
        self.T_m = nx.DiGraph(Tm_edges)
        # draw(self.T_m, "test/Tm.png")
        self.T2 = self.CT2.get_tree()
        # draw(self.T1, "test/T1_pre.png")
        # draw(self.T2, "test/T2_pre.png")
        self.genotypes1 = self.CT1.get_genotypes()
        self.genotypes2 = self.CT2.get_genotypes()
        self.snvs1 = self.CT1.get_all_muts()
        self.snvs2 = self.CT2.get_all_muts()
        self.old_to_new, self.new_to_old = {}, {}
        self.mut_to_segs = self.CT1.mut_to_seg.copy() | self.CT2.mut_to_seg.copy()
        self.rho = self.CT1.rho | self.CT2.rho
        self.k = max(self.T_m)
        self.cell_threshold = cell_threshold
        self.maxtrees = maxtrees


        # assert self.k <= 8
        starting_node = max(self.T1)
        for i,u in enumerate(self.T2):
            if u not in self.T_m:
                self.old_to_new[u] = starting_node+ i + 1
                self.new_to_old[starting_node+i +1]  =u

        self.T2 = nx.relabel_nodes(self.T2, self.old_to_new)
        self.all_nodes  = set([n for tree in [self.T1, self.T2] for n in tree])

        self.verbose = verbose 
        # if self.verbose:
        # draw(self.T1, "test/T1.png")
        # draw(self.T2, "test/T2.png")
        
        # if self.CT2.get_segments() == {19}:
        #     print("Here")
        # if not self.CT1.check_genotypes():
        #     self.CT1.png("test/CT1.png")
        #     print(self.CT1.get_segments())
        
        # if not self.CT2.check_genotypes():
        #     self.CT2.png("test/CT2.png")
        #     print(self.CT2.get_segments())
        # for n in self.T2:
        #     if n > 8 and n in self.T1:
        #         print("here")

    def is_compatible(self):
        Tm_nodes = set(self.T_m.nodes)
        for u in self.T1:
            if u not in Tm_nodes:
                desc_u = nx.descendants(self.T1, u).intersection(Tm_nodes)
                for v in self.T2:
                    if v not in Tm_nodes:
                        desc_v = nx.descendants(self.T2, v).intersection(Tm_nodes)
                        if len(desc_u.intersection(desc_v)) > 0 and not (desc_u <= desc_v or desc_v <= desc_u):
                                print("Warning: Trees are incompatible!")
                                print(f"T1: {self.CT1.get_segments()} T2: {self.CT2.get_segments()}")
                                return False 
        return True





    # @utils.timeit_decorator     
    def construct_clonal_tree(self, T, data, lamb):
        # utils.draw(T, "test/T.png")
        if not all(n in T for n in self.all_nodes):
            for n in self.all_nodes:
                if n not in T:
                    print(f" Node {n} not in T")
            print("ASSERT FAILURE!")
        # if not all(n in T for n in self.all_nodes):
        #     utils.pickle_object(self, "test/cnmerge.pkl")
        #     utils.draw(T, "test/bad_T.png")
        #     utils.draw(self.T1, "test/T1.png")
        #     utils.draw(self.T2, "test/T2.png")
           


        cache = {}

        def get_predecessor(u, T_orig):
            '''
            Get the predecessor of node u in T that is also in T_orig
            '''
            if (u, id(T_orig)) in cache:
                return cache[(u, id(T_orig))]
            
            parent = list(T.predecessors(u))
            while parent:
                parent = parent[0]
                if parent in T_orig:
                    cache[(u, id(T_orig))] = parent
                    return parent
                parent = list(T.predecessors(parent))
            
            result = u  # It must be the root if no predecessor in T_orig is found
            cache[(u, id(T_orig))] = result
            return result
 
            
        genos = {u: {} for u in T}
        root = [u for u in T if T.in_degree[u]==0][0] 

        for n in nx.dfs_preorder_nodes(T, root):
            
            #it's a mutation cluster node so concatentae genotypes
            if n in self.T1 and n in self.T2:
            
                genos[n] = self.genotypes1[n].copy()| self.genotypes2[n].copy()
          
            elif n in self.T1:
                genos[n] = self.genotypes1[n].copy()
                parent =get_predecessor(n, self.T2)
                if parent in self.T_m:
                    v = parent 
                else:
                    v = self.new_to_old[parent]
                    # genos[n].update({j: self.genotypes2[parent][j] for j in self.snvs2})
                for j in self.snvs2:
                    genos[n][j] = self.genotypes2[v][j]
                    
                    # geno = self.genotypes2[v][j].to_tuple()
                    # genos[n][j] = genotype(*geno)
          
            
            else:  #n is in T2
                parent = get_predecessor(n, self.T1)
    
                genos[n] = self.genotypes2[self.new_to_old[n]].copy()
                # genos[n].update({j: self.genotypes1[parent][j] for j in self.snvs1})
                for j in self.snvs1:
                    genos[n][j] =self.genotypes1[parent][j]
                    # geno = self.genotypes1[parent][j].to_tuple()
                    # genos[n][j] = genotype(*geno)
 
        
        ct = ClonalTree(T, genos, utils.inverse_dict(self.mut_to_segs), rho=self.rho)
        # assert ct.check_genotypes()
            # print("fail")
            # gt = ct.get_genotypes()
            # utils.pickle_object(self, "test/cnmerge19.pkl")
            # utils.pickle_object(data, "test/data19.pkl")
            # utils.pickle_object(T, "test/T19.png")
        

        # cost, ca = ct.optimize(data, lamb,cell_threshold=self.cell_threshold )
        cost , ca = ct.assign_cells_by_likelihood(data, lamb=lamb)

        return Solution(cost, ct, ca)


    # @utils.timeit_decorator
    def fit(self, data, lamb, top_n=3):

        all_results = []
        if not self.is_compatible():
            return all_results
        all_trees = self.enumerate_trees()
        if self.verbose:
            print(f"Enumerated {len(all_trees)} integrated clonal trees")
        if len(all_trees) > self.maxtrees:
            rng = np.random.default_rng(self.k)
            rng.shuffle(all_trees)
            all_trees = all_trees[:self.maxtrees]
            print(f"Downsampling to {len(all_trees)} integrated clonal trees")
        for i,tree in enumerate(all_trees):
        
            if i % 250 == 0 and self.verbose:
                print(f"Scoring tree {i} of {len(all_trees)}")

            sol= self.construct_clonal_tree(tree, data, lamb)
       

            #post process
            
       
            all_results.append(sol)

        all_results = sorted(all_results, key= lambda x: x.cost)
        if len(all_results) >= top_n:
            all_results = all_results[:top_n]
        
        for sol in all_results:
            sol.post_process(data, lamb, self.k, cell_threshold = self.cell_threshold)
        if len(all_results) > 0:
            sol = all_results[0]
            # for ell, S in self.init_cna_trees.items():
            #     infS = set(sol.ct.get_cna_tree(ell).edges)
        
            #     if S != infS:
            #         print(f"Inferred CNA tree not consistent with input CNA tree for segment {ell}")
            #         print(S)
            #         print(infS)
            #         print(self.CT1_segments)
            #         print(self.CT2_segments)
            #         # pickle_object(self, "test/s11/badmerge.pkl")
            #         # pickle_object(data, "test/s11/badmerge_data.pkl")
            #         assert S == infS

            return [sol]
        else:
            return []
        # if top_n <= len(all_results):
        #     return all_results[:top_n]
        # else:
        #     return all_results
        
                # self.data = data
        # self.lamb = lamb
        # utils.pickle_object(self, "test/cnmerge.pkl")

            # print(k)
            # utils.draw(tree, "test/int_tree.png")
            # assert all(tree.in_degree[n] <=1 for n in tree)

                # continue
                # draw(tree, "test/tree.png")
                # draw(self.T1, "test/T1.png")
                # draw(self.T2, "test/T2.png")
                # draw(self.T_m, "test/Tm.png")
                # print("Warning: clonal trees are incompatible!")
        
                # dd = self.get_desc_dict()
                # for key, val in dd.items():
                #     print(f"{key}:{val}")
                # foo = self.enumerate_trees()


            # sol.png("test/pre_collapse.png")
            # sol.collapse(5)
            # sol.png("test/post_collapse.png")
        
        

        # if len(all_results) ==0:
        #     print("pause")

        



    # @utils.timeit_decorator
    def enumerate_trees(self):
        desc_nodes   = self.get_desc_dict()
        # if self.verbose:
        #     for key,val in desc_nodes.items():
        #         print(f"{key}: {val}")

        trees_by_node = {u: self.construct_subgraphs(u, desc_nodes)  for u in desc_nodes }
        

        all_tree_lists = [trees_by_node[u] for u in trees_by_node]
        all_trees = []
        for trees in product(*all_tree_lists):
            tree = nx.compose_all(list(trees))
            all_trees.append(tree)

        
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
   
        all_subgraphs  = []
        prev_key = None
        # if u in [2]:
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

            # if u ==1:
            #     if key == (3,2):
            #         print("here")
            #     print(key)
            #     print(prev_key)
            #     for tr, par_dict in all_subgraphs:
            #         utils.draw(tr, "test/tr.png")
            #         print(par_dict)
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
                    # if u == 3:
                    #     draw(tree, "test/start.png")
                    for T, par_dict2 in tree_list:
                        new_par_dict = par_dict.copy()
                        # if u ==3:
                        #     draw(T, "test/subtree.png")
                
                    
                        if key == (-1,-1):
                            tree_comp = nx.compose(tree, T)
                      
                        else:
                            pt = T.copy()
                        
                        
                            new_root= list(pt.successors(u))[0]
                            pt.remove_node(u)

                            # if u ==3:
                            #     draw(pt, "test/pt.png")
                            #     draw(tree, "test/tree.png")
                        
                            tree_comp  = nx.compose(tree, pt)
                            # if u ==3:
                            #     draw(tree_comp, "test/tree_comp.png")
                            #     draw(tree, "test/tree.png")
                    
                            for v in key:
                                tree_comp.add_edge(par_dict[v], new_root)
                            
                                #if the child cn state is not the remaining keys
                                if all(w != v for k in all_keys for w in k if k !=key):
                                    tree_comp.add_edge(par_dict2[v], v)
                                else:
                                    new_par_dict[v] = par_dict2[v]
                    
                        # if u==3:
                        
                        #     draw(tree_comp, "test/subtree.png")
                        new_subgraphs.append((tree_comp, new_par_dict))
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
  

  
    @staticmethod
    def decode(tree, root):
        tr = nx.DiGraph()
        troot = (root, root)
        nodedict = {troot: root}
        for u in nx.dfs_preorder_nodes(tree, source=(root,root)):
            u1, u2 = u
            if u in nodedict:
                par = nodedict[u]
            for v in tree.successors(u):
                v1, v2 = v
                if u1 != v1:
                    tr.add_edge(par, v1)
                    nodedict[v] = v1 
                else:
                    tr.add_edge(par, v2)
                    nodedict[v] = v2
        return tr
                # v1, v2 = v 
                # if u1 != v1 and u2 == v2 and v2 == root:
                #     tr.add_edge(u1, v1)
                
                # if u1 != v1 and u2 == v2 and v2 != root:
                #     tr.add_edge(u2, v1)
                
             
                # if u2 != v2 and u1 == v1 and v1 == root:
                #     tr.add_edge(u2, v2)
                
                # if u2 != v2 and u1 == v1 and v1 != root:
                #     tr.add_edge(u1, v2)
      
    

             


    
    def get_cna_only_subtrees(self, n, perm):
        if n in self.T1:
            tree = self.T1 
        else:
            tree =self.T2 
        cna_only_subtrees = []
        for u in tree.successors(n):
            if u not in perm and u > self.k:
                if all(v not in self.T_m for v in nx.descendants(tree, u)):
              
                    subtree = tree.subgraph(nx.descendants(tree, u) | {u}).copy()
                    subtree.add_edge(n,u)
                    cna_only_subtrees.append(subtree)
        return cna_only_subtrees




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
                    subtree_list = self.get_cna_only_subtrees(n, perm)
                    for subtree in subtree_list:
                        tree = nx.compose(tree, subtree)
                    

                parent_dict = {v: parent for v in key if v is not None}

                
                
                tree_list.append((tree, parent_dict))

        else:
            # if u==2:
            #     print("here")
            subgraph_list = []
     
            for cn_set, tr in zip(cn_node_lists, [self.T1,self.T2]):
               
                G = nx.DiGraph()
                G.add_node(u)
            
                for n in cn_set:
                    G.add_edge(u,n)
                    descendants = nx.descendants(tr, n)
            
                    descendants.add(n)

                    # Create a subgraph containing the subtree rooted at node `u`
                    subtree = tr.subgraph(descendants)
                    G = nx.compose(G, subtree)
                # draw(G, "test/s11/G.png")
                subgraph_list.append(G)
            if len(subgraph_list) > 1:
                # for i,t in enumerate(subgraph_list):
                #     draw(t, f"test/s11/t{i}.png")
                refinements = Enumerate(subgraph_list[0], subgraph_list[1], same_root=True).solve()
                for ref in refinements:
                    # draw(ref, "test/s11/ref.png")
                    tr = self.decode(ref, u)
                    # draw(tr, "test/s11/tr.png")
                    tree_list.append((tr, {}))
            elif len(subgraph_list) ==1:
                    tree_list.append((subgraph_list[0], {}))
            else:
                pass


                
            # G = nx.DiGraph()
            # G.add_node(u)
            # for cn_set, tr in zip(cn_node_lists, [self.T1,self.T2]):
            
            
            #     for n in cn_set:
            #         G.add_edge(u,n)
            #         descendants = nx.descendants(tr, n)
            
            #         descendants.add(n)

            #         # Create a subgraph containing the subtree rooted at node `u`
            #         subtree = tr.subgraph(descendants)
            #         G = nx.compose(G, subtree)
     
            # for n1 in cn_node_lists[0]:
            #     for n2 in cn_node_lists[1]:
            #         G.add_edge(n1,n2)
            #         G.add_edge(n2, n1)
            # for n1,n2 in G.edges:
            #     G[n1][n2]["weight"] = 1
            # if u ==2:
            #     draw(G, "test/s11/G2.png") 
            # new_trees = nx.algorithms.tree.branchings.ArborescenceIterator(G)
            # for tr in new_trees:
            #     if u ==2:
            #         draw(tr, "test/s11/tr.png") 
            #     tree_list.append((tr, {}))
     


        return tree_list
    


#testing 
# instance = "s11_m5000_k25_l7"
# # instance = "s12_m5000_k25_l7"
# folder = "n1000_c0.05_e0" 
# pth = f"simulation_study/input"

# gtpth = "test"
# from data import Data, load_from_pickle
# dat = load_from_pickle( f"{pth}/{instance}/{folder}/data.pkl")

# inst = load_from_pickle(f"{gtpth}/cnmerge.pkl")
# inst.verbose = True
# lamb = 1e3
# foo = inst.fit(dat, lamb)


                            



                
                
                    

                

             



    
    
             
