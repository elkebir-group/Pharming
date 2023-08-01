import networkx as nx
import numpy as np
import itertools
# from copy import deepcopy
from scipy.stats import binom, mode 
from clonal_tree import ClonalTree
from sklearn import metrics
from sklearn.cluster import KMeans
import time
import seaborn as sns
import matplotlib.pyplot as plt

DCF_THESHOLD=1e-9
ANC = "ancestral"
CLUST = "clustered"   

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Function '{func.__name__}' took {elapsed_time:.6f} seconds to execute.")
        return result
    return wrapper

def likelihood_function( a,d,y,c, alpha):
        
        # if d ==0:
        #     return 1e-10
        
        vaf = (1/c)*y
        # elif y ==0:

        #     val =  binom.pmf(a,d,alpha)
        #     return val
          
        # else:
            # vaf = np.arange(1, c)/c
            # vaf  = 1/c
            # vaf = 0.5
        adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
        val = binom.logpmf(a,d,adjusted_vaf)

        # val[np.isnan(val)] = np.log(1e-10)
            # return binom.pmf(a,d,adjusted_vaf)
        return val.sum(axis=1)
# like_func = np.vectorize(likelihood_function)



'''
Input: 

 - T_CNA: A networkx CNA tree that represents copy number states in all samples
 - T_SNV: A dictionary (index: SNV) of SNV trees (nx DiGraphs) that are refinements of the CNA tree
 - Data a data object containing the input data

Output:
 - A segment tree: nx:DiGraph that is a refinment of the CNA tree
 - SNV genotypes for the segment tree (edge assignment of each SNV to the tree)
 - Cell clustering for the segment tree (MAP assignment of cell to the segment tree)
 - The Pharming likelihood score the segment tree
'''



class FitSegmentTree:
    def __init__(self, T_CNA, seed=1026, max_clusters=4, silhouette_min=0.6, silhouette_improve=0.2, verbose=False):
        self.verbose =verbose 
        self.T_CNA = T_CNA
        self.max_clusters = max_clusters
        # self.max_clusters = 10
        self.sil_min = silhouette_min
        self.sil_improve = silhouette_improve
        self.min_snv_cluster_size = 10
        T_Seg = nx.DiGraph()
        self.rng= np.random.RandomState(seed)
        self.root = self.T_CNA.root
        if len(self.T_CNA.tree) > 1:

            for u,v in self.T_CNA.tree.edges:
                u_x, u_y= self.T_CNA.tree.nodes[u]["genotype"]
                v_x, v_y= self.T_CNA.tree.nodes[v]["genotype"]
                T_Seg.add_node(u, genotype= (u_x, u_y))
                T_Seg.add_node(v, genotype= (v_x, v_y))
                if u_x + u_y < v_x + v_y:
                    etype = "gain"
                elif u_x + u_y >  v_x + v_y:
                    etype = "loss"

                T_Seg.add_edge(u,v ,event=etype)
        else:
            T_Seg.add_nodes_from(self.T_CNA.tree.nodes)
            for n in T_Seg:
               T_Seg.nodes[n]["genotype"] = self.T_CNA.tree.nodes[n]["genotype"]

        self.T_Seg = T_Seg
        #dictionary with snv as keys and genotype trees as values
        self.T_SNV_Clusters, self.T_SNVs = self.T_CNA.enumerate_snv_trees()  
        if self.verbose:
            for tree in self.T_SNVs:
                    print(f"id:{tree.id}")
                    print(tree)
        self.tree_to_snvs = {}
        self.id_to_tree = {tree.id: tree for  tree in self.T_SNVs}


        # for tree_clust in self.T_SNV_Clusters:
        #     for id in tree_clust:
        #         print(f"id:{id}")
        #         print(self.id_to_tree[id])
                    
        self.likelihood= 0

        self.mut_mapping = {}
        self.cell_mapping = None
        self.mut_loss_mapping = {}
 
        self.alpha = 0.001

        self.T_Seg_List = []
        self.G = nx.DiGraph()
        self.G.add_nodes_from(self.T_CNA.tree)
        self.G.add_edges_from(self.T_CNA.tree.edges)
        self.node_dcfs = {}
        self.genotypes = {n: self.T_CNA.tree.nodes[n]["genotype"] for n in self.T_CNA.tree}

        # self.T_SNV_Clusters = self.group_snv_trees()
        # for s, tree in self.T_SNVs.items():
        #     if tree.id not in self.tree_to_snvs:
        #         self.tree_to_snvs[tree.id] = [s]
        #     else:
        #         self.tree_to_snvs[tree.id].append(s)
     

        # if clusters is not None:
        #     self.snvs_by_cluster = {k: [] for k in self.cluster_ids}
        #     for s,k in zip(snvs, clusters):
        #         self.snvs_by_cluster[k].append(s)
   
 
        



    
        






    # def group_snv_trees(self):
    #     tree_clusters = []
    #     visited_trees = []
    #     # leaves = [n for n in self.T_CNA if self.T_CNA.out_degree[n]==0]
    #     # root = [n for n in self.T_CNA if self.T_CNA.in_degree[n]==0 ][0]
    #     for u_geno,v_geno in self.T_CNA.get_edge_genotypes(cna=True):
    #         trees = []
         
    #         for t in self.T_SNVs:
    #             if t.occurs_within(u_geno, v_geno) and t.id not in visited_trees:
    #                 trees.append(t.id)
    #                 visited_trees.append(t.id)
    #         tree_clusters.append(trees)
    #     for t in self.T_SNVs:
    #         if t.id not in visited_trees:
    #             tree_clusters.append([t.id])
    #     return tree_clusters
            


       
    def tree_to_string(self, T_Seg):
        mystr = ""
        for u in nx.dfs_preorder_nodes(T_Seg, source=self.root):
                for v in T_Seg.neighbors(u):
        
                    etype = T_Seg[u][v]["event"]
                    # if etype == "mutation":
                    #     # dcf = T_Seg[u][v]["dcf"]
                    #     # cluster = T_Seg[u][v]["cluster"]
                    # else:
                    #     dcf = "NA"
                    #     cluster = "NA"
                    mystr+= f"Node {u}: {T_Seg.nodes[u]['genotype']} -> Node {v}: {T_Seg.nodes[v]['genotype']}: {etype}\n"
                    # mystr+= f"{etype}: cluster: {cluster} dcf: {dcf}\n"
        return mystr 
       
    def __str__(self) -> str:
        mystr =f"List with {len(self.T_Seg_List)} Segment Trees\n"
        for T_Seg in self.T_Seg_List:
            mystr += self.tree_to_string(T_Seg)
            mystr += "\n"
        return mystr

    # def order_mutation_edges(self, T_Seg):
    #     # pass 
    #     #collapse edges that are all in the same cluster
    #     #reassign clusters to be in accordance with DCF values
    #     leaves = [n for n in T_Seg if T_Seg.out_degree[n]==0]
   
    #     for l in leaves:
    #         dcf_vals = {}
    #         mut_edges = []
    #         root_to_leaf_path = nx.shortest_path(T_Seg, source=self.root, target=l)
    #         for i, parent in enumerate(root_to_leaf_path):
                
    #             if i >= len(root_to_leaf_path)-1:
    #                 break 
    #             child = root_to_leaf_path[i+1]
           
    #             if T_Seg[parent][child]["event"] == "mutation":
    #                 mut_edges.append((parent, child))
    #                 k = T_Seg[parent][child]["cluster"]
    #                 dcf_vals[k] = T_Seg[parent][child]["dcf"]
 
    #         # Sort the dictionary by value in descending order
    #         sorted_dcfs = dict(sorted(dcf_vals.items(), key=lambda item: item[1], reverse=True))
  
    #         #order the edges by decreasing DCF order and update the assigned cluster and dcf values
    #         for k, edge in zip(sorted_dcfs, mut_edges):
    #             u,v = edge 
    #             T_Seg[u][v]["cluster"]=k 
    #             T_Seg[u][v]["dcf"] = sorted_dcfs[k]
    #             self.mut_mapping[v] = self.snvs_by_cluster[k]
            
    



  



    def print_verb(self, mystr ):
        if self.verbose:
            print(mystr)  


    @staticmethod
    def construct_overlap_graph(snvs, cells_by_snvs):
        G = nx.Graph()
        G.add_nodes_from(snvs)
        for i,j in itertools.combinations(snvs, 2):
            cell_overlap = np.intersect1d(cells_by_snvs[i], cells_by_snvs[j])
            if len(cell_overlap) > 0:
                G.add_edge(i,j)
        return G


    def cluster_snvs(self, clust_snvs, snv_index, alt, total):

        dcfs = np.zeros(shape=(len(clust_snvs)))
        #todo
        for i,s in enumerate(clust_snvs):
            # for j,cn in enumerate(alt):
            
                tree = self.id_to_tree[self.T_SNV_dict[s]]
                
    
                n1, n2 =tree.find_split_pairs()
                n, geno = n1 
                x,y, _ = geno 
                cn = x+y

                a = alt[cn][snv_index[s]]
                d = total[cn][snv_index[s]]

                dcfs[i] = tree.v_to_dcf(a/d,cn)
        zero_indices = np.where(dcfs <= DCF_THESHOLD)[0].tolist()
        indices_to_cluster = np.where(dcfs > DCF_THESHOLD)[0]


        dcfs = dcfs[indices_to_cluster]
        dcfs = dcfs.reshape(-1,1)
        clust_dcfs = np.mean(dcfs,axis=0)
  

        # snv_clusters =  np.full_like(clust_snvs, 0)
        snv_clusters = []
        best_k = 0
        intertia_dict= {}
        sil_score_dict = {}



        if len(np.unique(dcfs)) > 1 or len(zero_indices) < len(clust_snvs):
        
        
            max_score = -1
        
          
            bin_width = 0.025  # Set your desired bin width here
            sns.histplot(dcfs, bins=np.arange(min(dcfs), max(dcfs) + bin_width, bin_width), kde=True)

            # sns.histplot(dcfs, kde=True)

            plt.savefig("test/dcf_histogram.png", dpi=300)
            plt.clf()
            for k in range(1, self.max_clusters+1):
            
                km = KMeans(k, random_state=self.rng)
                
                clusters = km.fit_predict(dcfs)
                intertia_dict[k] =km.inertia_
                unique_values, counts = np.unique(clusters, return_counts=True)
                if np.any(counts < self.min_snv_cluster_size):
                    continue
                # for c in np.unique(clusters):
                #     print(f"{c}: {clusters[clusters==c].shape[0]}")
                cluster_centers = km.cluster_centers_  #nclusters x nfeatures
                try:
                    sil_score = metrics.silhouette_score(dcfs, clusters)
                  
                except:
                    sil_score=0.0
                sil_score_dict[k]=sil_score

                if max_score > 0:
                    sil_improve = 1-((1-sil_score)/(1-max_score))
                else:
                    sil_improve = 1
                
                # print(f"{k}: {sil_score}")
                if k==1 or (sil_score > self.sil_min and sil_score > max_score and sil_improve > self.sil_improve):
                    
                    max_score = sil_score
                    best_k = k
                    clust_dcfs =cluster_centers
                    snv_clusters = clusters
            for key, val in intertia_dict.items():
                print(f"k:{key} intertia: {val}")
            for key, val in sil_score_dict.items():
                print(f"k:{key} sil_score: {val}")
            cluster_map = {j: [] for j in range(best_k)}
            for i,j in zip(indices_to_cluster, snv_clusters):
                # s = 
                cluster_map[j].append(clust_snvs[i])
            
            clust_dcfs = clust_dcfs.reshape(-1)
            clust_snvs = np.array(clust_snvs)
            if len(zero_indices) > 0:
                cluster_map[best_k] = clust_snvs[zero_indices].tolist()
                best_k += 1
                clust_dcfs = np.append(clust_dcfs, 0)
        else:
            best_k = 1
            cluster_map = {0: clust_snvs}
            clust_dcfs= np.array([0.0])

        

        return best_k, cluster_map, clust_dcfs
    
    def optimal_clonal_tree(self, data, g):
        best_like = np.NINF
        # print(self)
        for T_Seg in self.T_Seg_List:
            # print(self.tree_to_string(T_Seg))
            segtree = ClonalTree(g, T_Seg, self.mut_mapping, mutated_copies= self.mutated_copies)
            loglike = segtree.compute_likelihood(data, g, self.alpha,attachment="map")
            # cell_mapping, loglike= self.map_assign_cells(segtree.tree, data,g)
            # segtree.set_cell_mapping(cell_mapping) 
            # loglike = segtree.compute_likelihood(data, g)
            if loglike > best_like:
                best_like = loglike
                opt_segtree = segtree 
      
        opt_segtree.loglikelihood = best_like
        return opt_segtree
    
    def add_to_graph(self, tree_clust, k, snv_clusters, dcfs):
        tree = self.id_to_tree[tree_clust[0]]
        
        sorted_clusters = np.argsort(dcfs)[::-1]
        node_id = np.array(list(snv_clusters.keys()))
        node_id.sort()
        # node_id=node_id[::-1]
        node_order = node_id[sorted_clusters]
        x,y, m = tree.split_node_parent_geno
        u,v, m_start = tree.split_geno
        root_id = self.T_CNA.node_mapping[(x,y)]
        if tree.is_leaf(tree.split_node):
            # x,y, m = tree.split_node_parent_geno
            # u,v, m_start = tree.split_node_geno
        
            # init_node_to_dcfs = {n: self.node_dcfs[n] for n in node_order }
            for n in node_order:
                self.G.add_edge(root_id, n)
                self.genotypes[n] = (u,v)
                for c in node_order:
                 
                    
                    if c != n and self.node_dcfs[n] >= self.node_dcfs[c]:
                        self.G.add_edge(n,c)
                
       
            for child in node_order:
   
                for parent in self.node_dcfs:
                    if parent in node_order:
                            continue
                    if parent != child and self.node_dcfs[parent] >= self.node_dcfs[child] and self.genotypes[parent]==self.genotypes[child]:
                        self.G.add_edge(parent, child)
               
        else:
            for n in node_order:

                self.genotypes[n] = (u,v)

            for i in range(1,len(node_order)):
                self.G.add_edge(node_order[i-1], node_order[i])
            
            cna_node = list(tree.tree.neighbors(tree.split_node))[0]
            c_x, c_y, c_m = tree.tree.nodes[cna_node]["genotype"]
            for n,g in self.genotypes.items():
                if g==(c_x, c_y):
                    end_node = n
                    break 
            #fix mutatated copies here! 
            if (root_id, end_node) in self.G.edges:
                self.G.remove_edge(root_id, end_node)
            self.G.add_edge(root_id, node_order[0])
            self.G.add_edge(node_order[-1], end_node)
        
        
        
       
            for child in self.node_dcfs:
                if child in node_order:
                    continue
                for parent in node_order:
                    if parent == child:
                        continue
                    if self.node_dcfs[parent] >= self.node_dcfs[child] and self.genotypes[parent]==self.genotypes[child]:
                        self.G.add_edge(parent, child)
                
                
            # for j,n in zip(sorted_clusters, node_order):
            #     self.G.
            #     self.add_mutation_edge(T_Seg, tree_clust[0], j, dcfs[j],n)
                





    
    # def build_subtrees(self, tree_clust, k, snv_clusters,dcfs):
    #     subtrees = []
    #     tree = self.id_to_tree[tree_clust[0]]
        
    #     sorted_clusters = np.argsort(dcfs)[::-1]
    #     node_id = np.array(list(snv_clusters.keys()))
    #     node_id.sort()
    #     # node_id=node_id[::-1]
    #     node_order = node_id[sorted_clusters]
    #     #if we only have one cluster or the mutation cluster on the path from the root to a 
    #     #to a CNA event, we only have 1 possible subtree
    #     if not tree.is_leaf(tree.split_node):
            
            
    #         T_Seg = deepcopy(self.T_Seg)
    #         for j,n in zip(sorted_clusters, node_order):
    #             self.add_mutation_edge(T_Seg, tree_clust[0], j, dcfs[j],n)
    #                         # self.print_verb(self)
                        
    #             # self.mut_mapping[n]= snv_clusters[n]
    #         subtrees.append(T_Seg)
        
    #     #if mutation clusters between the root and a CNA event the ordering must be linear
    #     else:
    #             x,y, m = tree.split_node_parent_geno
    #             root_id = self.T_CNA.node_mapping[(x,y)]
    #             init_node_to_dcfs = {n: self.node_dcfs[n] for n in node_order }
    #             subtrees = self.enumerate_subtrees(init_node_to_dcfs, root_id=root_id, geno=(x,y))
    #     print(f"length of subtrees: {len(subtrees)}")
    #     return subtrees
    
    def draw_graph(self, fname)-> None:
        G_pg = nx.nx_agraph.to_agraph(self.G)
        node_labels = {}
        for n in self.G:
            if n in self.node_dcfs:
                node_labels[n] = f"{n}\n{self.node_dcfs[n]}\n{self.genotypes[n]}"
            else:
                node_labels[n] =f"{n}\n{self.genotypes[n]}"

        for node, label in node_labels.items():
            G_pg.get_node(node).attr['label'] = label

        # Set graph attributes (optional)
        G_pg.graph_attr['label'] = 'Ancestral Graph'
        G_pg.graph_attr['rankdir'] = 'LR'  # Left to right layout

        G_pg.draw(fname, format='png', prog='dot')

    # Draw the graph to a file (optional)


    @staticmethod
    def find_root(tree)-> int:
        return [n for n in tree if tree.in_degree[n]==0][0]
    

    # def find_root_candidates(self,base_tree, ext_tree):
    #     ext_root= self.find_root(ext_tree)
    #     ext_geno = ext_tree.nodes[ext_root]["genotype"]
    #     ext_root_children= list(ext_tree.neighbors(ext_root))
    #     max_dcf = max([self.node_dcfs[c] for c in ext_root_children if c in self.node_dcfs])
    #     root_cand = []
    #     for n in base_tree:
    #         if n in self.node_dcfs:
    #             if base_tree.nodes[n]["genotype"] ==ext_geno and self.node_dcfs[n] >= max_dcf:
    #                 root_cand.append(n)

    #     return root_cand


    # def merge_and_rename(trees, root_property):
    #     if len(trees) == 1:
    #         return trees

    #     merged_trees = []
    #     current_tree = trees[0]
    #     remaining_trees = trees[1:]

    #     #get all possible root nodes 
    #     for t in 

    #     for root_label in root_property[len(trees) - 2]:
    #         # Rename the root node based on the provided root_property
    #         current_tree = nx.relabel_nodes(current_tree, {1: root_label})
    #         merged_subtrees = merge_and_rename(remaining_trees, root_property)
    #         merged_trees.extend(nx.compose_all([current_tree] + subtrees) for subtrees in merged_subtrees)

    #     return merged_trees

    # def get_all_compositions(self, combos):
    #     all_compositions = []
   
    #     tree_list= []
    #     for i,t1 in enumerate(combos):
    #         t1_roots = [t1]
    #         t1_root= self.find_root(t1)
    #         for j,t2 in enumerate(combos):
    #             if j  != i:
    #                 cand_roots = self.find_root_candidates(t2, t1)
    #                 for r in cand_roots:
    #                         t1_relab = nx.relabel_nodes(t1,{t1_root: r})
    #                         t1_roots.append(t1_relab)
    #         tree_list.append(t1_roots)
    #     tree_list.append([self.T_Seg.copy()])
    #     combinations = list(itertools.product(*tree_list))
    #     #this combines all trees by matching the copy number states
    #     for combo in combinations:
    #         new_trees =  nx.compose_all(combo)
    #         all_compositions.append(new_trees)
    #     return all_compositions

        

            






        

                  
             



            # for c in combos:
                   
                  
                
            
            #         subtree = deepcopy(c)
                    
            #         for r in join_cand:
            #             subtree_root = self.find_root(subtree)
            #             subtree =subtree.relabel_nodes({subtree_root: r})
            #             T_Seg = nx.compose(subtree, T_Seg)
    # def combine_subtrees(self, subtrees_list):
        
        


    #     combinations = list(itertools.product(*subtrees_list))

    #     #this combines all trees by matching the copy number states
    #     for combo in combinations:
    #         new_trees =  self.get_all_compositions(combo)

        
         
    #         self.T_Seg_List.extend(new_trees)
        # self.combine_non_roots()
        
        #however we might also want to join on mutation nodes
        # for combo in combinations:
        #     T_Seg = deepcopy(self.T_Seg)
        #     for c in combo:
                #find all nodes in T_Seg that could be a root for c 
                #relabel c with the node id of that mutation node then compose



            
                

        

        # new_list = []
        # for T_Seg in T_Seg_list:
          
        #     next_node =max(T_Seg_base.nodes) + 1
        #     for i, subtrees in enumerate(subtrees_list):

              
        #         for tree in subtrees:
        #             T_Seg = nx.compose(T_Seg_base, tree)
        #             new_list.append(T_Seg)
                
                

                    
                    # tree_root_node = [n for n in tree if tree.in_degree[n]==0][0]
                    # tree_root_geno =tree.nodes[tree_root_node]["genotype"]  #(x,y)
                    # root_ids = [n for n in T_Seg_base if T_Seg_base.nodes[n]["genotype"]==tree_root_geno ]
                    
                    # for r in root_ids:
                    #     node_mapping = {-1: r}             
              
                    #     node_to_dcfs= {}
                   
                    #     for i,dcf in enumerate(dcfs):
                    #         node_id = next_node + i
                    #         node_mapping[i] = node_id
                    #         self.mut_mapping[node_id] = snv_clusters[i]
                    #         node_to_dcfs[node_id] = dcf
                    #     for subtree in subtrees:
                    #         subtree = nx.relabel_nodes(subtree, node_mapping)
                        
                    #         new_list.append(T_Seg)
                    # self.T_Seg_List= new_list     

                                    #identify compatible nodes for joining

              
                # for n in nx.dfs_postorder_nodes(T_Seg_base, self.root):

       


            
    

                

      
    
    def fit(self, data, segment):

        self.data = data
        self.cells_by_cn = self.data.cells_by_cn(segment)
        
        #not sure I need this, but let's keep it for now
        snvs, alt, total = self.data.count_marginals(segment)
        self.snvs = snvs 
        snv_index = {s: i for i,s in enumerate(snvs)}
             
        self.m = len(self.snvs)
        
        cell_counts_by_snv, cells_by_snvs = self.data.count_cells_by_snv(segment)

        
        self.print_verb(f"Average cell count per SNV: {cell_counts_by_snv.mean()}")
        #construct overlap graph and check if it is fully connected
        # G = self.construct_overlap_graph(snvs, cells_by_snvs)
        # num_components = nx.number_connected_components(G)
   
        # self.print_verb(f"Overlap graph for segment {segment} has {num_components} component(s)")


        #initialize T_SNV assignment 
        self.T_SNV_dict = {}

        for s in self.snvs:
            max_like= np.NINF
            # if s in [1703, 3158, 3526]:
            #     print("check")
            cells_by_cn = {cn: np.intersect1d(cells_by_snvs[s], self.cells_by_cn[cn]) for cn in self.cells_by_cn}
      
            for tree in self.T_SNVs:
                # print(tree)
                like, cell_assign =tree.likelihood(cells_by_cn, data.var[:,s], data.total[:,s],0.001)
                if like > max_like:
                    max_like = like
                    self.T_SNV_dict[s] = tree.id
        self.mutated_copies = {}
        for s in self.snvs:
            tree = self.id_to_tree[self.T_SNV_dict[s]]
            self.mutated_copies[s]=tree.m_star
        # T_SNV_Clusters = [[0], [1,2], [3]]
        starting_cluster = max(self.T_CNA.tree.nodes) + 1
        subtrees_list = []
        for tree_clust in self.T_SNV_Clusters:
            clust_snvs = [s for s in snvs if self.T_SNV_dict[s] in  tree_clust]
            if len(clust_snvs) ==0:
                continue
       
            k, snv_clusters, dcfs = self.cluster_snvs(clust_snvs, snv_index, alt, total)
    
            self.print_verb(f"{tree_clust}: k={k}, dcfs={dcfs}")
            # dcfs = dcfs.reshape(-1)
            # sorted_clusters = np.argsort(dcfs)[::-1]
            # sorted_dcfs = np.sort(dcfs)[::-1]
            # for c in np.unique(snv_clusters):
            #     self.print_verb(f"{c}: {snv_clusters[snv_clusters==c].shape[0]}")
            #sort the 
            tree = self.id_to_tree[tree_clust[0]] 
            # print(tree)
            
            #this ensures that the node and cluster id always match up
            snv_clusters = {j+starting_cluster: clust for j, clust in snv_clusters.items()}
            for j in range(k):
                self.node_dcfs[j+ starting_cluster] = dcfs[j]
                self.mut_mapping[j + starting_cluster]  = snv_clusters[j+starting_cluster]
            # subtrees =self.build_subtrees(tree_clust, k, snv_clusters,dcfs) 
            self.add_to_graph(tree_clust, k, snv_clusters, dcfs) 
            # subtrees_list.append(subtrees)
            starting_cluster += k
        self.draw_graph("test/ancestral_graph.png")
        self.T_Seg_List = self.enumerate_segment_trees()
        print(f"Segment {segment}: Found {len(self.T_Seg_List)} candidate segment trees, using likelihood to prioritize....")
        # self.combine_subtrees(subtrees_list)
             

            # if k==1 or not tree.is_leaf(tree.split_node):
            #     for j in sorted_clusters:
            #         #modify in place
            #         for T_Seg in self.T_Seg_List:
                
            #             node_id = self.add_mutation_edge(T_Seg, tree_clust[0], j, dcfs[j])
            #             # self.print_verb(self)
                    
            #             self.mut_mapping[node_id]= snv_clusters[j]
            # else:
    
               


                    #

                    

                #we need to enumerate all possible subtrees
                



          
                        

        self.cell_mapping = {}
        self.total_cn_by_sample = {}
        for s in self.cells_by_cn:
        
            self.total_cn_by_sample[s] = self.data.copy_numbers[self.cells_by_cn[s], :][0,0]
    

  
        SegTree = self.optimal_clonal_tree(data, segment)

        return SegTree
                       

            
            

            

               
            
            


            


            
            #initially cluster snvs 
     

  
    
      
        # self.total_cn_by_sample = {}
        

  
        
        # tree_id_to_indices= {}
        # for i,s in enumerate(self.snvs):
        #     id = self.T_SNV_dict[s]
        #     if id in tree_id_to_indices:
        #         tree_id_to_indices[id].append(i)
        #     else:
        #         tree_id_to_indices[id] = [i]

        # for k in np.unique(self.clusters):
            
        #     state_tree_snvs_pairs = self.get_unique_state_trees(self.snvs_by_cluster[k])
        #     for id, _ in state_tree_snvs_pairs:
        #         print(self.id_to_tree[id])
        #     if len(state_tree_snvs_pairs) > 0:
        #         id, tr = state_tree_snvs_pairs[0]
        #         print(tr)
          
        #     print(self)
                #   self.mut_mapping[node_id] = muts
       



    @staticmethod
    def enumerate_subtrees(node_to_dcfs, root_id, geno):
   
        #construct the ancestry graph
        G = nx.DiGraph()
        G.add_node(root_id)

        for i,d_par in node_to_dcfs.items():
            G.add_node(i)
            G.add_edge(root_id, i)
            for j,d_child in node_to_dcfs.items():
                if i != j and d_par >= d_child:
                    G.add_edge(i, j)
        
        iterator = nx.algorithms.tree.branchings.ArborescenceIterator(G)
        subtrees = []
        for arbor in iterator:
            nx.set_node_attributes(arbor, geno, name="genotype")
            nx.set_edge_attributes(arbor, "mutation", name="event")
            subtrees.append(arbor)

        return subtrees
    

    def enumerate_segment_trees(self):
   
        iterator = nx.algorithms.tree.branchings.ArborescenceIterator(self.G)
        subtrees = []
        for arbor in iterator:
            nx.set_node_attributes(arbor, self.genotypes, 'genotype')
            # nx.set_node_attributes(arbor, geno, name="genotype")
            # nx.set_edge_attributes(arbor, "mutation", name="event")
            subtrees.append(arbor)

        return subtrees


            

        


        
    

     
    # @time_it
    def map_assign_cells(self, T_Seg, data, seg):
        # print(self.tree_to_string(T_Seg))
        clone_order = list(nx.dfs_preorder_nodes(T_Seg, source=self.root))
        # clone_order = [c for c in clone_order if c != self.root]

        total_cn_states = {}
        like_list =[]
        seg_snvs = data.seg_to_snvs[seg]
        cell_mapping = {}
        cell_mapping[self.root] = []
        for n in clone_order:
            x,y = T_Seg.nodes[n]["genotype"]
            total_cn = x+y 
            total_cn_states[n]= total_cn
            cna_geno = np.full(self.m, total_cn, dtype=int).reshape(1,-1)
            clone = nx.shortest_path(T_Seg, source=self.root, target=n)
            # if T_Seg[clone[-2]][clone[-1]]["event"]== "mutation":
            snvs = []
            for c in clone:
                if c in self.mut_mapping:
                    snvs += self.mut_mapping[c]
                
            y_vec = np.zeros(shape=data.M, dtype=int)
           
            y_vec[snvs] = 1
            y_vec= y_vec[seg_snvs]
            y_vec = y_vec.reshape(1,-1)

            # cell_by_snv_like = np.zeros((data.var.shape[0],self.m))
            # for i in range(self.data.var.shape[0]):
            #     for j in range(self.m):
            #         a = self.data.var[i,j]
            #         d = self.data.total[i,j]
            #         c = cna_geno[:,j][0]
            #         alpha= self.alpha 
            #         y = y_vec[:,j][0]
                    
            #         out = likelihood_function(a,d,y,c,alpha)
            #         if np.log(out) > 0:
            #             print(f"a:{a} d:{d} y:{y} c:{c}: prob: {out} logprob:{np.log(out)}")
            #         cell_by_snv_like[i,j] =out
            # print(cell_by_snv_like)
            cell_like = likelihood_function(data.var[:, seg_snvs], data.total[:,seg_snvs], y_vec, cna_geno, 0.001)
         
            # assert(np.array_equal(cell_by_snv_like, cell_by_snv_like2))
            # print(cell_by_snv_like2)
            # cell_by_snv_like = np.log(cell_by_snv_like)
            # print(cell_by_snv_like)
            # print(f"node: {n} :{cell_by_snv_like.sum(axis=1)}")
            like_list.append(cell_like)
     

        #rows are nodes and columns are cells
        cell_likes = np.vstack(like_list)

        self.likelihood =0
        cell_assign = {}
 
        #now only consider clone assignments for cells with the correct CN
        for s in self.cells_by_cn:
            # cn = self.total_cn_by_sample[s]
            clone_filter = [total_cn_states[n]==s for n in clone_order]

            clones = [c for i,c in enumerate(clone_order) if clone_filter[i]]
            for c in clones:
                cell_mapping[c] = []
            cells = self.cells_by_cn[s]
            sample_likes = cell_likes[clone_filter][:, cells]

            like = np.max(sample_likes, axis=0).sum()
            self.likelihood += like
            map_clone = np.argmax(sample_likes, axis=0)
            for c, m in zip(cells, map_clone):
                cell_assign[c]= clones[m]
   
                cell_mapping[clones[m]].append(c)
          
            

      
      
        self.print_verb(f"Log Likelihood: {self.likelihood}")
        if self.verbose:
            for n in T_Seg:

                print(f"Node{n}: {len(cell_mapping[n])} cells assigned")

        # cell_mapping = {n: np.array(cell_mapping[n], dtype=int) for n in T_Seg}
        # self.print_verb(self)
        return  cell_mapping, self.likelihood
 
            






        




    def add_mutation_edge(self, T_Seg, id, cluster, edge_dcf, n):
        tree = self.id_to_tree[id]

                    
 
        # next_node  = max([n for n  in T_Seg]) + 1
        next_node = n
        # edge_dcf = tuple([d for d in self.DCF[:,cluster]])
        m_0, m_m = tree.find_split_pairs()
        

        u, u_geno = m_0
        v, v_geno = m_m

       

     
        leaves = [n for n in T_Seg if T_Seg.out_degree[n]==0]
        for l in leaves:
            root_to_leaf_path = nx.shortest_path(T_Seg, source=self.root, target=l)
            #does  the mutation edge occur within or after the path
            start = root_to_leaf_path[0]
            end = root_to_leaf_path[-1]
            start_geno =   T_Seg.nodes[start]["genotype"]
            end_geno = T_Seg.nodes[end]["genotype"]
            node_path, snv_genotypes = tree.find_path(start_geno, end_geno)

            #insert 
            if v_geno in snv_genotypes:
                T_Seg.add_node(next_node, genotype=(v_geno[0], v_geno[1] ))
      

                if v_geno != snv_genotypes[-1]:
                    for i,s in enumerate(snv_genotypes):
                        if s ==v_geno:
                            p_x, p_y, p_m= snv_genotypes[i-1]
                            c_x, c_y, c_m = snv_genotypes[i+1]
                            break
            
                    for i,parent in enumerate(root_to_leaf_path):
                        if i >= len(root_to_leaf_path) -1:
                            break
                        cp_x, cp_y=  T_Seg.nodes[parent]["genotype"]
                        child_node = root_to_leaf_path[i+1]
                        cc_x, cc_y =T_Seg.nodes[child_node]["genotype"]

                        if p_x == cp_x and p_y==cp_y and cc_x == c_x and cc_y==c_y:
                            etype = T_Seg[parent][child_node]["event"]
                    
                            T_Seg.remove_edge(parent, child_node)
                            T_Seg.add_edge(parent, next_node, event="mutation", cluster=cluster, dcf =edge_dcf)
                  
                            T_Seg.add_edge(next_node, child_node, event=etype)
                            return next_node
                else:
                    parent = root_to_leaf_path[-1]
                    T_Seg.add_edge(parent, next_node,event="mutation", dcf =edge_dcf, cluster=cluster)
                    return next_node 
     
        else:
            T_Seg.add_node(next_node, genotype=(v_geno[0], v_geno[1] ))
            T_Seg.add_edge(self.root, next_node,event="mutation", dcf =edge_dcf, cluster=cluster)
            return next_node

        #if it wasn't added anywhere add as a child of the root

            


  




        # if len(T_Seg) ==0:
                
        #         self.T_Seg = deepcopy(tree.tree)
             
        #         parent = list(self.T_Seg.predecessors(node_id))[0]
              
        #         self.T_Seg[parent][node_id]["DCF"] = edge_dcf
             
        # else:
        #     preorder_nodes = nx.dfs_preorder_nodes(self.T_Seg)
        
        #     for n in preorder_nodes:
        #         x,y,m = self.T_Seg.nodes[n]["genotype"]
        #         # if m > 0:

        #         # print(tree.tree.nodes[n]['genotype'])

        # return node_id 
                    


    def get_unique_state_trees(self, snvs):
        '''
        Takes in a set of SNVs and returns a tuple pairs of state tree and snv sets
        '''
        trees = [self.T_SNVs[s].id for s in snvs]

        unique_trees = []
        pairs = []
        for i,t in enumerate(trees):
            if t not in unique_trees:
                unique_trees.append(t)
                muts = self.tree_to_snvs[t]     
                pairs.append((t, np.intersect1d(muts, snvs)))
        return pairs

    


    # def add_mutation_edge(self, id, snvs, cluster):
    #     mut_tree = self.id_to_tree[id]
    #     DCF = self.DCF[:,cluster]
    #     preorder = nx.dfs_preorder_nodes(self.T_Seg)
    #     split_cand = []
    #     #find split copy number edge
    #     geno= self.find_split_edge(mut_tree)
    #     start_node = None
    #     for n, geno in preorder:
    #         x,y,m = n["genotype"]
    #         if (x,y) == (g1[0], g1[0]) and start_node is None:
    #             start_node =n 
    #         else if (x,y) =(g1[0], g1[0]):
    #             end_node = n 
    #      path  = nx.shortest_path(self.T_Seg, start_node, end_node)
    #      #check if a mutation edge already exists in the path
    #      for p in path:
    #         x,y,m = p["genotype"]
    #         if m ==1:
    #             p["DCF"]        
         
        
            
    class TreeMerge:
        def __init__(self, jaccard_threshold=0.05) -> None:
            self.jaccard_threshold = jaccard_threshold 
 

        def merge(self, tree1, tree2, data):
            self.T1 = tree1
            self.T2 = tree2
            self.data = data 
            self.trees = [self.T1, self.T2]
            self.snvs1 = self.T1.get_all_muts()
            self.snvs2 = self.T2.get_all_muts()

            connected_components = self.constuct_bipartite_graph()
   
            '''
            identify any singleton components
            '''

            


            merge_graph = self.construct_merge_graph(connected_components)

            clonal_trees = self.generate_candidate_trees(merge_graph)

            opt_clonal_tree = self.find_optimal(clonal_trees)

            return opt_clonal_tree
        
        def generate_candidate_trees(self, G):
             
            cand_trees = []
            iterator = nx.algorithms.tree.branchings.ArborescenceIterator(G)
            for arbor in iterator:
         
                mut_mapping = {}
                tree= nx.DiGraph()
                node_mapping = {}
                index = 0
                #traverse the tree in preorder starting from the root
                # root = -1
                # index = 0 
                # root_neighbors  = tree.neighbors(-1)
                # if len(root_neighbors) ==1:
                #     root = root_neighbors[0]

                # node_mapping[root] = index 
                # index +=1 
                # tree.add_node(index )
                # mut_mapping[index] = 

                #we need to first add all clustered nodes
                clust_nodes = []
                for u,v, attr in arbor.edges(data=True):
                
                    
                    if attr["etype"] == CLUST:
                        '''
                        We only neeed to add one node to the that has 
                        '''
                        snvs = []
                        for gnode in [u,v]:
                            clust_nodes.append(gnode)
                            if gnode not in node_mapping:
                                node_mapping[gnode] =index 
                            t, n = gnode
                            snvs += self.tree[t].mut_mapping[n]
                        node_id = node_mapping[u]
                        tree.add_node(node_id)
                        mut_mapping[node_id] = snvs 
                for u,v, attr in arbor.edges(data=True):
                    if attr["etype"] == ANC:
                        if u not in node_mapping:
                                node_mapping[u] = index 
                                index+=1
                        if v not in clust_nodes:  #incoming mutations need to be added to the merged tree
                            node_mapping[v] = index 
                            index += 1
                            t, n = v
                            mut_mapping[node_mapping[v]]= self.tree[t].mut_mapping[n]
                      
                        tree.add_edge(node_mapping[u], node_mapping[v])
                
                cand_trees.append(ClonalTree(0, tree, mut_mapping = mut_mapping))
            return cand_trees 

                    


        def find_optimal(self, clonal_trees):
            likelihoods = [ct.compute_likelihood(self.data) for ct in clonal_trees]
            opt_index = np.argmax(likelihoods)
            return clonal_trees[opt_index]
        
        def ancestral_weight(self, t_parent, parent, t_child, child):

            #what are the snvs that are introduced in the parent node
            parent_snvs = t_parent.get_muts(parent)
            if len(parent_snvs) ==0:
                return np.Inf
            
            desc_cells = t_child.get_desc_cells(child)



            #now compute the observed VAF for each parent snvs across pooled desc cells
            obs_vafs = self.data.compute_vaf(parent_snvs, desc_cells)
            exp_vaf = []
            for i,s in enumerate(parent_snvs):
                g = self.data.seg_by_snv[s]
                mut_copies = [t_parent.get_mut_copies(g,s, c) for c in desc_cells]
                mut_copy = mode(mut_copies)
                tot_copies = self.data.copy_profiles[desc_cells,g]
                exp_vaf.append(mut_copy/tot_copies)
            
            return np.mean(np.abs(exp_vaf- obs_vafs))
        

        def clustered_weight(self, t0, n0, t1, n1):
            clusters = [(t0, n0), (t1, n1)]
            #what are the snvs that are introduced in the parent node
            snvs = []
            cells = []
            for i,c in enumerate(clusters):
                t, n =c 
                snvs[i] =t.get_muts(n)
                cells[i] = t.get_desc_cells(n)
            
         
            



            if len(snvs[0] + snvs[1]) ==0:
                return np.Inf
            
            
            # cell_union = list(set(cells[0]).union(set(cells[1])))


            obs_vafs, exp_vafs  = []
         
            #now compute the observed VAF for each parent snvs across pooled desc cells
            rev_cells = cells[::-1]
            for i,node_snvs in enumerate(snvs):
                tree = clusters[i][0]
                evaf = []
                for dcells in rev_cells:


                    obs_vafs.append(self.data.compute_vaf(node_snvs, dcells))
              
                    for i,s in enumerate(node_snvs):
                        g = self.data.seg_by_snv[s]
                        #TODO need to fix segment tree to properly store mut_copies in seg tree
                        # mut_copies = [tree.get_mut_copies(g, s, c) for c in dcells]
                        #for now, use 1 mutated copy but this is just a placehold for debugging
                        mut_copies = [1 for c in dcells]
                        mut_copy = mode(mut_copies)
                        tot_copies = self.data.copy_profiles[dcells,g]
                        evaf.append(mut_copy/tot_copies)
                exp_vafs.append(evaf)
            obs_vafs = np.concatenate(obs_vafs)
            exp_vafs = np.concatenate(exp_vafs)
            
            return np.mean(np.abs(exp_vafs- obs_vafs))


            
        def compute_edge_weight(self, parent, child):
            t_p, n_p = parent
            t_c, n_c = child 

            
            snvs = self.trees[t_p].get_muts(n_p)
            if len(snvs) ==0:
                return 0, None

            #if parent is ancestral to child
            anc_weight = self.ancestral_weight(self.trees[t_p], n_p, self.trees[t_c], n_c)

            #if parent and child snvs should be clustered togeter
            clust_weight = self.clustered_weight(self.trees[t_p], n_p, self.trees[t_c], n_c)

            if anc_weight < clust_weight:
                return anc_weight, ANC
            else:
                return clust_weight, CLUST
            



        def construct_merge_graph(self, connected_components):
            merge_graph = nx.DiGraph()
        
            #
            
            #first figure out what to do with the root
     
            merge_root =  np.all([len(self.trees[t].get_muts(self.trees[t].root)) ==0 for t in [0,1] ])
            merge_graph.add_node(-1)
            if not merge_root:
              
                for t in [0,1]:
                    rt = self.trees[0].root
                    self.add_edge(-1, (t,rt), weight=0, etype=ANC)
            else:
                r0 = self.trees[0].root
                r1 = self.trees[1].root
                self.add_edge(-1, (0,r0), weight=0, etype=ANC)
                self.add_edge((0,r0), (1, r1), weight=0, etype=ANC)
            

         


            #now add all ancestral edges within a tree
            for t in [0,1]:
                tree = self.trees[t]
                for u,v in tree.tree.edges:
                    child = (t,v)
                    if u == tree.root and merge_root:
                        merge_graph.add_edge((1,r1), child, weight=0, etype=ANC)
                    else:
                        par = (t,u)
                        w, _ = self.compute_edge_weight((par, child))
                        merge_graph.add_edge(par, child, weight=w, etype=ANC)

            #now add cross edges 
            for comp in connected_components:
                merge_graph.add_nodes_from(comp)
                t_nodes = { tree: [n for t,n in comp if t==tree] for tree in [0,1]}
        

                for n0  in t_nodes[0]:
                    if n0 == self.trees[0].root and merge_root:
                        node = (1,r1)
                    else:
                        node = (0,n0)
                
                 
                    for n1 in t_nodes[1]:
                        if n1 == self.trees[1].root  and n0 == self.trees[0].root and merge_root:
                            continue
                  
        
                        w01, type01 = self.compute_edge_weight(node, (1,n1))
                        w10, type10 = self.compute_edge_weight((1,n1),node)
                        if type01 is not None:
                            merge_graph.add_edge((0,n0), (1,n1), weight=w01, etype=type01)
                        if type10 is not None:
                            merge_graph.add_edge((1,n1), (0,n0), weight=w10, etype=type10)
            return merge_graph



                


                        
                           




                      

                #group nodes into tree 1 and tree2"
                

                   

               





        @staticmethod
        def jaccard_similarity(set1, set2):
            # Calculate the size of the intersection
            intersection_size = len(set1.intersection(set2))
            
            # Calculate the size of the union
            union_size = len(set1.union(set2))
            
            # Calculate the Jaccard similarity
            jaccard_similarity = intersection_size / union_size
            return jaccard_similarity

                        
                    
    
        def constuct_bipartite_graph(self):
            self.B = nx.Graph()
            T_nodes = [[(i,n) for n in tree.cell_mapping] for i,tree in enumerate([self.T1, self.T2])]
            for i, nodes in enumerate(T_nodes):
                self.B.add_nodes_from(nodes, bipartite=i)  # Nodes in set1 have bipartite=0

            for n0 in T_nodes[0]:
                cells1 = self.T1.cell_mapping[n0]
                for n1 in T_nodes[1]:
                    cells2 = self.T2.cell_mapping[n1]
                    
                    if len(cells1) >0 or len(cells2) >0:
                            jacc = self.jaccard_similarity(set(cells1), set(cells2))
                            if jacc > self.jaccard_threshold:
                                self.B.add_edge((0,n0), (1,n1), weight= jacc)
            
            singletons = {t: [] for t in [0,1]}
            connected_components = list(nx.connected_components(self.B))
            for comp in connected_components:
                if len(comp) ==1:
                    t,n = comp[0]
                    singletons[t].append(n)



            '''
            preprocess the connected components, so that there are no singleton nodes
            are any singletons descendents of subtree_parent?
            if so, we are going to add them to be a member of that connected component
            '''
            for comp in connected_components:
                
                t_nodes = { tree: [n for t,n in comp if t==tree] for tree in [0,1]}
                for t, nodes in t_nodes.items():
                    tree = self.trees[t].tree
                    subtree_parents =self.find_closest_nodes_to_root(tree, self.trees[t].root, nodes)


                    singles = singletons[t]
                    for par in subtree_parents:
                        for s in singles:
                            if nx.has_path(tree, par, s):
                                comp.append((t,s))
            return connected_components

        @staticmethod    
        def find_closest_nodes_to_root(tree, root, nodes_list):
            # Perform BFS and get the shortest path lengths from the root
            shortest_paths = nx.single_source_shortest_path_length(tree, root)

            # Find the minimum distance from the root to the nodes in the input list
            min_distance = min(shortest_paths[node] for node in nodes_list)

            # Filter the nodes in the input list that are at the minimum distance from the root
            closest_nodes = {node for node in nodes_list if shortest_paths[node] == min_distance}

            return closest_nodes


            



