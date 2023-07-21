import networkx as nx
import numpy as np
import itertools
from copy import deepcopy
from scipy.stats import binom
from clonal_tree import ClonalTree
from sklearn import metrics
from sklearn.cluster import KMeans
import time

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
        
        if d ==0:
            return 1e-10
        elif y ==0:
            val =  binom.pmf(a,d,alpha)
            return val
          
        else:
            vaf = np.arange(1, c)/c
            # vaf  = 1/c
            # vaf = 0.5
            adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
            val = (1/vaf.shape[0])*np.sum(binom.pmf(a,d,adjusted_vaf))
            # return binom.pmf(a,d,adjusted_vaf)
            return val
like_func = np.vectorize(likelihood_function)



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



class BuildSegmentTree:
    def __init__(self, T_CNA, seed=1026, max_clusters=4, silhouette_min=0.6, silhouette_improve=0.2, verbose=False):
        self.verbose =verbose 
        self.T_CNA = T_CNA
        self.max_clusters = max_clusters
        self.sil_min = silhouette_min
        self.sil_improve = silhouette_improve
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

        self.T_Seg_List = [T_Seg]

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
        dcfs = dcfs.reshape(-1,1)
        clust_dcfs = np.mean(dcfs,axis=0)
  

        snv_clusters =  np.full_like(clust_snvs, 0)
        best_k = 1
        if len(np.unique(dcfs)) > 1:
        
        
            max_score = -1
           
            for k in range(2, self.max_clusters+1):
            
                km = KMeans(k, random_state=self.rng)
                clusters = km.fit_predict(dcfs)
                # for c in np.unique(clusters):
                #     print(f"{c}: {clusters[clusters==c].shape[0]}")
                cluster_centers = km.cluster_centers_  #nclusters x nfeatures
                sil_score = metrics.silhouette_score(dcfs, clusters)
                if max_score > 0:
                    sil_improve = 1-((1-sil_score)/(1-max_score))
                else:
                    sil_improve = 1
                
                # print(f"{k}: {sil_score}")
                if sil_score > self.sil_min and sil_score > max_score and sil_improve > self.sil_improve:
                    max_score = sil_score
                    best_k = k
                    clust_dcfs =cluster_centers
                    snv_clusters = clusters
        cluster_map = {j: [] for j in range(best_k)}
        for s,j in zip(clust_snvs, snv_clusters):
            cluster_map[j].append(s)

        

        return best_k, cluster_map, clust_dcfs
    
    def optimal_clonal_tree(self, data, g):
        best_like = np.NINF
        # print(self)
        for T_Seg in self.T_Seg_List:
            # print(self.tree_to_string(T_Seg))
            segtree = ClonalTree(g, T_Seg, self.mut_mapping, mutated_copies= self.mutated_copies)
            cell_mapping, loglike= self.map_assign_cells(segtree.tree, data,g)
            segtree.set_cell_mapping(cell_mapping) 
            # loglike = segtree.compute_likelihood(data, g)
            if loglike > best_like:
                best_like = loglike
                opt_segtree = segtree 
      
        opt_segtree.loglikelihood = best_like
        return opt_segtree
    
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
            cells_by_cn = {cn: np.intersect1d(cells_by_snvs[s], self.cells_by_cn[cn]) for cn in self.cells_by_cn}
      
            for tree in self.T_SNVs:
                like, cell_assign =tree.likelihood(cells_by_cn, data.var[:,s], data.total[:,s],0.001)
                if like > max_like:
                    max_like = like
                    self.T_SNV_dict[s] = tree.id
        self.mutated_copies = {}
        for s in self.snvs:
            tree = self.id_to_tree[self.T_SNV_dict[s]]
            self.mutated_copies[s]=tree.m_star
        # T_SNV_Clusters = [[0], [1,2], [3]]
        for tree_clust in self.T_SNV_Clusters:
            clust_snvs = [s for s in snvs if self.T_SNV_dict[s] in  tree_clust]
            if len(clust_snvs) ==0:
                continue
       
            k, snv_clusters, dcfs = self.cluster_snvs(clust_snvs, snv_index, alt, total)
    
            self.print_verb(f"{tree_clust}: k={k}, dcfs={dcfs}")
            dcfs = dcfs.reshape(-1)
            sorted_clusters = np.argsort(dcfs)[::-1]
            # sorted_dcfs = np.sort(dcfs)[::-1]
            # for c in np.unique(snv_clusters):
            #     self.print_verb(f"{c}: {snv_clusters[snv_clusters==c].shape[0]}")
            #sort the 
            tree = self.id_to_tree[tree_clust[0]]    
            if k==1 or not tree.is_leaf(tree.split_node):
                for j in sorted_clusters:
                    #modify in place
                    for T_Seg in self.T_Seg_List:
                
                        node_id = self.add_mutation_edge(T_Seg, tree_clust[0], j, dcfs[j])
                        # self.print_verb(self)
                    
                        self.mut_mapping[node_id]= snv_clusters[j]
            else:
    
                x,y, m = tree.split_node_parent_geno
                init_node_to_dcfs = {i: dcf for i,dcf in enumerate(dcfs) }
                subtrees = self.enumerate_subtrees(init_node_to_dcfs, root_id=-1, geno=(x,y))
                new_list = []
                for T_Seg_base in self.T_Seg_List:
                   
                    next_node =max(T_Seg_base.nodes) + 1
                    # for n in nx.dfs_postorder_nodes(T_Seg_base, self.root):

                    root_ids = [n for n in T_Seg_base if (T_Seg_base.out_degree[n]==0 or n==self.root) and T_Seg_base.nodes[n]["genotype"]==(x,y) ]
                    root_id = root_ids[0]
                    node_to_dcfs= {}
                    node_mapping = {-1: root_id}
                    for i,dcf in enumerate(dcfs):
                        node_id = next_node + i
                        node_mapping[i] = node_id
                        self.mut_mapping[node_id] = snv_clusters[i]
                        node_to_dcfs[node_id] = dcf
                    for subtree in subtrees:
                        subtree = nx.relabel_nodes(subtree, node_mapping)
                        T_Seg = nx.compose(T_Seg_base, subtree)
                        new_list.append(T_Seg)
                self.T_Seg_List= new_list


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
    



            

        


        
    

     
    @time_it
    def map_assign_cells(self, T_Seg, data, seg):
        # print(self.tree_to_string(T_Seg))
        clone_order = list(nx.dfs_preorder_nodes(T_Seg, source=self.root))
        clone_order = [c for c in clone_order if c != self.root]

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
            cell_by_snv_like = like_func(data.var[:, seg_snvs], data.total[:,seg_snvs], y_vec, cna_geno, 0.001)
         
            # assert(np.array_equal(cell_by_snv_like, cell_by_snv_like2))
            # print(cell_by_snv_like2)
            cell_by_snv_like = np.log(cell_by_snv_like)
            # print(cell_by_snv_like)
            # print(f"node: {n} :{cell_by_snv_like.sum(axis=1)}")
            like_list.append(cell_by_snv_like.sum(axis=1))
     

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
 
            






        




    def add_mutation_edge(self, T_Seg, id, cluster, edge_dcf):
        tree = self.id_to_tree[id]

                    
 
        next_node  = max([n for n  in T_Seg]) + 1
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
         
        
            
        
                        
                        
                    
    





            



