import networkx as nx
import numpy as np
import itertools
from copy import deepcopy
from scipy.stats import binom



def likelihood_function( a,d,y,c, alpha):
        
        if d ==0:
            return 0
        elif y ==0:
        
            return binom.pmf(a,d,alpha)
          
        else:
            vaf = np.arange(1, c)/c
            # vaf  = 1/c
            # vaf = 0.5
            adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
            # return binom.pmf(a,d,adjusted_vaf)
            return (1/vaf.shape[0])*np.sum(binom.pmf(a,d,adjusted_vaf))
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
    def __init__(self, id, T_CNA, T_SNVs, DCF,clusters,cells_by_sample, data=None):
        self.id = id
        self.DCF = DCF
        self.T_CNA = T_CNA
        self.T_SNVs = T_SNVs
        self.tree_to_snvs = {}
        self.clusters = clusters 
        self.id_to_tree = {tree.id: tree for key, tree in self.T_SNVs.items()}
        for s, tree in self.T_SNVs.items():
            if tree.id not in self.tree_to_snvs:
                self.tree_to_snvs[tree.id] = [s]
            else:
                self.tree_to_snvs[tree.id].append(s)
        
        self.snvs = list(T_SNVs.keys())
        self.m = len(self.snvs)
        self.data = data
        self.T_Seg = nx.DiGraph()
        self.root = self.T_CNA.root
        for u,v in self.T_CNA.tree.edges:
            u_x, u_y, _ = self.T_CNA.tree.nodes[u]["genotype"]
            v_x, v_y, _ = self.T_CNA.tree.nodes[v]["genotype"]
            self.T_Seg.add_node(u, genotype= (u_x, u_y))
            self.T_Seg.add_node(v, genotype= (v_x, v_y))
            if u_x + u_y < v_x + v_y:
                etype = "gain"
            elif u_x + u_y >  v_x + v_y:
                etype = "loss"

            self.T_Seg.add_edge(u,v ,event=etype)
        
        self.root = [n for n in self.T_Seg if self.T_Seg.in_degree[n]==0][0]

        self.snvs_by_cluster = {k: np.where(self.clusters == k)[0] for k in np.unique(self.clusters)}

        self.likelihood= 0
        self.cluster_ids = np.unique(self.clusters)
        self.mut_mapping = {}
        self.cell_mapping = None
        self.mut_loss_mapping = {}
        if self.data is not None:
            self.cells = self.data.cell_lookup.index.to_list()
        else:
            self.cells = []
        self.alpha = 0.001
        self.total_cn_by_sample = {}
        self.cells_by_sample = cells_by_sample
        for s in self.cells_by_sample:
            
            self.total_cn_by_sample[s] = self.data.copy_numbers[self.cells_by_sample[s], :][0,0]
        






    
    def __str__(self) -> str:
        mystr =""

        for u in nx.dfs_preorder_nodes(self.T_Seg, source=self.root):
            for v in self.T_Seg.neighbors(u):
    
                etype = self.T_Seg[u][v]["event"]
                if etype == "mutation":
                    dcf = self.T_Seg[u][v]["dcf"]
                    cluster = self.T_Seg[u][v]["cluster"]
                else:
                    dcf = "NA"
                    cluster = "NA"
                mystr+= f"Node {u}: {self.T_Seg.nodes[u]['genotype']} -> Node {v}: {self.T_Seg.nodes[v]['genotype']}:"
                mystr+= f"{etype}: cluster: {cluster} dcf: {dcf}\n"
        return mystr

    def order_mutation_edges(self):
        # pass 
        #collapse edges that are all in the same cluster
        #reassign clusters to be in accordance with DCF values
        leaves = [n for n in self.T_Seg if self.T_Seg.out_degree[n]==0]
   
        for l in leaves:
            dcf_vals = {}
            mut_edges = []
            root_to_leaf_path = nx.shortest_path(self.T_Seg, source=self.root, target=l)
            for i, parent in enumerate(root_to_leaf_path):
                
                if i >= len(root_to_leaf_path)-1:
                    break 
                child = root_to_leaf_path[i+1]
           
                if self.T_Seg[parent][child]["event"] == "mutation":
                    mut_edges.append((parent, child))
                    k = self.T_Seg[parent][child]["cluster"]
                    dcf_vals[k] = self.T_Seg[parent][child]["dcf"]
 
            # Sort the dictionary by value in descending order
            sorted_dcfs = dict(sorted(dcf_vals.items(), key=lambda item: item[1], reverse=True))
  
            #order the edges by decreasing DCF order and update the assigned cluster and dcf values
            for k, edge in zip(sorted_dcfs, mut_edges):
                u,v = edge 
                self.T_Seg[u][v]["cluster"]=k 
                self.T_Seg[u][v]["dcf"] = sorted_dcfs[k]
                self.mut_mapping[v] = self.snvs_by_cluster[k]
            
    



  



                





        
    
    def fit(self):
        
        for k in np.unique(self.clusters):
            
            state_tree_snvs_pairs = self.get_unique_state_trees(self.snvs_by_cluster[k])
            if len(state_tree_snvs_pairs) > 0:
                id, _ = state_tree_snvs_pairs[0]
                node_id = self.add_mutation_edge(id, k)
                #   self.mut_mapping[node_id] = muts
       
        
        self.order_mutation_edges()
        self.cell_mapping = {n: [] for n in self.T_Seg}
        cell_assign = self.map_assign_cells()
        print(self)
        return cell_assign


 

    

     

    def map_assign_cells(self):
 
        alpha = np.full(self.m, self.alpha).reshape(1,-1)
        clone_order = list(nx.dfs_preorder_nodes(self.T_Seg, source=self.root))
        total_cn_states = {}
        like_list =[]
        for n in clone_order:
            x,y = self.T_Seg.nodes[n]["genotype"]
            total_cn = x+y 
            total_cn_states[n]= total_cn
            cna_geno = np.full(self.m, total_cn, dtype=int).reshape(1,-1)
            clone = nx.shortest_path(self.T_Seg, source=self.root, target=n)
            # if self.T_Seg[clone[-2]][clone[-1]]["event"]== "mutation":
            snvs = []
            for c in clone:
                if c in self.mut_mapping:
                    snvs += self.mut_mapping[c].tolist()
                
            y_vec = np.zeros(shape=self.m, dtype=int)
           
            y_vec[snvs] = 1
            y_vec = y_vec.reshape(1,-1)

            # cell_by_snv_like = np.zeros((self.data.var.shape[0],self.m))
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
            cell_by_snv_like = like_func(self.data.var[:, self.snvs], self.data.total[:,self.snvs], y_vec, cna_geno, alpha)
         
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
        for s in self.cells_by_sample:
            cn = self.total_cn_by_sample[s]
            clone_filter = [total_cn_states[n]==cn for n in clone_order]

            clones = np.array([clone_order[i] for i in range(len(clone_filter)) if clone_filter[i]],dtype=int)
            
            cells = self.cells_by_sample[s]
            sample_likes = cell_likes[clone_filter][:, cells]

            like = np.max(sample_likes, axis=0).sum()
            self.likelihood += like
            map_clone = np.argmax(sample_likes, axis=0)
            for c, m in zip(cells, map_clone):
                cell_assign[c]= clones[m]
                self.cell_mapping[clones[m]].append(c)
            

      
      
        print(f"Log Likelihood: {self.likelihood}")
        for n in self.T_Seg:
            print(f"Node{n}: {len(self.cell_mapping[n])} cells assigned")

        self.cell_mapping = {n: np.array(self.cell_mapping[n], dtype=int) for n in self.T_Seg}
    
        return cell_assign
 
            






        




    def add_mutation_edge(self, id, cluster):
        tree = self.id_to_tree[id]
                    
 
        next_node  = max([n for n  in self.T_Seg]) + 1
        edge_dcf = tuple([d for d in self.DCF[:,cluster]])
        m_0, m_m = tree.find_split_pairs()

        u, u_geno = m_0
        v, v_geno = m_m
     
     
        leaves = [n for n in self.T_Seg if self.T_Seg.out_degree[n]==0]
        for l in leaves:
            root_to_leaf_path = nx.shortest_path(self.T_Seg, source=self.root, target=l)
            #does  the mutation edge occur within or after the path
            start = root_to_leaf_path[0]
            end = root_to_leaf_path[-1]
            start_geno =   self.T_Seg.nodes[start]["genotype"]
            end_geno = self.T_Seg.nodes[end]["genotype"]
            node_path, snv_genotypes = tree.find_path(start_geno, end_geno)

            #insert 
            if v_geno in snv_genotypes:
                self.T_Seg.add_node(next_node, genotype=(v_geno[0], v_geno[1] ))
      

                if v_geno != snv_genotypes[-1]:
                    for i,s in enumerate(snv_genotypes):
                        if s ==v_geno:
                            p_x, p_y, p_m= snv_genotypes[i-1]
                            c_x, c_y, c_m = snv_genotypes[i+1]
                            break
            
                    for i,parent in enumerate(root_to_leaf_path):
                        if i >= len(root_to_leaf_path) -1:
                            break
                        cp_x, cp_y=  self.T_Seg.nodes[parent]["genotype"]
                        child_node = root_to_leaf_path[i+1]
                        cc_x, cc_y =self.T_Seg.nodes[child_node]["genotype"]

                        if p_x == cp_x and p_y==cp_y and cc_x == c_x and cc_y==c_y:
                            etype = self.T_Seg[parent][child_node]["event"]
                    
                            self.T_Seg.remove_edge(parent, child_node)
                            self.T_Seg.add_edge(parent, next_node, event="mutation", cluster=cluster, dcf =edge_dcf)
                  
                            self.T_Seg.add_edge(next_node, child_node, event=etype)
                            return next_node
                else:
                    parent = root_to_leaf_path[-1]
                    self.T_Seg.add_edge(parent, next_node,event="mutation", dcf =edge_dcf, cluster=cluster)
                    return next_node 
     
        else:
            self.T_Seg.add_node(next_node, genotype=(v_geno[0], v_geno[1] ))
            self.T_Seg.add_edge(self.root, next_node,event="mutation", dcf =edge_dcf, cluster=cluster)
            return next_node

        #if it wasn't added anywhere add as a child of the root

            


  




        # if len(self.T_Seg) ==0:
                
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
         
        
            
        
                        
                        
                    
    





            



