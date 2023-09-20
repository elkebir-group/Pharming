import networkx as nx
from dataclasses import dataclass
import numpy as np
from scipy.stats import binom
from scipy.stats import beta
from copy import deepcopy
import itertools
import cnatrees
EPSILON = -10e10

class BaseTree:
    def __init__(self,  tree, node_mapping, id=0) -> None:
        self.tree: nx.DiGraph = tree
        self.node_mapping: dict = node_mapping
        for key, val in self.node_mapping.items():
            self.tree.nodes[val]["genotype"] = key
        self.id = id
        for n in self.tree:
            if self.tree.in_degree[n] ==0:
                self.root = n 
                break
        self.node_to_geno = {val: key for key, val in self.node_mapping.items()}
    
    def __str__(self):
        mystr = ""
        
        for u,v in self.tree.edges:
        
            mystr += f"Node {u}: {self.tree.nodes[u]['genotype']} -> Node {v}: {self.tree.nodes[v]['genotype']}\n"
        return(mystr)
    
    def preorder(self):
        return list(nx.dfs_preorder_nodes(self.tree, source=self.root))
    
    def node_genotypes(self):
        
        return [self.tree.nodes[n]["genotype"] for n in self.tree]
    
    def node_desc_genotypes(self, u):
        geno_list = []
        if self.tree.out_degree[u] > 0:
   
          
            for v in self.tree.neighbors(u):
                geno_list.append(self.node_to_geno[v])
        return geno_list

    def edge_genotypes(self):
        edge_genos =[]
        for u,v in self.tree.edges:
            u_g = self.tree.nodes[u]["genotype"]
            v_g = self.tree.nodes[v]["genotype"]
            
            edge_genos.append((u_g, v_g))
        return set(edge_genos)
    
    def set_id(self, id):
        self.id = id 
    
    def is_leaf(self,node):
        return self.tree.out_degree[node] ==0
    

class CNATree(BaseTree):
    '''
    tree: nx.DiGraph
    node_mapping: provides a mapping CNA gentoypes (key) and node ids 
    '''
    def __init__(self, tree, node_mapping, id=0):
        super().__init__(tree, node_mapping, id)
    
    
    @staticmethod
    def assign_permutations(tree, node, current_genos):
            node_mapping_list =[]
            # Get the descendant nodes of the input node
       
            descendants =  list(nx.dfs_preorder_nodes(tree, node))
            descendants = [n for n in descendants if n != node]
          

            # Generate all permutations of [0, 1] with the length equal to the number of descendants
            permutations = list(itertools.product([0, 1], repeat=len(descendants)))
            for perm in permutations:
        
                geno = {key: val for key,val in current_genos.items()}
                for p,n in zip(perm, descendants):
          
                    v_x, v_y, v_m = current_genos[n]
                    parent = list(tree.predecessors(n))[0]
                    u_x, u_y, u_m = geno[parent]
                    #the non-mutated copy is amplified and m remains at u_m
                    if p ==0:
                        geno[n] = (v_x, v_y, u_m)
                    else:
                      #the mutated copy is amplified
                      diff = max([v_x- u_x, v_y-u_y])
                      geno[n] = (v_x, v_y, u_m+diff)
                nm = {val: key for key,val in geno.items()}
                node_mapping_list.append(nm)
            return node_mapping_list
    def enumerate_genotype_trees(self):
        tcna = list(self.edge_genotypes())
        T_SNVS = []
        if len(tcna) > 0:
            tsnvs = cnatrees.get_genotype_trees(tcna)
         
            for j,t in enumerate(tsnvs):
                geno_tree = nx.DiGraph()
                geno_tree.add_edges_from(t)
                node_mapping = {}
                for i,u in enumerate(geno_tree):
                    node_mapping[u] = i
                geno_tree = nx.relabel_nodes(geno_tree, node_mapping)
                T_SNVS.append(SNVTree(geno_tree, node_mapping, id=j))
        else:
            geno_tree = nx.DiGraph()
            geno_tree.add_edge(0, 1)
            node_mapping = {}
            node_mapping[(1,1,0,0)] = 0
            node_mapping[(1,1,1,0)] = 1
            T_SNVS.append(SNVTree(geno_tree, node_mapping,0))
        return T_SNVS

    def enumerate_snv_trees(self):
        '''
        returns list of list of tree clusters by tree id  and full list of TSNVs
        '''
    

        tree_clusters = []
        T_SNVS = []

        node_mapping = {}
        for geno in self.node_mapping:
            x,y = geno
            node_mapping[(x,y,0)] = self.node_mapping[geno]
        genos = {val: key for key, val in node_mapping.items()}
        id = 0
        

        num_nodes = len(self.tree)
        for node in self.tree:
            x,y,m = genos[node]
            #always add a mutation edge as the direct descendent of node
            children = list(self.tree.neighbors(node))
            for i in range(len(children)+1):
            
                n_id = num_nodes
                snv_node_map = deepcopy(node_mapping)
   
                snv_node_map[(x,y,1)] =n_id
                genos[n_id]=(x,y,1)


                if i ==0:
                    snv_tree =  nx.DiGraph()
                    snv_tree.add_edges_from(self.tree.edges())
                    snv_tree.add_edge(node, n_id)
                    T_SNV = SNVTree(snv_tree, snv_node_map,id)
                    # print(id)
                    # print(T_SNV)
                    T_SNVS.append(T_SNV)
                    tree_clusters.append( [id])
                    id += 1
                else:
                    for kids in itertools.combinations(children, i):
                        #insert a mutation edge before this subet of children
                        snv_tree =  nx.DiGraph()
                        snv_tree.add_edges_from(self.tree.edges())
                        snv_tree.add_edge(node, n_id)

                        for v in kids:
                          
                            snv_tree.remove_edge(node, v)
                            snv_tree.add_edge(n_id, v)
                        
                        node_mapping_list = self.assign_permutations(snv_tree,n_id, genos)
                        clust_tree = []
                        for nm in node_mapping_list:
                            # print(id)
             
                            tree= SNVTree(snv_tree, deepcopy(nm),id)
                            # print(tree)
                            T_SNVS.append(deepcopy(tree))
                            clust_tree.append(id)
                            id+=1
                        tree_clusters.append(clust_tree)
            


             
                                    
        return tree_clusters, T_SNVS            
                

                           




    def get_nodes_by_cn(self, cn):
        return [self.node_mapping[(x,y)] for x,y in self.node_mapping if x+y==cn]
    
class SNVTree(BaseTree):
    def __init__(self, tree, node_mapping, id):
        super().__init__(tree, node_mapping, id)  # Call the parent class constructor
        

   
        self.EPSILON: float= 10e-10
        split_pair = self.find_split_pairs()
        if split_pair is not None:
            u,v = split_pair
            self.split_node_parent, self.split_node_parent_geno = u
            self.split_node, self.split_geno = v
            self.x_star,self.y_star,x_bar, y_bar = self.split_geno
            self.m_star =x_bar + y_bar
            self.cn_star = self.x_star + self.y_star
            self.gamma = self.node_genotypes()
            self.desc_genotypes = self.get_desc_genotypes(self.split_node)

        


    

    
    def generate_CNA_tree(self):

        node_id = -1
        T_CNA = nx.DiGraph()
            
        node_mapping = {}
        for u in nx.dfs_preorder_nodes(self.tree, source=self.root):

            u_x, u_y, u_m  = self.tree.nodes[u]["genotype"]
            if (u_x, u_y) not in node_mapping:
                node_id +=1
                u_node = node_id
                T_CNA.add_node(u_node, genotype=(u_x, u_y,0))
                node_mapping[(u_x, u_y)] = u_node
         
            else:
                u_node = node_mapping[(u_x,u_y)]
            for v in self.tree.neighbors(u):
                v_x, v_y, v_m  = self.tree.nodes[v]["genotype"]
                if (v_x, v_y) not in node_mapping:
                    node_id +=1
                    v_node =node_id
                    T_CNA.add_node(v_node, genotype=(v_x, v_y,0))
                    T_CNA.add_edge(u_node, v_node)

                    node_mapping[(v_x, v_y)] = v_node
        return CNATree(T_CNA, node_mapping=node_mapping)

    def find_split_pairs(self):
        split_cand =[]
        for u, geno in self.tree.nodes("genotype"):
            x, y, x_bar, y_bar = geno
            m = x_bar + y_bar
            if m==0 and geno not in split_cand:
                split_cand.append((u,geno))
            else:
                for v, v_geno in split_cand:
                    v_x, v_y, v_x_bar, v_y_bar = v_geno
                    if x==v_x and y==v_y:
                        return (v,v_geno), (u, geno)
    

    def occurs_within( self, start_cna_geno, target_cna_geno):

        visited = set()  # Set to keep track of visited nodes
        x,y = start_cna_geno
        start_node = self.node_mapping[(x,y,0,0)]
        stack = [start_node]  # Stack to store nodes to visit
        split_node_passed = False 
        while stack:
            node = stack.pop()  # Get the next node from the stack
            if node == self.split_node:
                split_node_passed = True 
            x,y,_,_ = self.tree.nodes[node]["genotype"]
            if target_cna_geno == (x,y):
            
                # Terminate the search when the target node is found
                # print("Target node found!")
                return split_node_passed

            visited.add(node)  # Mark the node as visited

            # Add unvisited neighbors to the stack
            neighbors = self.tree.neighbors(node)
            unvisited_neighbors = [n for n in neighbors if n not in visited]
            stack.extend(unvisited_neighbors)

        # If the loop finishes without finding the target node
        return False 


 
        

    def get_desc_genotypes(self, node=None, geno=None):
        if node is not None:
            node_id = node
        elif geno is not None:
            node_id = self.node_mapping[geno]
        else:
            _,v = self.find_split_pairs()
            node_id, _ = v
   
        desc_nodes = nx.dfs_preorder_nodes(self.tree, source=node_id)
        return [self.tree.nodes[n]["genotype"] for n in desc_nodes if n != node_id]

    def is_identical(self, tree):
        if type(tree) is SNVTree:
            nodes1 = sorted(self.node_genotypes())

 
            nodes2 = sorted(tree.node_genotypes())

            if nodes1 != nodes2:
                return False

            edges1 = sorted(self.edge_genotypes())
            edges2 = sorted(tree.edge_genotypes())

            if edges1 != edges2:
                return False

            return True
        else:
            return False 
    
 


    def find_path(self, parent_cna_geno, child_cna_geno):
        p_x, p_y = parent_cna_geno
        c_x, c_y = child_cna_geno
        parent_node = self.node_mapping[(p_x, p_y,0)]
        for x,y,m in self.node_mapping:
            if x==c_x and y==c_y:
                child_node = self.node_mapping[(x,y,m)]
        shortest_path = nx.shortest_path(self.tree, source=parent_node, target=child_node)
        genotypes_in_shortest_path = [self.tree.nodes[node]["genotype"] for node in shortest_path]
        return shortest_path, genotypes_in_shortest_path

    #is self a refinment of a CNATree object   
    def is_refinement(self, cna_tree):
        if isinstance(cna_tree, CNATree):
            #every cna gentoype in tree 2 must be in tree 1
            tree1_genotypes = self.node_genotypes()
            tree1_cna_genotypes = set([(x,y) for x,y,_ in tree1_genotypes])
            for node, data in cna_tree.tree.nodes(data=True):
                n_x, n_y, _ = data["genotype"]
                if (n_x, n_y) not in tree1_cna_genotypes:
                    return False
        

            #if we remove mutation edges where the CNA does change
            #the trees should be identical
            t1_edges  = self.edge_genotypes(cna=True)
            t1_filt= set([(u,v) for u,v in t1_edges if u !=v])
    
            t2_edges  = cna_tree.edge_genotypes(cna=True)
            return(t1_filt==t2_edges)
        else:
            return False 
    
    def get_nodes_by_cn(self, cn):
        return [self.node_mapping[(x,y,x_bar, y_bar)] for x,y,x_bar, y_bar in self.node_mapping if x==cn[0] and y==cn[1]]


    def likelihood(self, cells_by_cn, a,d, alpha):
        '''
        cells_by_cn: a dict of lists (np:arrays?) that map the cn number to the cells with that total_cn
        a, d: vectors representing the column in the data of a particular snv

        '''

        all_cells =[]
        for cn ,cells in cells_by_cn.items():
            all_cells += cells.tolist()
        all_cells.sort()
        # index_to_cells= {i: c for i,c in enumerate(all_cells)}
        # cells_to_index = {c: i for i,c in enumerate(all_cells)}
        #for each node in the tree, calculate the likelihood of cells being assign to that node
        # m0,m1 = self.find_split_pairs()
        node_assign = {}
        # mut_node, mut_geno = m1
        # desc_genos = self.get_desc_genotypes(mut_node)
        total_likelihood = 0
        for cn, cells in cells_by_cn.items(): 

            # cn_indices = [cells_to_index[c] for c in cells]
            a_vec, d_vec  = a[cells],d[cells]
            cand_nodes = self.get_nodes_by_cn(cn)
            like_list = []
            for n in cand_nodes:
                x,y,x_bar, y_bar = self.tree.nodes[n]["genotype"]
                m = x_bar + y_bar
                vaf = m/(x+y)
                adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
          
                like_list.append(binom.logpmf(a_vec, d_vec,adjusted_vaf))
            cell_like_by_cand_node = np.vstack(like_list)
            assign = cell_like_by_cand_node.argmax(axis=0)
            cn_like = cell_like_by_cand_node.max(axis=0).sum()
            total_likelihood += cn_like
            for c,ass in zip(cells, assign):
                node_assign[c] = cand_nodes[ass]
        return total_likelihood, node_assign
            
        


    def posterior_dcf(self, dcf,a,d,cn ):

        if d==0:
            posterior = self.EPSILON
        else:
            
            cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
            #v = self.dcf(dcf, cn)
            #copy code so that posterior 
            v = (dcf*self.m_star)/cn \
                + (1/cn)*sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes])
    
            posterior = max(beta.logpdf(v, a+1, d-a + 1), self.EPSILON)
        return posterior
        


    def vectorized_posterior(self, dcf_vec, a_vec, d_vec, cn):
         alpha =0.001
         cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
         const_part = sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes])
    
         vaf = dcf_vec*self.m_star/cn + const_part/cn
         if vaf> 1:
             print(self)
         adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
         logpost = binom.logpmf(a_vec, d_vec, adjusted_vaf)
        #  logpost = beta.logpdf(v_vec, a_vec+1, d_vec-a_vec +1)
         if np.any(logpost==np.NINF):
             print("contains -NINF")
        #  logpost[logpost==np.NINF] = EPSILON

         return logpost 
        # Call the non-static function for a single element


    def dcf_to_v(self,dcf,cn):
       
        cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
        v = (dcf*self.m_star)/cn \
            + (1/cn)*sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes])
    
        if v < 0:
            return 0
        elif v > 1:
            return 1
        
        return v
    
    

    def v_to_dcf(self,v,cn, trunc=True):
        cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,_,_ in self.gamma}
        if cn == self.cn_star:
            d = (1/self.m_star)*(v*cn -sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes]))
        else:
            d = sum([cn_prop[x+y] for x,y,_, _ in self.desc_genotypes if x+y==cn])
        if trunc:
            if d < 0: return 0
            elif d > 1: return 1
        return d
    


    def v_minus(self, cn):
        cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
        return (1/cn)*sum([ m*cn_prop[x+y] for x,y,m in self.desc_genotypes if x==self.x_star and y==self.y_star])
    
    def v_plus(self, cn):
        cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
        return (self.v_min(cn) + self.m_star * cn_prop[self.x_star + self.y_star]/cn)

 


    












    
    


    