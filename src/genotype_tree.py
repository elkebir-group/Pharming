import networkx as nx
from dataclasses import dataclass
import numpy as np
from scipy.stats import binom
from scipy.stats import beta
EPSILON = -10e10
from clonal_tree import ClonalTree
from genotype import genotype, CNAgenotype 

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
    
 

    
    


    
class GenotypeTree(ClonalTree):
    def __init__(self, tree, node_mapping, snv=-1, id=0):
        self.snv  = snv
        genotypes = {v: {} for v in tree}
        for geno, v in node_mapping.items():
            genotypes[v][snv] = genotype(*geno)
        self.node_mapping: dict = node_mapping
   
        self.id = id
        self.node_to_geno = {val: key for key, val in self.node_mapping.items()}
    
        

        super().__init__(tree, genotypes, key=self.snv)  # Call the parent class constructor
        
        for key, val in self.node_mapping.items():
            self.tree.nodes[val]["genotype"] = key
   
        self.EPSILON: float= 10e-10
        split_pair = self.find_split_pairs()
        if split_pair is not None:
            u,v = split_pair
            self.split_node_parent, self.split_node_parent_geno = u
            self.split_node, self.split_geno = v
            # self.x_star,self.y_star,x_bar, y_bar = self.split_geno
            self.m_star = self.split_geno.z
            self.cn_star = self.split_geno.w
            self.gamma = self.node_genotypes()
            self.desc_genotypes = self.get_desc_genotypes(self.split_node)


    def node_genotype(self, u):
        return self.genotypes[u][self.snv]
    
    def set_node_genotype(self, u, geno):
        self.genotypes[u][self.snv] = geno

    def get_node_by_genotype(self, geno):
        for v in self.tree:
                if self.node_genotype(v) == geno:
                    return v
        return None
    

    def get_node_by_cna_genotype(self, cna_geno):
        
        for v in self.tree:
            geno = self.node_genotype(v)
            if geno.x == cna_geno.x and geno.y == cna_geno.y:
                return v
        return None
            
        
    def node_genotypes(self):
        '''
        returns the set of node genotypes as 4-tuples
        TODO: simply return the list of genotype objects 
        '''
        return [self.node_genotype(v) for v in self.tree]
        
    def __str__(self):
        mystr = ""
        
        for u,v in self.tree.edges:
            u_geno = self.node_genotype(u)
            v_geno =self.node_genotype(v)
        
            mystr += f"Node {u}: {u_geno} -> Node {v}: {v_geno}\n"
        return(mystr)

    def relabel_snv(self, new_snv_label):
        old_key = self.snv 
        self.snv = new_snv_label 

        for v in self.tree:
            if old_key in self.genotypes[v]:
                geno = self.genotypes[v].pop(old_key)  
                self.set_node_genotype(v, geno)
        for v, snvs in self.mut_mapping.items():
            if old_key in snvs:
                self.mut_mapping[v].remove(old_key)
                self.mut_mapping[v].append(self.snv)
                break
        for v, snvs in self.mut_loss_mapping.items():
            if old_key in snvs:
                self.mut_loss_mapping[v].remove(old_key)
                self.mut_loss_mapping[v].append(self.snv)
                break
        self.key = new_snv_label


    def find_split_pairs(self):
        split_cand =[]
        for u in self.tree:
            u_geno = self.node_genotype(u)
            if u_geno.z==0 and u_geno not in split_cand:
                split_cand.append((u,u_geno))
            else:
                for v, v_geno in split_cand:
                    if u_geno.cna_eq(v_geno):
                        return (v,v_geno), (u, u_geno)
    
    

    def get_desc_genotypes(self, node=None, geno=None):
        if node is not None:
            node_id = node
        elif geno is not None:
            if not isinstance(geno, genotype):
                geno = genotype(*geno)
            node_id = self.get_node_by_genotype(geno)
        else:
            _,v = self.find_split_pairs()
            node_id, _ = v
   
        desc_nodes = nx.dfs_preorder_nodes(self.tree, source=node_id)
        return [self.node_genotype(n) for n in desc_nodes if n != node_id]


    def edge_genotypes(self):
        edge_genos =[]
        for u,v in self.tree.edges:
            u_g = self.node_genotype(u)
            v_g = self.node_genotype(v)
            
            edge_genos.append((u_g, v_g))
        return set(edge_genos)
 

    def node_desc_genotypes(self, u):
        geno_list = []
        if self.tree.out_degree[u] > 0:
            for v in self.tree.neighbors(u):
                geno_list.append(self.node_genotype(v))
        return geno_list

 
 
    

    def find_path(self, parent_cna_geno, child_cna_geno):
    
        if not isinstance(parent_cna_geno, CNAgenotype):
            parent_cna_geno = CNAgenotype(*parent_cna_geno)
        if not isinstance(child_cna_geno, CNAgenotype):
            child_cna_geno = CNAgenotype(*child_cna_geno)

        parent_node = self.get_node_by_cna_genotype(parent_cna_geno)
        child_node = self.get_node_by_cna_genotype(child_cna_geno)
        shortest_path = nx.shortest_path(self.tree, source=parent_node, target=child_node)
        genotypes_in_shortest_path = [self.node_genotype(node) for node in shortest_path]
        return shortest_path, genotypes_in_shortest_path


    
    def get_nodes_by_cn(self, cn):
        if not isinstance(cn, CNAgenotype):
            cn = CNAgenotype(*cn)
        return [self.node_mapping[(x,y,x_bar, y_bar)] for x,y,x_bar, y_bar in self.node_mapping if x==cn.x and y==cn.y]


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
            
        


    def posterior_dcf(self, dcf,a,d,cn_prop ):
        '''
        dcf: float dcf cluster center
        a: int  total variant reads for a single SNV
        d: int  total reads for a single SNV
        cn_prop: dict   cn states with cn proportions 
        '''
        if d==0:
            posterior = self.EPSILON
        else:
            

            #compute fractional copy number 
            F = sum([(cn[0] + cn[1])*cn_prop[cn]] for cn in cn_prop)
            v = (dcf*self.m_star)/F\
                + (1/F)*sum([(m-self.m_star)*cn_prop[(x,y)] for x,y,m in self.desc_genotypes])
    
            # posterior = max(beta.logpdf(v, a+1, d-a + 1), self.EPSILON)
            posterior = max(binom.logpmf(a,d,v), self.EPSILON)
        return posterior
  
        


    def vectorized_posterior(self, dcf_vec, a_vec, d_vec, cn):
         raise NotImplemented
         ######## USE non-vectorized version above####
         ### Needs to be udpated#####
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


    def dcf_to_v(self,dcf,cn_prop):
        F = sum([(cn[0] + cn[1])*cn_prop[cn]] for cn in cn_prop)

        v = (dcf*self.m_star)/F \
            + (1/F)*sum([(g.z-self.m_star)*cn_prop[g.x+g.y] for g in self.desc_genotypes])
    
        if v < 0:
            return 0
        elif v > 1:
            return 1
        
        return v
    
    

    def v_to_dcf(self,v,cn, trunc=True):
        ###### needs to be updated ###
        raise NotImplemented
        cn_prop = {g.w: 1.0 if g.w==cn else 0.0 for g in self.gamma}
        if cn == self.cn_star:
            d = (1/self.m_star)*(v*cn -sum([(g.m-self.m_star)*cn_prop[g.w] for g in self.desc_genotypes]))
        else:
            d = sum([cn_prop[g.w] for g in self.desc_genotypes if g.w==cn])
        if trunc:
            if d < 0: return 0
            elif d > 1: return 1
        return d
    


    # def v_minus(self, cn):
    #     cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
    #     return (1/cn)*sum([ m*cn_prop[x+y] for x,y,m in self.desc_genotypes if x==self.x_star and y==self.y_star])
    
    # def v_plus(self, cn):
    #     cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
    #     return (self.v_min(cn) + self.m_star * cn_prop[self.x_star + self.y_star]/cn)
    



 


    












    
    

    # def generate_CNA_tree(self):

    #     node_id = -1
    #     T_CNA = nx.DiGraph()
            
    #     node_mapping = {}
    #     for u in nx.dfs_preorder_nodes(self.tree, source=self.root):

    #         u_x, u_y, u_m  = self.tree.nodes[u]["genotype"]
    #         if (u_x, u_y) not in node_mapping:
    #             node_id +=1
    #             u_node = node_id
    #             T_CNA.add_node(u_node, genotype=(u_x, u_y,0))
    #             node_mapping[(u_x, u_y)] = u_node
         
    #         else:
    #             u_node = node_mapping[(u_x,u_y)]
    #         for v in self.tree.neighbors(u):
    #             v_x, v_y, v_m  = self.tree.nodes[v]["genotype"]
    #             if (v_x, v_y) not in node_mapping:
    #                 node_id +=1
    #                 v_node =node_id
    #                 T_CNA.add_node(v_node, genotype=(v_x, v_y,0))
    #                 T_CNA.add_edge(u_node, v_node)

    #                 node_mapping[(v_x, v_y)] = v_node
    #     return CNATree(T_CNA, node_mapping=node_mapping)
    


    # def find_split_pairs(self):
    #     split_cand =[]
    #     for u, geno in self.tree.nodes("genotype"):
    #         x, y, x_bar, y_bar = geno
    #         m = x_bar + y_bar
    #         if m==0 and geno not in split_cand:
    #             split_cand.append((u,geno))
    #         else:
    #             for v, v_geno in split_cand:
    #                 v_x, v_y, v_x_bar, v_y_bar = v_geno
    #                 if x==v_x and y==v_y:
    #                     return (v,v_geno), (u, geno)
    

        #is self a refinment of a CNATree object   
    # def is_refinement(self, cna_tree):
    #     if isinstance(cna_tree, CNATree):
    #         #every cna gentoype in tree 2 must be in tree 1
    #         tree1_genotypes = self.node_genotypes()
    #         tree1_cna_genotypes = set([(x,y) for x,y,_ in tree1_genotypes])
    #         for node, data in cna_tree.tree.nodes(data=True):
    #             n_x, n_y, _ = data["genotype"]
    #             if (n_x, n_y) not in tree1_cna_genotypes:
    #                 return False
        

    #         #if we remove mutation edges where the CNA does change
    #         #the trees should be identical
    #         t1_edges  = self.edge_genotypes(cna=True)
    #         t1_filt= set([(u,v) for u,v in t1_edges if u !=v])
    
    #         t2_edges  = cna_tree.edge_genotypes(cna=True)
    #         return(t1_filt==t2_edges)
    #     else:
    #         return False 


        # def occurs_within( self, start_cna_geno, target_cna_geno):

    #     visited = set()  # Set to keep track of visited nodes
    #     x,y = start_cna_geno
    #     start_node = self.node_mapping[(x,y,0,0)]
    #     stack = [start_node]  # Stack to store nodes to visit
    #     split_node_passed = False 
    #     while stack:
    #         node = stack.pop()  # Get the next node from the stack
    #         if node == self.split_node:
    #             split_node_passed = True 
    #         x,y,_,_ = self.tree.nodes[node]["genotype"]
    #         if target_cna_geno == (x,y):
            
    #             # Terminate the search when the target node is found
    #             # print("Target node found!")
    #             return split_node_passed

    #         visited.add(node)  # Mark the node as visited

    #         # Add unvisited neighbors to the stack
    #         neighbors = self.tree.neighbors(node)
    #         unvisited_neighbors = [n for n in neighbors if n not in visited]
    #         stack.extend(unvisited_neighbors)

    #     # If the loop finishes without finding the target node
    #     return False 


    # def is_identical(self, tree):
    #     if type(tree) is SNVTree:
    #         nodes1 = sorted(self.node_genotypes())

 
    #         nodes2 = sorted(tree.node_genotypes())

    #         if nodes1 != nodes2:
    #             return False

    #         edges1 = sorted(self.edge_genotypes())
    #         edges2 = sorted(tree.edge_genotypes())

    #         if edges1 != edges2:
    #             return False

    #         return True
    #     else:
    #         return False 