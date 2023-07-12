import networkx as nx
from dataclasses import dataclass
import numpy as np
from scipy.stats import binom
from scipy.stats import beta

EPSILON = -10e10
@dataclass
class GenotypeTree:

    tree: nx.DiGraph
    node_mapping: dict
    id: int = 0
    EPSILON: float= 10e-10

    # def __init__(self, tree:nx.DiGraph=None, node_mapping: dict=None, id=0) -> None:


    def __post_init__(self):


        #dynamically find the root node
        for n in self.tree:
            if self.tree.in_degree[n] ==0:
                self.root = n 
                break
        split_pair = self.find_split_pairs()
        if split_pair is not None:
            u,v = split_pair
            self.split_node, self.split_geno = v
            self.x_star,self.y_star,self.m_star  = self.split_geno
            self.cn_star = self.x_star + self.y_star
            self.gamma = self.get_node_genotypes()
            self.desc_genotypes = self.get_desc_genotypes(self.split_node)

        

    def __str__(self):
        mystr = ""
        
        for u,v in self.tree.edges:
        
            mystr += f"Node {u}: {self.tree.nodes[u]['genotype']} -> Node {v}: {self.tree.nodes[v]['genotype']}\n"
        return(mystr)
    
    def preorder(self):
        return list(nx.dfs_preorder_nodes(self.tree, source=self.root))
    
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
        return GenotypeTree(T_CNA, node_mapping=node_mapping)

    def find_split_pairs(self):
        split_cand =[]
        for u, geno in self.tree.nodes("genotype"):
            x, y, m = geno
            if m==0 and geno not in split_cand:
                split_cand.append((u,geno))
            else:
                for v, v_geno in split_cand:
                    v_x, v_y, v_m = v_geno
                    if x==v_x and y==v_y:
                        return (v,v_geno), (u, geno)
    
    def get_node_genotypes(self):
        
        return [self.tree.nodes[n]["genotype"] for n in self.tree]

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


    
    def get_edge_genotypes(self, cna=False):
        edge_genos =[]
        for u,v in self.tree.edges:
            u_g = self.tree.nodes[u]["genotype"]
            v_g = self.tree.nodes[v]["genotype"]
            
            edge_genos.append((u_g, v_g))
        
        if cna:
            edge_genos =[((u_g[0], u_g[1]), (v_g[0], v_g[1])) for u_g, v_g in edge_genos]

        return set(edge_genos)
    
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

    #is self a refinment of tree2?   
    def is_refinement(self, cna_tree):

        #every cna gentoype in tree 2 must be in tree 1
        tree1_genotypes = self.get_node_genotypes()
        tree1_cna_genotypes = set([(x,y) for x,y,_ in tree1_genotypes])
        for node, data in cna_tree.tree.nodes(data=True):
            n_x, n_y, _ = data["genotype"]
            if (n_x, n_y) not in tree1_cna_genotypes:
                return False
    

        #if we remove mutation edges where the CNA does change
        #the trees should be identical
        t1_edges  = self.get_edge_genotypes(cna=True)
        t1_filt= set([(u,v) for u,v in t1_edges if u !=v])
  
        t2_edges  = cna_tree.get_edge_genotypes(cna=True)
        return(t1_filt==t2_edges)
    
    def get_nodes_by_cn(self, cn):
        return [self.node_mapping[(x,y,m)] for x,y,m in self.node_mapping if x+y==cn]

    # def likelihood_function( a,d,y,c, alpha):
        
    #     if d ==0:
    #         return 0
    #     elif y ==0:
        
    #         return binom.pmf(a,d,alpha)
          
    #     else:
    #         vaf = np.arange(1, c)/c
    #         # vaf  = 1/c
    #         # vaf = 0.5
    #         adjusted_vaf =  vaf*(1- alpha) + (1-vaf)*(alpha/3)
    #         # return binom.pmf(a,d,adjusted_vaf)
    #         return (1/vaf.shape[0])*np.sum(binom.pmf(a,d,adjusted_vaf))
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
                x,y,m = self.tree.nodes[n]["genotype"]

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
        cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
        if cn == self.cn_star:
            d = (1/self.m_star)*(v*cn -sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes]))
        else:
            d = sum([cn_prop[x+y] for x,y,m in self.desc_genotypes if x+y==cn])
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

    # def calc_lamb(self, v,cn):
    #     cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
    #     return  (1/self.m_star)*(v*cn - sum([]))
    
     
  
    # def geno_prop(self, geno, v):
    #     '''
    #     we make the assumption that a sample is composed of cells with the same cn number state
    #     so that mu_{x,y} is either 1 or 0
    #     '''
    #     x, y, m = geno
    #     if x != self.x_star and y != self.y_star:
    #         return 1
    #     else:
    #         self.lamb = self.calc_lamb()


    












    
    


    