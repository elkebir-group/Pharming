import networkx as nx
from dataclasses import dataclass


from scipy.stats import beta
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
         cn_prop = {x+y: 1.0 if x+y==cn else 0.0 for x,y,m in self.gamma}
         const_part = sum([(m-self.m_star)*cn_prop[x+y] for x,y,m in self.desc_genotypes])
    
         v_vec = dcf_vec*self.m_star/cn + const_part
         logpost = -1*beta.logpdf(v_vec, a_vec, d_vec, cn)
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


    












    
    


    