import networkx as nx
import numpy as np
from scipy.stats import binom
from scipy.stats import beta
EPSILON = -10e10
from genotype import genotype, CNAgenotype 


 


    
class GenotypeTree:
    def __init__(self, edges, id=0):
        """
        list edges: list of edges in the form of genotype tuples
        int: id: an id number for the genotype tree
        """

        self.tree = nx.DiGraph(edges)
        
        #recode the tree so the nodes are labeled by integers
        self.node_mapping:dict ={u: i for i,u in enumerate(self.tree)}
        self.tree= nx.relabel_nodes(self.tree, self.node_mapping)
        
        self.id = id


        for n in self.tree:
            if self.tree.in_degree[n] ==0:
                self.root = n 
                break
    
   
        self.id = id
        self.node_to_geno = {val: genotype(*key) for key, val in self.node_mapping.items()}


    
        
   

        split_par=  self.get_split_nodes()
        if split_par is not None:
            self.split_node_parent, self.split_node, self.split_node_parent_geno, self.split_geno  = split_par 
        
        
            self.m_star = self.split_geno.z
        
            self.gamma = [self.node_to_geno[u] for u in self.tree]
            self.desc_genotypes = [self.node_to_geno[u] for u in nx.descendants(self.tree, self.split_node)]



    def get_node_by_genotype(self, geno):
        for v in self.tree:
                if self.node_to_geno[v] == geno:
                    return v
        return None
    

    def get_node_by_cna_genotype(self, cna_geno):
        
        for v in self.tree:
            geno = self.node_to_geno[v]
            if geno.x == cna_geno.x and geno.y == cna_geno.y:
                return v
        return None
            
        
        
    def __str__(self):
        mystr = ""
        
        for u,v in self.tree.edges:
            u_geno = self.node_to_geno[u]
            v_geno = self.node_to_geno[v]
        
            mystr += f"Node {u}: {u_geno} -> Node {v}: {v_geno}\n"
        return(mystr)

    
    def get_nodes_by_cn(self, cn):
        if not isinstance(cn, CNAgenotype):
            cn = CNAgenotype(*cn)
        return [self.node_mapping[(x,y,x_bar, y_bar)] for x,y,x_bar, y_bar in self.node_mapping if x==cn.x and y==cn.y]

    def get_split_nodes(self):
        
        for u in nx.dfs_preorder_nodes(self.tree, self.root):
            if u != self.root:
                parent= list(self.tree.predecessors(u))[0]
                p_geno =self.node_to_geno[parent]
                u_geno = self.node_to_geno[u]
                if p_geno.cna_eq(u_geno) and u_geno.z > 0:
                    return parent, u, p_geno, u_geno


        


    def posterior_dcf(self, dcf,a,d,cn_prop ):
        '''
        dcf: float dcf cluster center
        a: int  total variant reads for a single SNV
        d: int  total reads for a single SNV
        cn_prop: dict   cn states with cn proportions 
        '''
        if d==0:
            posterior =  EPSILON
        else:
            

            v = self.dcf_to_v(dcf, cn_prop)
            # posterior = max(beta.logpdf(v, a+1, d-a + 1), self.EPSILON)
            posterior = max(binom.logpmf(a,d,v),EPSILON)
        return posterior
  
        


    def vectorized_posterior(self, dcf, a_vec, d_vec, cn_prop):
        '''
        dcf: int representing the cluster dcf value
        a_vec: np:array of length snvs containing the total variant reads for each SNV 
        d_vec: np:array of length snvs containing the total reads for each SNV 
        cn_prop: dict containing CN states and proportions (mu in decifer paper)

        returns np:array containing the posterior probability of the read counts for each SNV
        '''

        v = self.dcf_to_v(dcf, cn_prop)
    
        if np.isnan(v):
            print("Warning, DCF not valid for cn proportions and genotype tree!")

        logpost = binom.logpmf(a_vec, d_vec, v)
        #
        #  logpost = beta.logpdf(v_vec, a_vec+1, d_vec-a_vec +1)
        #if v (converted vaf) is outside [0,1] will return np.NINF
        if np.any(logpost==0):
             print("at least one log probability 0")

        logpost[np.isnan(logpost)] = EPSILON


        # logpost[logpost == np.NINF] = EPSILON

        return logpost 


    @staticmethod
    def compute_F(cn_prop):
        return  sum([(cn[0] + cn[1])*cn_prop[cn] for cn in cn_prop])

    def dcf_to_v(self,dcf,cn_prop):
        F = self.compute_F(cn_prop)



        v = (dcf*self.m_star)/F \
            + (1/F)*sum([(g.z-self.m_star)*
                         cn_prop[(g.x,g.y)] for g in self.desc_genotypes 
                         if not g.cna_eq(self.split_geno)])
        
        vmin =self.v_minus(cn_prop)
        vplus = self.v_plus(cn_prop, vmin)

        if v >= vmin and v <= vplus:
            return v 
        else:
            return np.NAN
        

    
    

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
    


    def v_minus(self, cn_prop):
        F= self.compute_F(cn_prop)
        return sum(g.z*cn_prop[(g.x, g.y)] for g in self.gamma if not g.cna_eq(self.split_geno))/ F

    
    def v_plus(self, cn_prop, vmin=None):
        F = self.compute_F(cn_prop)
        if vmin is None:
            vmin = self.v_minus(cn_prop)

        v_plus = vmin + self.m_star* cn_prop[(self.split_geno.x, self.split_geno.y)] / F 

        return v_plus 

    



 


    












    
    

   