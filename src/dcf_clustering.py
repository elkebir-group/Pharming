from genotype_tree import GenotypeTree
import numpy as np
from dataclasses import dataclass
from scipy.stats import beta
from scipy.special import betaln
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
import time


TOLERANCE = 1e-03
EPSILON = -1e40
SEQERROR = 1e-40



class DCF_Clustering:
    def __init__(self, T_CNAs, T_SNVs, clusters= 3, nrestarts=10, rng=None):
        self.nrestarts =nrestarts 
        if rng is None:
            self.rng = rng 
        else:
            self.rng = np.random.default_rng(1026)
        # self.clusters= [i+1 for i in range(max_clusters)]
        
        self.k = clusters
        self.T_CNAs = T_CNAs

        self.T_SNVs = T_SNVs
    
    

    def init_cluster_centers(self):
        ''' 
        randomly initialize cluster centers 
        '''
        return self.rng.uniform(low=0.0, high=1.0, size=(self.nsamples, self.k))


    def optimize_snv_assignments(self, T_SNVs, dcfs, snvs, alt, total):
        #assignment of snvs to T_SNV_tree
        T_SNV_dict= {}
        snv_lookup= {i: snvs[i] for i in range(snvs.shape[0])}
        clusters = np.full_like(snvs, fill_value=1)
        CNs = list(alt.keys())
        #for each homogenous sample
        combos = []
        for k in range(self.dcfs.shape[1]):
                all_arrs = []
                for tree in T_SNVs:
                    combos.append((k,tree.id))

                    post_arr =[]
                    for sample, cn in enumerate(CNs):
                        a_vec = alt[cn]
                        d_vec=total[cn]
                 
                        dcf = dcfs[sample,k]
        
                        post_arr.append(tree.vectorized_posterior(dcf, a_vec,d_vec, cn))
                        all_arrs.append(np.vstack(post_arr).sum(1))
        cand_assign = np.vstack(all_arrs)
        likelihood = cand_assign.max(axis=1).sum()
        combo_assign = np.argmax(all_arrs, axis=1)
        for i,c in enumerate(combo_assign):
            clusters[i], T_SNV_dict[snv_lookup[i]] = combos[c]
        return clusters, T_SNV_dict, likelihood

                






        # for i in range(self.k):
        #    for tree in T_SNVs:

            #compute the probability of the given cluster dcf given tree t and the read counts
            #compute the likelihood
            
     
   

    def optimize_cluter_centers(self):
        pass 

    def compute_likelihood(self):
        pass 

    def fit_cna_tree(self, T_CNA,  snvs, alt, total):


        cand_genotype_trees = self.filter_genotype_trees(T_CNA)
        for j in range(self.max_iterations):

            snv_clusters, T_SNVs, likelihood = self.optimize_snv_assignments(cand_genotype_trees, dcfs, snvs,alt, total)
            dcfs = self.optimize_cluster_centers(snv_clusters, T_SNVs)
            likelihood = self.compute_likelihood(dcfs, snv_clusters, T_SNVs)
            if likelihood > best_likelihood:
                best_likelihood = likelihood
        return DCF_Data(best_likelihood, snv_clusters, T_SNVs, T_CNA, dcfs)

    def fit(self, snvs, alt, total):
        self.nsamples = len(alt)
        CNs = list(alt.keys())
        for i in range(self.nrestarts):
            
            dcfs = self.init_cluster_centers()
            best_likelihood = np.NINF
    
            cna_tree_likes = []
            res = []
            for cna_tree in self.T_CNAs:
                cand_T_SNVs = [tree for tree in self.T_SNVs if tree.is_refinement(cna_tree)]
                tree_results= self.fit_CNA(cand_T_SNVs, snvs, alt, total)
                cna_tree_likes.append(tree_results.likelihood)
                res.append(tree_results)
            max_index = cna_tree_likes.index(max(cna_tree_likes))

        if best_likelihood > max(cna_tree_likes):
            best_likelihood = max(cna_tree_likes)
            best_res = res[max_index]
        return best_res
            
                

            
@dataclass 
class DCF_Data:
    likelihood: float 
    clusters: np.array 
    T_SNVs: dict
    T_CNA: GenotypeTree 
    DCFs: np.array


        
#equation 28
#purity rho = 1 always
#Fraction copy number of a sample = cn
# def convert_dcf_to_vaf( dcf,cn, a, d, tree):
#     v, u = tree.find_split_pairs()
   
#     tree.get_descendents()
#     vaf = d*m_star/cn
    
#cn prop is a dictionary
# def posterior_dcf(dcf, cn, alt, total, T_SNV):
#     #equation 15
  
#     v, u = T_SNV.find_split_pairs()
#     split_node, split_geno = v
#     x_star,y_star,m_star  = split_geno
#     gamma = T_SNV.get_edge_genotypes()
#     cn_prop = {x+y: x+y==cn for x,y,m in gamma}


     

#     desc_genotypes = T_SNV.get_desc_genotypes(split_node)


#     def v_minus():
#         return (1/cn)*sum([ m*cn_prop[x+y] for x,y,m in desc_genotypes if x==x_star and y==y_star])
    
#     def v_plus():
#         return v_minus() + (m_star * cn_prop[x_star + y_star]/cn)
  
#     return v_minus(), v_plus()
#     # def convert_dcf_to_vaf():
       


    



# def posterior_c(c, a, d, config, sample, mut=None):
#     try:
#         v = config.c_to_v(c, sample)
#     except IndexError:
#         print("Mutation {}".format(mut.label))
#         raise
#     if v:
#         post = max(beta.pdf(v, a[sample]+1, (d[sample]-a[sample]) + 1), EPSILON)
#     else:
#         post = EPSILON

#     return math.log(post)



