from genotype_tree import GenotypeTree
import numpy as np
from dataclasses import dataclass


class DCF_Clustering:
    def __init__(self, clusters= 2, nrestarts=10, rng=None, tol=):
        self.nrestarts =nrestarts 
        if rng is None:
            self.rng = rng 
        else:
            self.rng = np.random.default_rng(1026)
        self.k= clusters 
        self.likelihood = np.NINF
    

    def initilize_cluster_centers(self):
        ''' 
        randomly initialize cluster centers 
        '''
        pass 


    def optimize_snv_assignments(self):
        pass 

    def optimize_cluter_centers(self):
        pass 

    def compute_likelihood(self):
        pass 

    def fit_cna_tree(T_CNA):


        cand_genotype_trees = self.filter_genotype_trees(T_CNA)
        for j in range(self.max_iterations):
         
            snv_clusters, T_SNVs = self.optimize_snv_assignments(cand_genotype_trees, dcfs)
            dcfs = self.optimize_cluster_centers(snv_clusters, T_SNVs)
            likelihood = self.compute_likelihood(dcfs, snv_clusters, T_SNVs)
            if likelihood > best_likelihood:
                best_likelihood = likelihood
        return DCF_Data(best_likelihood, snv_clusters, T_SNVs, T_CNA, dcfs)

    def fit(self, data):
        for i in range(self.nrestarts):
            
            dcfs = self.initialize_cluster_centers()
            best_likelihood = np.NINF
            cna_likelihood = np.NINF
            cna_tree_likes = []
            res = []
            for cna_tree in cna_trees:
                tree_results= self.fit_CNA(cna_tree)
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


        

