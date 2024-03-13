from genotype_tree import GenotypeTree
import numpy as np
from dataclasses import dataclass
from scipy.optimize import minimize_scalar
import clonelib
import networkx as nx 

TOLERANCE = 1e-03
EPSILON = -1e40
SEQERROR = 1e-40






def scalar_obj( dcf, id_to_index, trees, alt_vec, total_vec,cn):
        obj = 0
        for id, tree in zip(id_to_index, trees): #could loop over snv/tree matchings and do non-vectoirzed posterior
            alt = alt_vec[cn][id_to_index[id]]
            total =total_vec[cn][id_to_index[id]]
            obj += tree.vectorized_posterior(dcf, alt, total,cn).sum()
        return -1*obj


def unpack_assignments(cluster, assignments):
    tree_id_to_indices ={}
    for index, assign in enumerate(assignments):
        k, tree_id = assign
        if k == cluster:
            if tree_id in tree_id_to_indices:
                    tree_id_to_indices[tree_id].append(index)
            else:
                    tree_id_to_indices[tree_id] = [index]
    return tree_id_to_indices
class DCF_Clustering:
    def __init__(self, nrestarts=25, seed=1026, verbose=False):
        self.nrestarts =nrestarts 
   
        self.rng = np.random.default_rng(seed)

        self.verbose=verbose


        # self.clusters= [i+1 for i in range(max_clusters)]
        
        # self.k = clusters
        # self.T_CNAs = T_CNAs

        # self.T_SNVs = T_SNVs
        # self.id_to_tree = {tree.id: tree for tree in T_SNVs}
        # self.max_iterations=100
        # self.combos = []
   
        # self.T_SNV_clusters = T_SNV_clusters
        # for k in range(self.k):
        #     for t,tree in enumerate(self.T_SNVs):
        #             self.combos.append((k,t))
        
        # self.verbose = False
    

    def init_cluster_centers(self):
        ''' 
        randomly initialize cluster centers 
        '''
        #return self.rng.uniform(low=0.0, high=1.0, size=self.k)
        return [0.056, 0.146, 0.179, 0.996, 0.617, 0.138, 0.382]


    def optimize_cluster_and_tree_assignments(self, tree_id_to_indices, dcfs, alt, total, ell):
        snvs = self.data.seg_to_snvs[ell]
        cn_prob = self.data.cn_proportions(ell)
        # get a_vec and d_vec for the snvs
        a_vec = alt.sum(axis=0)[snvs]
        d_vec = total.sum(axis=0)[snvs]
        best_prob = {}
        best_cluster = {}
        best_tree = {}
        likelihood = {}
        for snv in snvs:
            best_prob[snv] = -np.inf

        for cluster in range(len(dcfs)):
            for tree in tree_id_to_indices:
                #prob_arr = []
                #for a, d in zip(a_vec, d_vec):
                prob_arr = tree.vectorized_posterior(dcfs[cluster], a_vec, d_vec, cn_prob)
                #prob_arr.append(prob)
                idx = 0
                for prob in prob_arr:
                    if prob > best_prob[snvs[idx]]:
                        best_cluster[snvs[idx]] = cluster
                        best_tree[snvs[idx]] = tree
                        likelihood[snvs[idx]] = prob
                    idx += 1

        return best_cluster, best_tree, likelihood


    def scalar_obj_new(dcf, tree_assignments, segs_in_cluster, alt, total):
        obj = 0
        # for each SNV:
        # get alt, cn, total

        for seg in segs_in_cluster:
            cn_prob = self.data.cn_proportions(seg)
            tree = tree_assignments[seg]
            for snv in self.data.seg_to_snvs[seg]:
                obj += tree.posterior_dcf(dcf, alt[snv], total[snv], cn_prob)

        return obj

    def optimize_cluster_centers(self, dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, alt, total): #TO DO: Check if this should be updated?
        new_dcfs = []
        for k in range(len(dcfs)): #looping over clusters
            segs_in_cluster = []
            tree_in_cluster = []
            #find segments of clusters
            for seg, cluster in CLUSTER_ASSIGNMENTS.items():
                if cluster == k:
                    segs_in_cluster.append(seg)
                if len(segs_in_cluster) > 0:
                    new_dcf = minimize_scalar(scalar_obj_new, args=(TREE_ASSIGNMENTS, segs_in_cluster, alt, total), method='bounded', bounds=[0,1]).x
                else:
                    new_dcf = self.rng.random()
            new_dcfs.append(new_dcf)
        return np.array(new_dcfs)


    #     return np.array(new_dcfs)
    def compute_likelihood(self, dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, alt, total):

        #scalar_obj_new(dcf, tree_assignments, segs_in_cluster, alt, total):


        likelihood = 0
        #optimize the cluster dcf for each sample and cluster
        for cluster in range(len(dcfs)):
            segs_in_cluster = []
            for seg, clust in CLUSTER_ASSIGNMENTS.items():
                if clust == cluster:
                    segs_in_cluster.append(seg)
            if len(segs_in_cluster) > 0:
                likelihood += scalar_obj_new(dcfs[cluster], segs_in_cluster, alt, total)
        return likelihood


    def decifer(self, data, k=5):
        '''
        See data.py for the data object 
        k = # of SNV clusters 
        '''
        prev_likelihood = np.NINF
        self.k = k
        dcfs = self.init_cluster_centers()
        self.data = data 
        self.segments = list(self.data.seg_to_snvs.keys())
        #enumerat valid CNA trees for each segment
        S = {}
        cn_props ={}
        for ell in self.segments:
             cn_props[ell] = self.data.cn_proportions(ell)

             #find all copy number states (x,y) and proportions in segment ell
             states = set(cn_props[ell].keys())
             #enumerate CNA trees 
             S[ell] = clonelib.get_cna_trees(states, 1, 1)
             
        #S = {\ell: [CNA trees]}
       
        #for j in range(self.max_iterations): #TO DO: make max_iterations a parameter
        for j in range(100):
            CLUSTER_ASSIGNMENTS = {}
            TREE_ASSIGNMENTS = {}
            CNA_tree ={}
            for ell in self.segments: #TO DO: Swap these segments and CNA tree assignment
                snvs = self.data.seg_to_snvs[ell]
                seg_like = -np.Inf
                for s in S[ell]:
                    #here, get optimal cluster assignments, then save the optimal CNA tree + assignments of SNVs to trees + assignments of SNVs to clusters
                    #after doing this for all the segments, then update the clusters (using assignments & dcf values)
                    scriptT = clonelib.get_genotype_trees(s)
                    T_SNVs = [GenotypeTree(edges= edge_list, id=i) for i, edge_list in enumerate(scriptT) ]
                    # for snv_edges in scriptT:

                        # T = nx.DiGraph(snv_edges)
                        # #recode the tree so the nodes are labeled by integers
                        # node_mapping ={u: i for i,u in enumerate(T)}
                        # T= nx.relabel_nodes(T, node_mapping)

                        # T_SNVs.append(GenotypeTree(T, node_mapping))
                    
                    #find an SNV tree assignment and DCF cluster assignment
                        
                    '''
                    TO DO: Optimize cluster assignments needs to be updated.
                    You will need to compute the posterior probability of each dcf q and T pair in T_SNVs for each SNV
                    Set omega equal to the tree with max posterior prob
                    Set alpha equal cluster id q with max posterior prob
                    the likelihood is the sum of the  log posterior probabilities for all optimal assignments
                    '''   
                    cluster, tree, likelihood  = self.optimize_cluster_and_tree_assignments(T_SNVs, dcfs, self.data.var, self.data.total, ell)

                    new_seg_likelihood = 0
                    for key, val in likelihood.items():
                        new_seg_likelihood += val

                    if new_seg_likelihood > seg_like:
                        CNA_tree[ell] = s
                        CLUSTER_ASSIGNMENTS[ell] = cluster
                        TREE_ASSIGNMENTS[ell] = tree
                        seg_like = new_seg_likelihood
                
                                    
        
            
            '''
            TO DO: Optimize cluster centers need to be updated to account the assignments (omega, alpha) 
            being a dictionary of dictionaries
            


            '''  
            old_dcfs = dcfs.copy()
            dcfs = self.optimize_cluster_centers(dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, self.data.var, self.data.total)
            dcfs[dcfs > 0.99] =1.0
            dcfs[dcfs < 1e-3] =0

            '''
            TO DO: The likelihood computation needs to be updated as well 
            '''
            new_likelihood = self.compute_likelihood(dcfs,CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, self.data.var, self.data.total)
            #check for covergence:
            diff = new_likelihood -prev_likelihood
            if self.verbose:
                print(f"Previous likelihood: {prev_likelihood} New likelihood: {new_likelihood} Diff: {diff}")
            if diff < 0 or abs(diff) <  TOLERANCE:
                 dcfs = old_dcfs 
                 break
            else:
                 prev_likelihood = new_likelihood
            

            return new_likelihood, dcfs, CNA_tree, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS
            
                            

    def run(self, data, k_vals= [5,6]):
       
        results = {}
        for k in k_vals:
            for i in range(self.nrestarts):
               results[k,i] = self.decifer(data, k)


        #add model selection 
               
        return #best 
        



if __name__ == "__main__":
     
     #load in pickled data object 
     from data import Data, load_from_pickle
     dat = load_from_pickle("/Users/annahart/CLionProjects/Pharming/s11_m5000_k25_l7/n1000_c0.25_e0/data.pkl")
     
     #initialize object 
     dec = DCF_Clustering(nrestarts=50,seed=21,verbose=True)
     #all_results = dec.run(dat, k_vals=[i+1 for i in range(8)])
     results = dec.run(dat, k_vals=[7])

     

