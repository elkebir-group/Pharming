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
        for id, tree in zip(id_to_index, trees):
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
        return self.rng.uniform(low=0.0, high=1.0, size=self.k)


    def optimize_snv_assignments(self, T_SNVs, dcfs, alt, total):
        #assignment of snvs to T_SNV_tree
        # T_SNV_dict= {}
        CNs = list(alt.keys())
        #for each homogenous sample
        # combos = []
        all_arrs = []
        snvClust_toTree = {}
        total_like = 0
        tree_assigns ={}
        tree_like = {}
        for k in range(dcfs.shape[1]):
     
     
               
                best_like = np.NINF
                for j, tree_clust in enumerate(self.T_SNV_clusters):
                    tree_clust_arrs = []
                    for id in tree_clust:
                        tree= self.id_to_tree[id]
                        # combos.append((k, j, id))
                        # tree_combos.append(j, id)

                        post_arr =[]
                        for sample, cn in enumerate(CNs):
                            a_vec = alt[cn]
                            d_vec=total[cn]
                            dcf = dcfs[sample,k]
                            post_arr.append(tree.vectorized_posterior(dcf, a_vec,d_vec, cn))
                        sample_arr = np.vstack(post_arr).sum(axis=0)
                        sample_arr= sample_arr.reshape(1,-1)

                        tree_clust_arrs.append(sample_arr)

                    #pick best tree for that cluster
                    temp = np.vstack(tree_clust_arrs).max(axis=0)
                    argtemp = np.vstack(tree_clust_arrs).argmax(axis=0)
                    tree_assigns[(k,j)] = temp.copy()
                    # [tree_clust[a] for a in argtemp]
                    # if temp.sum() > best_like:
                    #      snvClust_toTree[k]=j
                    #      best_like = temp.sum()
                    #      tree_like[k] = temp.copy()
                         #this is the tree id of the best tree in the cluster
        
                # all_arrs.append(tree_like[k])

                
                
                    
        
              
        clust_assign = np.vstack(all_arrs)
        likelihood = clust_assign.max(axis=0).sum()
        combo_assign = np.argmax(clust_assign, axis=0)
        assignments = []
        #assignments is list of length snvs containing a two-tuple of cluster and tree id to which the 

        for i, c in enumerate(combo_assign):
    
             assignments.append((c, tree_assigns[c][i]))
             

        # assignments = [combos[c] for c in combo_assign ]
 
        return assignments, likelihood

                
    def optimize_cluster_assignments(self, tree_id_to_indices, dcfs, alt, total):
        CNs = list(alt.keys())
    #     #for each homogenous sample
    #     combos = []
        like = 0
        assignments = []
        clusters = np.empty(shape= alt[CNs[0]].shape[0], dtype=int)
        for id, snv_indices in tree_id_to_indices.items():
            all_arrs = []

            for k in range(dcfs.shape[1]):


    
                post_arr =[]
                for sample, cn in enumerate(CNs):
                        dcf = dcfs[sample,k]
                        a_vec = alt[cn][snv_indices]
                        d_vec=total[cn][snv_indices]
                 

                 
                   
        
                        post_arr.append(self.id_to_tree[id].vectorized_posterior(dcf, a_vec,d_vec, cn))
                sample_arr = np.vstack(post_arr).sum(axis=0)

                all_arrs.append(sample_arr)
            cluster_by_snv_index =np.vstack(all_arrs)
            cluster_choice= cluster_by_snv_index.argmax(axis=0)
            like += cluster_by_snv_index.max(axis=0).sum()
            for c,s in zip(cluster_choice, snv_indices):
                 clusters[s] = c
        for i in range(clusters.shape[0]):
             assignments.append((clusters[i], self.snv_index_to_tree_id[i]))
             
        return assignments, like
        

                  
                





        # for i in range(self.k):
        #    for tree in T_SNVs:

            #compute the probability of the given cluster dcf given tree t and the read counts
            #compute the likelihood
            
     


    def optimize_cluster_centers(self, assignments, alt, total):
        
        new_dcfs = []
        CNs = list(alt.keys())
        #optimize the cluster dcf for each sample and cluster
        for c in CNs:
             sample_dcf =[]
             for k in range(self.k):
                tree_id_to_indices= unpack_assignments(k, assignments)
                trees =[self.id_to_tree[i] for i in tree_id_to_indices]  
                if len(trees) >0:
                    new_dcf= minimize_scalar(scalar_obj, args=(tree_id_to_indices, trees, alt, total, c), method='bounded', bounds=[0,1]).x
                    
                else:
                    new_dcf = self.rng.random()    
                sample_dcf.append(new_dcf)
                    
             new_dcfs.append(sample_dcf)
        
        
        return np.array(new_dcfs)
    
    # def optimize_cluster_centers(self, snv_clusters, alt, total):
        
    #     new_dcfs = []
    #     CNs = list(alt.keys())
    #     #optimize the cluster dcf for each sample and cluster
    #     for c in CNs:
    #          sample_dcf =[]
    #          for k in range(self.k):
    #             snv_indices = np.where(snv_clusters==k)[0]
    #             tree_id_to_indices= {}
    #             for id, all_indices in self.tree_id_to_indices.items():
    #                  inters = np.intersect1d(snv_indices, all_indices)
    #                  if len(inters) > 0:
    #                       tree_id_to_indices[id] = inters
                

    #             trees =[self.id_to_tree[i] for i in tree_id_to_indices]  
    #             if len(trees) >0:
    #                 new_dcf= minimize_scalar(scalar_obj, args=(tree_id_to_indices, trees, alt, total, c), method='bounded', bounds=[0,1]).x
                    
    #             else:
    #                 new_dcf = self.rng.random()    
    #             sample_dcf.append(new_dcf)
                    
    #          new_dcfs.append(sample_dcf)
        
        
    #     return np.array(new_dcfs)
    def compute_likelihood(self,dcfs, assignments, alt, total):
   
        CNs = list(alt.keys())
        likelihood = 0
        #optimize the cluster dcf for each sample and cluster
        for j,c in enumerate(CNs):
             for k in range(self.k):
                tree_id_to_indices= unpack_assignments(k, assignments)
                trees =[self.id_to_tree[i] for i in tree_id_to_indices]  
                if len(trees) >0:
                    likelihood += -1*scalar_obj(dcfs[j,k],tree_id_to_indices, trees, alt, total, c)
                    
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
       
        for j in range(self.max_iterations):
            OMEGA = {}
            ALPHA = {}
            CNA_tree ={}
            for ell in self.segments:
                snvs = self.data.seg_to_snvs[ell]
                seg_like = np.Inf
                for S in S[ell]:
                    scriptT = clonelib.get_genotype_trees(S)
                    T_SNVs = []
                    for snv_edges in scriptT:

                        T = nx.DiGraph(snv_edges)
                        #recode the tree so the nodes are labeled by integers
                        node_mapping ={u: i for i,u in enumerate(T)}
                        T= nx.relabel_nodes(T, node_mapping)

                        T_SNVs.append(GenotypeTree(T, node_mapping))
                    
                    #find an SNV tree assignment and DCF cluster assignment
                        
                    '''
                    Optimize cluster assignments needs to be updated.
                    You will need to compute the posterior probability of each dcf q and T pair in T_SNVs for each SNV
                    Set omega equal to the tree with max posterior prob
                    Set alpha equal cluster id q with max posterior prob
                    the likelihood is the sum of the  log posterior probabilities for all optimal assignments
                    '''   
                    omega, alpha, likelihood  = self.optimize_cluster_assignments(T_SNVs, dcfs, snvs)
        
                    if likelihood > seg_like:
                        CNA_tree[ell] = S 
                        OMEGA[ell] = omega 
                        ALPHA[ell] = alpha 
                        seg_like = likelihood
                
                                    
        
            
            '''
            Optimize cluster centers need to be updated to account the assignments (omega, alpha) 
            being a dictionary of dictionaries
            


            '''  
            old_dcfs = dcfs.copy()
            dcfs, new_likelihood = self.optimize_cluster_centers(dcfs, OMEGA, ALPHA)
            dcfs[dcfs > 0.99] =1.0
            dcfs[dcfs < 1e-3] =0

            '''
            The likelihood computation needs to be updated as well 
            '''
            new_likelihood = self.compute_likelihood(dcfs,OMEGA, ALPHA)
            #check for covergence:
            diff = new_likelihood -prev_likelihood
            if self.verbose:
                print(f"Previous likelihood: {prev_likelihood} New likelihood: {new_likelihood} Diff: {diff}")
            if diff < 0 or abs(diff) <  TOLERANCE:
                 dcfs = old_dcfs 
                 break
            else:
                 prev_likelihood = new_likelihood
            

            return new_likelihood, dcfs, CNA_tree, ALPHA, OMEGA,
            
                            

    def run(self, data, k_vals= [5,6]):
       
        results = {}
        for k in k_vals:
            for i in self.nrestarts:
               results[k,i] = self.decifer(data, k)


        #add model selection 
               
        return #best 
        



if __name__ == "__main__":
     
     #load in pickled data object 
     from data import Data, load_from_pickle
     dat = load_from_pickle("pth/to_pickled_object.pkl")
     
     #initialize object 
     dec = DCF_Clustering(nrestarts=50,seed=21)
     all_results = dec.run(dat, k_vals=[i+1 for i in range(8)])

     

#     def fit(self, tree_id_to_indices, snvs, alt, total):
#         self.snv_lookup= {i: snvs[i] for i in range(len(snvs))}
#         self.snv_index_to_tree_id = {}
#         self.tree_id_to_indices = tree_id_to_indices
#         for id, indices in tree_id_to_indices.items():
#              for i in indices:
#                   self.snv_index_to_tree_id[i] = id
                  
            

#         self.nsamples = len(alt)

#         for i in range(self.nrestarts):
#             print(f"Restart: {i}")
            
#             dcfs = self.init_cluster_centers()
#             best_likelihood = np.NINF
    
#             cna_tree_likes = []
#             res = []
#             for cna_tree in self.T_CNAs:
#                 cand_T_SNVs = [tree for tree in self.T_SNVs if tree.is_refinement(cna_tree)]
#                 # tree_id_to_indices
#                 tree_results= self.fit_cna_tree(cna_tree, tree_id_to_indices, cand_T_SNVs, dcfs, alt, total)
#                 cna_tree_likes.append(tree_results.likelihood)
#                 res.append(tree_results)
#             max_index = cna_tree_likes.index(max(cna_tree_likes))

#         if best_likelihood < max(cna_tree_likes):
#             best_likelihood = max(cna_tree_likes)
#             best_res = res[max_index]
#         return best_res
            
                

            
# @dataclass 
# class DCF_Data:
#     likelihood: float 
#     assignments: list
#     T_CNA: GenotypeTree 
#     DCFs: np.array



#     def get_clusters(self):
#          return np.array([k for k, _ in self.assignments],dtype=int)
    
#     def snv_to_tree(self, snvs, T_SNVs):
#          mydict = {}
#          index = 0
#          for _, tree_id in self.assignments:
#               for t in T_SNVs:
#                    if t.id == tree_id:
#                         mydict[snvs[index]] = t
#                         index +=1
#                         break
#          return mydict

        
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



