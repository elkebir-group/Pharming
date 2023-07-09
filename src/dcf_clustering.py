from genotype_tree import GenotypeTree
import numpy as np
from dataclasses import dataclass
from scipy.optimize import minimize_scalar



TOLERANCE = 1e-04
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
    def __init__(self, T_CNAs, T_SNVs, clusters= 3, nrestarts=10, rng=None, verbose=False):
        self.nrestarts =nrestarts 
        if rng is None:
            self.rng = rng 
        else:
            self.rng = np.random.default_rng(1026)
        # self.clusters= [i+1 for i in range(max_clusters)]
        
        self.k = clusters
        self.T_CNAs = T_CNAs

        self.T_SNVs = T_SNVs
        self.id_to_tree = {tree.id: tree for tree in T_SNVs}
        self.max_iterations=100
        self.combos = []
   
        for k in range(self.k):
            for t,tree in enumerate(self.T_SNVs):
                    self.combos.append((k,t))
        
        self.verbose = verbose
    

    def init_cluster_centers(self):
        ''' 
        randomly initialize cluster centers 
        '''
        return self.rng.uniform(low=0.0, high=1.0, size=(self.nsamples, self.k))


    def optimize_snv_assignments(self, T_SNVs, dcfs, alt, total):
        #assignment of snvs to T_SNV_tree
        # T_SNV_dict= {}
        CNs = list(alt.keys())
        #for each homogenous sample
        combos = []
        all_arrs = []
        for k in range(dcfs.shape[1]):
                
                for tree in T_SNVs:
                    combos.append((k, tree.id))
          

                    post_arr =[]
                    for sample, cn in enumerate(CNs):
                        a_vec = alt[cn]
                        d_vec=total[cn]
                 
                        dcf = dcfs[sample,k]
        
                        post_arr.append(tree.vectorized_posterior(dcf, a_vec,d_vec, cn))
                    sample_arr = np.vstack(post_arr).sum(axis=0)
                    all_arrs.append(sample_arr)
        cand_assign = np.vstack(all_arrs)
        likelihood = cand_assign.max(axis=0).sum()
        combo_assign = np.argmax(cand_assign, axis=0)
        #assignments is list of length snvs containing a two-tuple of cluster and tree id to which the 
        assignments = [combos[c] for c in combo_assign ]
 
        return assignments, likelihood

                






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
        
        

    def fit_cna_tree(self, T_CNA, T_SNVs, dcfs,  alt, total, ):


    
        for j in range(self.max_iterations):

            assignments, likelihood = self.optimize_snv_assignments(T_SNVs, dcfs, alt, total)
            old_dcfs = dcfs.copy()
            dcfs = self.optimize_cluster_centers(assignments, alt, total)
            new_likelihood = self.compute_likelihood(dcfs, assignments, alt, total)
            diff = new_likelihood -likelihood
            if self.verbose:
                print(f"Previous likelihood: {likelihood} New likelihood: {new_likelihood} Diff: {diff}")
            if diff < 0 or abs(diff) <  TOLERANCE:
                 dcfs = old_dcfs 
                 break
            # if likelihood > best_likelihood:
            #     best_likelihood = likelihood
        
        return DCF_Data(likelihood, assignments, T_CNA, dcfs)

    def fit(self, snvs, alt, total):
        self.snv_lookup= {i: snvs[i] for i in range(len(snvs))}

        self.nsamples = len(alt)

        for i in range(self.nrestarts):
            
            dcfs = self.init_cluster_centers()
            best_likelihood = np.NINF
    
            cna_tree_likes = []
            res = []
            for cna_tree in self.T_CNAs:
                cand_T_SNVs = [tree for tree in self.T_SNVs if tree.is_refinement(cna_tree)]
                tree_results= self.fit_cna_tree(cna_tree, cand_T_SNVs,dcfs, alt, total)
                cna_tree_likes.append(tree_results.likelihood)
                res.append(tree_results)
            max_index = cna_tree_likes.index(max(cna_tree_likes))

        if best_likelihood < max(cna_tree_likes):
            best_likelihood = max(cna_tree_likes)
            best_res = res[max_index]
        return best_res
            
                

            
@dataclass 
class DCF_Data:
    likelihood: float 
    assignments: list
    T_CNA: GenotypeTree 
    DCFs: np.array



    def get_clusters(self):
         return np.array([k for k, _ in self.assignments],dtype=int)
    
    def snv_to_tree(self, snvs, T_SNVs):
         mydict = {}
         index = 0
         for _, tree_id in self.assignments:
              for t in T_SNVs:
                   if t.id == tree_id:
                        mydict[snvs[index]] = t
                        index +=1
                        break
         return mydict

        
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



