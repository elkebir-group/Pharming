from genotype_tree import GenotypeTree
import numpy as np
from scipy.optimize import minimize_scalar, minimize
import clonelib
import networkx as nx
from collections import defaultdict 
import argparse
from scipy.optimize import linear_sum_assignment
import csv
import pickle
import pandas as pd
from utils import timeit_decorator
import  multiprocessing
from decifer_post_process import post_process
TOLERANCE = 2
EPSILON = -1e40
SEQERROR = 1e-40






def scalar_obj( dcf, id_to_index, trees, alt_vec, total_vec,cn):
        obj = 0
        for id, tree in zip(id_to_index, trees): #could loop over snv/tree matchings and do non-vectoirzed posterior
            alt = alt_vec[cn][id_to_index[id]]
            total =total_vec[cn][id_to_index[id]]
            obj += tree.vectorized_posterior(dcf, alt, total,cn).sum()
        return -1*obj

def objective_function(x, clust_group_map, q, data):
    return scalar_obj_val(x[0], clust_group_map, q, data)

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

def scalar_obj_val(dcf, clust_group_map, q, data):
    obj = 0
    for ell in clust_group_map:
        cn_prob = data.cn_proportions(ell)
        if (1,1) not in cn_prob:
            cn_prob[(1,1)] = 0
        if q in clust_group_map[ell]:
            for tree, snvs in clust_group_map[ell][q]:
                if len(snvs) > 0:
                    alt = data.var_marg(snvs)
                    total = data.total_marg(snvs)
                    obj +=  -1*tree.vectorized_posterior(dcf, alt, total, cn_prob).sum()
    return obj
# def scalar_obj_new(dcf, tree_assignments, snvs_in_cluster, alt, total, deciferObject):
#     obj = 0
#     # for each SNV:
#     # get alt, cn, total

#     for snv in snvs_in_cluster:
#         seg = deciferObject.data.snv_to_seg[snv]
#      
#         tree = tree_assignments[seg][snv]
#         alt_summed = sum(alt)
#         total_summed = sum(total)
#         obj += tree.posterior_dcf(dcf, alt_summed[snv], total_summed[snv], cn_prob)
#         if (obj <= EPSILON):
#             b = 5

#     return -obj

class DCF_Clustering:
    def __init__(self, nrestarts=25,  seed=1026, verbose=False, cna_restriction=True, rng=None, iterations=100 ):
        self.nrestarts =nrestarts 
        if rng is None:
            self.rng = np.random.default_rng(seed)
        else:
            self.rng = rng 

        self.iterations = iterations
        self.verbose=verbose

        self.cna_restriction = cna_restriction


    
    

    # def init_cluster_centers(self,k, include_truncal=True):
    #     ''' 
    #     randomly initialize cluster centers 
    #     '''
    #     init = self.rng.uniform(low=0.0, high=1.0, size=k)
    #     if include_truncal:
    #         init[np.argmax(init)] = 0.98
    #     return init
  
    def init_cluster_centers(self, k):
        block_size = 1.0 / k
        # Generate random values uniformly in [0, 1] for each block
        block_uniform_values = self.rng.uniform(size=k)
        # Calculate the starting point of each block
        block_starts = np.arange(k) * block_size
        # Calculate the random sample within each block
        samples = block_starts + block_uniform_values * block_size
        return samples


        #return [0.056, 0.146, 0.179, 0.996, 0.617, 0.138, 0.382]


    def optimize_cluster_and_tree_assignments(self, snv_trees, dcfs, alt, total, ell):
        snvs = self.data.seg_to_snvs[ell]
        cn_prob = self.data.cn_proportions(ell)
        if (1,1) not in cn_prob:
            cn_prob[(1,1)] = 0

        # get a_vec and d_vec for the snvs
        a_vec = alt.sum(axis=0)
        d_vec = total.sum(axis=0)
        best_cluster = {}
        best_tree = {}
        likelihood = {}
        # for snv in snvs:
            # likelihood[snv] = -np.inf

        all_clust_likes= []
        all_clust_tree_assign = []
        for cluster in range(len(dcfs)):
            # for tree in tree_id_to_indices: #FIX: allow SNVs to be assigned separately
            #     #prob_arr = []
            #     #for a, d in zip(a_vec, d_vec):
            #     if isinstance(tree, list):
            #         for t in tree:
            #             prob = t.posterior_dcf(dcfs[cluster], a_vec[snv], d_vec[snv], cn_prob)
            #             #prob_arr.append(prob)
            #             if prob > likelihood[snv]:
            #                 best_cluster[snv] = cluster
            #                 best_tree[snv] = t
            # #                 likelihood[snv] = prob
            #     else:
            tree_probs = np.vstack([tree.vectorized_posterior(
                            dcfs[cluster], a_vec[snvs], d_vec[snvs], cn_prob) for tree in snv_trees])
            clust_like = np.max(tree_probs, axis=0) 
            clust_tree_assign = np.argmax(tree_probs, axis=0)
            all_clust_likes.append(clust_like)
            all_clust_tree_assign.append(clust_tree_assign)

        clikes = np.vstack(all_clust_likes)
        opt_clust = np.argmax(clikes, axis=0)
        opt_like = np.max(clikes, axis=0)

        clust_to_snvs = defaultdict(list)
        for i, m in enumerate(snvs):
            q = opt_clust[i]
            best_cluster[m] = q
            clust_to_snvs[q].append(m)
            tree_index = all_clust_tree_assign[q][i]
            best_tree[m] = snv_trees[tree_index]
            likelihood[m] = opt_like[i]

        clust_group = defaultdict(list)
        for q, clust_snvs in clust_to_snvs.items():
            selected_tree = all_clust_tree_assign[q]
            indices = [i for i in range(len(snvs)) if snvs[i] in clust_snvs]  
            snv_to_tree_index = selected_tree[indices] 

            for  t in set(snv_to_tree_index):
        
                tree_indices = np.where(selected_tree==t)[0]
                best_tree_indices = set(indices).intersection(set(tree_indices))
                if len(best_tree_indices) > 0:
               
                    clust_group[q].append((snv_trees[t], [snvs[j] for j in best_tree_indices]))

                        # # prob_arr.append(prob)
                        # if prob > likelihood[snv]:
                        #     best_cluster[snv] = cluster
                        #     best_tree[snv] = tree
                        #     likelihood[snv] = prob
        all_snvs = []
        for q, mylist in clust_group.items():
            for tree, snv_list in mylist:
                all_snvs += snv_list
        assert len(all_snvs) == len(snvs)
        return best_cluster, best_tree, likelihood, clust_group




    def optimize_cluster_centers(self, dcfs, clust_group_map): 
        new_dcfs = np.zeros(len(dcfs))

        for q in range(len(dcfs)): #looping over clusters
            # new_dcfs[q] = minimize_scalar(scalar_obj_val, args=(clust_group_map, q, self.data), 
            #                               method='bounded', bounds=[0,1]).x
            result = minimize(objective_function, x0=[dcfs[q]], args=(clust_group_map, q, self.data), bounds=[(0.025,1)])
            if result ==0:
                new_dcfs[q] = self.rng.uniform()
            else:
                new_dcfs[q]= result.x[0]

           
        return new_dcfs




    def compute_likelihood(self, dcfs, clust_group_mapping):
        likelihood = 0
        #scalar_obj_new(dcf, tree_assignments, segs_in_cluster, alt, total):
        for q, dcf in enumerate(dcfs):
         likelihood += scalar_obj_val(dcf, clust_group_mapping, q, self.data)

        return likelihood

    # @timeit_decorator
    def decifer(self, data, dcfs):
        '''
        See data.py for the data object 
        k = # of SNV clusters 
        '''
        init_dcfs =dcfs.copy()
        prev_likelihood = np.NINF
        self.k = len(dcfs)
        # dcfs = self.init_cluster_centers()
        self.data = data 
        self.segments = list(self.data.seg_to_snvs.keys())
        #enumerat valid CNA trees for each segment
        S = {}
        cn_props ={}
        for ell in self.segments:
             cn_props[ell] = self.data.cn_proportions(ell)
             if (1,1) not in cn_props[ell]:
                 cn_props[ell][(1,1)] = 0

             #find all copy number states (x,y) and proportions in segment ell
             states = set(cn_props[ell].keys())
             #enumerate CNA trees 
             S[ell] = clonelib.get_cna_trees(states, 1, 1)
             
        #S = {\ell: [CNA trees]}
       
 
        for j in range(self.iterations):
            CLUSTER_ASSIGNMENTS = {}
            TREE_ASSIGNMENTS = {}
            CNA_tree ={}
            CLUST_GROUP_MAPPING= {}
            for ell in self.segments: #TO DO: Swap these segments and CNA tree assignment

                #skip segments with a lot of copy number states but a small number of SNVs
                thresh_cn_prop= self.data.thresholded_cn_prop(ell, 0.05)
                if len(thresh_cn_prop) > 4 and self.data.num_snvs(ell) < 100:
                    continue
          
                seg_like = np.NINF
                
                if self.cna_restriction:

                    for s in S[ell]:
                        #here, get optimal cluster assignments, then save the optimal CNA tree + assignments of SNVs to trees + assignments of SNVs to clusters
                        #after doing this for all the segments, then update the clusters (using assignments & dcf values)
                   
                            
                        if len(s) ==0:
                            states = list(cn_props[ell].keys())
                            assert len(states) ==1
                            cn = states[0]
                            scriptT=[[((*cn, 0, 0), (*cn, 1,0))]]
                        else:
                            scriptT = clonelib.get_genotype_trees(s)
                        # print(len(scriptT))
                        T_SNVs = [GenotypeTree(edges= edge_list, id=i) for i, edge_list in enumerate(scriptT) ]
         
                        #find an SNV tree assignment and DCF cluster assignment
                            
                        '''
                        TO DO: Optimize cluster assignments needs to be updated.
                        You will need to compute the posterior probability of each dcf q and T pair in T_SNVs for each SNV
                        Set omega equal to the tree with max posterior prob
                        Set alpha equal cluster id q with max posterior prob
                        the likelihood is the sum of the  log posterior probabilities for all optimal assignments
                        '''   
                        cluster, tree, likelihood, clust_group  = self.optimize_cluster_and_tree_assignments(T_SNVs, dcfs, self.data.var, self.data.total, ell)
                        new_seg_likelihood = sum(likelihood[m] for m in likelihood)
                        # new_seg_likelihood = 0
                        # for key, val in likelihood.items():
                        #     new_seg_likelihood += val

                        if new_seg_likelihood > seg_like:
                            CNA_tree[ell] = s
                            CLUSTER_ASSIGNMENTS[ell] = cluster
                            TREE_ASSIGNMENTS[ell] = tree
                            seg_like = new_seg_likelihood
                            CLUST_GROUP_MAPPING[ell] = clust_group
                    
                else:
                    T_SNVs = []
                    for s in S[ell]:
                        scriptT = clonelib.get_genotype_trees(s)
                        T_SNVs.append([GenotypeTree(edges= edge_list, id=i) for i, edge_list in enumerate(scriptT)])
                        
                    cluster, tree, likelihood, clust_group  = self.optimize_cluster_and_tree_assignments(T_SNVs, dcfs, self.data.var, self.data.total, ell)
                    # new_seg_likelihood = 0
                    new_seg_likelihood = sum(likelihood[m] for m in likelihood)
           
                    if new_seg_likelihood > seg_like:
                        CNA_tree[ell] = s
                        CLUSTER_ASSIGNMENTS[ell] = cluster
                        TREE_ASSIGNMENTS[ell] = tree
                        CLUST_GROUP_MAPPING[ell] = clust_group
                        seg_like = new_seg_likelihood
                    
                                    
        
            
            '''
            TO DO: Optimize cluster centers need to be updated to account the assignments (omega, alpha) 
            being a dictionary of dictionaries
            


            '''  
            # old_dcfs = dcfs.copy()
            # dcfs = self.optimize_cluster_centers(dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, self.data.var, self.data.total)
            # test_like = self.compute_likelihood(dcfs, CLUST_GROUP_MAPPING)
            dcfs = self.optimize_cluster_centers(dcfs, CLUST_GROUP_MAPPING)
      
            dcfs[dcfs > 0.999] =1.0
            dcfs[dcfs < 1e-3] =0
     
            '''
            TO DO: The likelihood computation needs to be updated as well 
            '''
            new_likelihood = self.compute_likelihood(dcfs,CLUST_GROUP_MAPPING)
            #check for covergence:
            diff = new_likelihood -prev_likelihood
            if self.verbose:
                print(f"{j}: Previous likelihood: {prev_likelihood} New likelihood: {new_likelihood} Diff: {diff}")
            if abs(diff) <  TOLERANCE:
                 break
            else:
                 prev_likelihood = new_likelihood

            
        if self.verbose:
            print(f"Restart complete with likelihood: {new_likelihood}")
            print(f"Initial DCFs: {init_dcfs}")
            print(f"Inferred DCFs: {dcfs}")
 
            clust_assign = defaultdict(list)
            for ell in CLUSTER_ASSIGNMENTS:
                for j, q in CLUSTER_ASSIGNMENTS[ell].items():
                    clust_assign[q].append(j)
            for q in clust_assign:
                print(f"{q}: {dcfs[q]}: {len(clust_assign[q])} SNVs")

        return new_likelihood, dcfs, CNA_tree, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS

    @staticmethod         
    def elbow_criteria_model_selection(objs, mink, maxk, ubleft=0.06):
        if mink == maxk:
            return mink, objs, {}
   
        for k in range(mink + 1, maxk + 1):
            if objs[k - 1] < objs[k]:
    
                objs[k] = objs[k - 1]
            

        chk = (lambda v: v if v != 0.0 else 0.01)
        left = (lambda k: min((objs[k - 1] - objs[k]) / abs(chk(objs[k - 1])), ubleft) if k > mink else ubleft)
        right = (lambda k: (objs[k] - objs[k + 1]) / abs(chk(objs[k])))
        # left_vals = {k: left(k) for k in range(mink, maxk) }
        # right_vals = {k: right(k) for k in range(mink, maxk) }
        elbow = {k: left(k) - right(k) for k in range(mink, maxk)}
   
    
        selected = max(range(mink, maxk), key=(lambda k: elbow[k]))

        return selected, objs, elbow
                            
    @timeit_decorator
    def run(self, data, k_vals, cores=1):
       
        #results = {}
        #for k in k_vals:
        #    for i in range(self.nrestarts):
        #       results[k,i] = self.decifer(data, k)
        results_by_k ={}
        likelihoods  = {}
        for k in k_vals:
            best_result = []
            best_likeli = np.inf
            all_results = []
            if cores <= 1:
                for i in range(self.nrestarts):
                    init_dcfs =self.init_cluster_centers(k)
                    results = self.decifer(data, init_dcfs)
                    all_results.append(results)
                    if self.verbose:
                        print(f"Restart {i}: {results[0]} ")
                        print(f"DCFS: {results[1]}")
            
            else:
                args = [(data, self.init_cluster_centers(k)) for _ in range(self.nrestarts)]
                with multiprocessing.Pool(cores) as pool:
                    all_results = pool.starmap(self.decifer, args)
                    # print(results)
                    # print("****")
        
            for results in all_results:
                likeli = results[0]
                if likeli < best_likeli:
                    best_result = results
                    best_likeli = likeli
            likelihoods[k] = best_likeli
            results_by_k[k] = best_result
        
        print("Starting model selection using elbow criteria..")
        selected, updated_obj, elbow =   self.elbow_criteria_model_selection(likelihoods, min(k_vals), max(k_vals))
        print(updated_obj)
        print(elbow)

        #add model selection 
               
        return results_by_k[selected]#best
        

def main(args):
    
    data = pd.read_pickle(args.pickle_path)
    print(data)
    
    dec = DCF_Clustering(nrestarts=args.num_restarts, seed=args.seed, 
                         verbose=True, cna_restriction=args.restrict_CNA_trees)

 

    if args.ground_truth is not None:

        ground_truth = read_ground_truth_text_file(args.ground_truth)
        k_vals = [len(ground_truth)]
        res  = dec.decifer(data, np.array(ground_truth))
        gt_like = res[0]
        print(f"GT Likelihood: {gt_like}")
        print(f"DCFS: {res[1]}")
    elif args.clusters is not None:
        k_vals  = [args.clusters]
    else:
        if args.mink > args.maxk:
            maxk = args.mink
            mink = args.maxk
        else:
            mink = args.mink
            maxk = args.maxk

        #do one extra value of k to check to allow maxk to be an elbow point.
        k_vals = [k for k in range(mink, maxk+2)]

    best = dec.run(data, k_vals=k_vals, cores=args.cores)



    dcfs = best[1]
    print("Best DCFs:")
    print(dcfs)


    if args.output_path is not None:
        with open(args.output_path, 'wb') as file:
            pickle.dump(best, file)

    if args.dcfs is not None:

        with open(args.dcfs, "w+") as file:
            for d in dcfs:
                file.write(f"{d}\n")
    
    if args.post_dcfs is not None:
        clust_assign, tree_assign = best[3], best[4]
        snv_assign = defaultdict(list)
        tassign = {}
        for ell in clust_assign:
            for j,q in clust_assign[ell].items():
                snv_assign[q].append(j)
            for j, t in tree_assign[ell].items():
                tassign[j] = t
        dcf_dict = {i: d for i, d in enumerate(dcfs)}
        post_dcfs = post_process(dcf_dict, snv_assign, tassign, data)
        with open(args.post_dcfs, "w+") as file:
            for q, d in post_dcfs.items():
                file.write(f"{d}\n")
    
    



def read_ground_truth_text_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            # Remove whitespace and convert to float
            value = float(line.strip())
            data.append(value)
    return data



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform clustering analysis on data from a pickled object")
    parser.add_argument("-d","--pickle_path", type=str, help="Path to the pickled object")
    parser.add_argument("-g", "--ground_truth", type=str, help="Path to the ground truth data")
    parser.add_argument("-o", "--output_path", type=str, help="Path to the output  data")
    parser.add_argument("-r", "--num-restarts", type=int, help="Number of restarts")
    parser.add_argument("-c", "--restrict_CNA_trees", action="store_true")
    parser.add_argument("-k", "--clusters", type=int, help="number of clusters")
    parser.add_argument("--mink", type=int, default=2, help="minimum number of clusters")
    parser.add_argument("--maxk", type=int, default= 5, help="maximum number of clusters")
    parser.add_argument("-D", "--dcfs", type=str, help="output file for inferred dcfs")
    parser.add_argument("-P", "--post-dcfs", type=str, help="output file for post-processed dcfs")
    parser.add_argument("-s", "--seed", type=int, default=21,  help="output file for inferred dcfs")
    parser.add_argument("-j", "--cores", default=1, type=int, help="number of cores to use")

    # Parse command line arguments
    args = parser.parse_args()

    # gtpth = "test"
    # seed = 10
    # cov = 0.01
    # instance = f"s{seed}_m5000_k25_l5_d2"
    # folder = f"n1000_c{cov}_e0" 
    # pth = f"simulation_study/sims"



    # args = parser.parse_args([

    #     "-d", f"{pth}/{instance}/{folder}/data.pkl",
    #     # "-k", "4",
    #     "--mink", "4",
    #     "--maxk", "6",
    #     "-j", "5",
    #     # "-g", f"{pth}/{instance}/{folder}/dcfs.txt",
    #     "-s", f"{seed}",
    #     "-r", "25",
    #     "-D", f"{gtpth}/dcfs.txt",
    #     "-P", f"{gtpth}/post_dcfs.txt",
    #     "-o",  f"{gtpth}/results.pkl",
        
    #     "-c"

    # ])


    # Call the main function with provided arguments
    main(args)



     

# def silhouette_model_selection(shared, best, mink, maxk):
#     silhouette_scores = {}
#     for k in range(mink, maxk + 1):
#         mutations = shared['bmut'][best[k]]
#         X, labels = [], []
#         for mut in mutations:
#             # don't use mutations that couldn't be assigned to any cluster
#             if mut.assigned_cluster > 0:
#                 labels.append( mut.assigned_cluster )
#                 var = mut.a
#                 tot = mut.d
#                 vaf = [float(v) / t if t > 0 else 0.5 for v, t in zip(var, tot)]
#                 estC = [mut.assigned_config.v_to_cf(vaf[sam], sam, truncate = False)/PURITY[sam] for sam in range(len(mut.a))]
#                 X.append(estC)
#         silhouette_scores[k] = silhouette_score(np.array(X), np.array(labels))
#     selected = max(range(mink, maxk+1), key=(lambda k: silhouette_scores[k]))
#     return selected, silhouette_scores

