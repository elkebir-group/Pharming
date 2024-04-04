from genotype_tree import GenotypeTree
import numpy as np
from dataclasses import dataclass
from scipy.optimize import minimize_scalar
import clonelib
import networkx as nx
import math
import argparse
from data import Data
from scipy.optimize import linear_sum_assignment
import csv
import pickle
import pandas as pd

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

def scalar_obj_new(dcf, tree_assignments, snvs_in_cluster, alt, total, deciferObject):
    obj = 0
    # for each SNV:
    # get alt, cn, total

    for snv in snvs_in_cluster:
        seg = deciferObject.data.snv_to_seg[snv]
        cn_prob = deciferObject.data.cn_proportions(seg)
        if (1,1) not in cn_prob:
            cn_prob[(1,1)] = 0
        tree = tree_assignments[seg][snv]
        alt_summed = sum(alt)
        total_summed = sum(total)
        obj += tree.posterior_dcf(dcf, alt_summed[snv], total_summed[snv], cn_prob)
        if (obj <= EPSILON):
            b = 5

    return -obj

class DCF_Clustering:
    def __init__(self, nrestarts=25, seed=1026, verbose=False, cna_restriction=1):
        self.nrestarts =nrestarts 
   
        self.rng = np.random.default_rng(seed)

        self.verbose=verbose

        self.cna_restriction = cna_restriction


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
        #return [0.056, 0.146, 0.179, 0.996, 0.617, 0.138, 0.382]


    def optimize_cluster_and_tree_assignments(self, tree_id_to_indices, dcfs, alt, total, ell):
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
        for snv in snvs:
            likelihood[snv] = -np.inf

        for snv in snvs:
            for cluster in range(len(dcfs)):
                for tree in tree_id_to_indices: #FIX: allow SNVs to be assigned separately
                    #prob_arr = []
                    #for a, d in zip(a_vec, d_vec):
                    prob = tree.posterior_dcf(dcfs[cluster], a_vec[snv], d_vec[snv], cn_prob)
                    #prob_arr.append(prob)
                    if prob > likelihood[snv]:
                        best_cluster[snv] = cluster
                        best_tree[snv] = tree
                        likelihood[snv] = prob

        return best_cluster, best_tree, likelihood




    def optimize_cluster_centers(self, dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, alt, total): #TO DO: Check if this should be updated?
        new_dcfs = []
        for k in range(len(dcfs)): #looping over clusters
            snvs_in_cluster = []
            tree_in_cluster = []
            #find segments of clusters

            for seg, cluster in CLUSTER_ASSIGNMENTS.items():
                snvs = self.data.seg_to_snvs[seg]
                for snv in snvs:
                    if cluster[snv] == k:
                        snvs_in_cluster.append(snv)
            if len(snvs_in_cluster) > 0:
                obj = scalar_obj_new(dcfs[k], TREE_ASSIGNMENTS, snvs_in_cluster, alt, total, self)
                new_dcf = minimize_scalar(scalar_obj_new, args=(TREE_ASSIGNMENTS, snvs_in_cluster, alt, total, self), method='bounded', bounds=[0,1]).x
            else:
                #new_dcf = self.rng.random()
                new_dcf = dcfs[k]
            new_dcfs.append(new_dcf)
        return np.array(new_dcfs)


    #     return np.array(new_dcfs)
    def compute_likelihood(self, dcfs, CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS, alt, total):

        #scalar_obj_new(dcf, tree_assignments, segs_in_cluster, alt, total):




        likelihood = 0
        #optimize the cluster dcf for each sample and cluster
        for k in range(len(dcfs)):
            snvs_in_cluster = []
            for seg, cluster in CLUSTER_ASSIGNMENTS.items():
                snvs = self.data.seg_to_snvs[seg]
                for snv in snvs:
                    if cluster[snv] == k:
                        snvs_in_cluster.append(snv)
            if len(snvs_in_cluster) > 0:
                likelihood += scalar_obj_new(dcfs[k], TREE_ASSIGNMENTS, snvs_in_cluster, alt, total, self)
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
             if (1,1) not in cn_props[ell]:
                 cn_props[ell][(1,1)] = 0

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
                
                if self.cna_restriction == 1:

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
                    
                else:
                    T_SNVs = []
                    for s in S[ell]:
                        scriptT = clonelib.get_genotype_trees(s)
                        T_SNVs.append(GenotypeTree(edges= edge_list, id=i) for i, edge_list in enumerate(scriptT))
                        
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
       
        #results = {}
        #for k in k_vals:
        #    for i in range(self.nrestarts):
        #       results[k,i] = self.decifer(data, k)
        for k in k_vals:
            best_result = []
            best_likeli = np.inf
            for i in range(self.nrestarts):
                results = self.decifer(data, k)
                print(results)
                print("****")
                likeli = results[0]
                if likeli < best_likeli:
                    best_result = results
                    best_likeli = likeli


        #add model selection 
               
        return best_result#best
        

def main (data_path, ground_truth_path, output_path, accuracy_path, num_restarts, cna_restriction):
    print(data_path)
    data = pd.read_pickle(data_path)
    ground_truth = read_ground_truth_text_file(ground_truth_path)
    k = [len(ground_truth)]
    dec = DCF_Clustering(nrestarts=num_restarts, seed=21, verbose=True, cna_restriction=cna_restriction)
    all_results = dec.run(data, k_vals=k)
    dcfs = all_results[1]
    mean_difference = compute_mean_difference(ground_truth, dcfs)/k

    with open(accuracy_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([mean_difference])

    with open(output_path, 'wb') as file:
        pickle.dump(all_results, file)


def compute_mean_difference(ground_truth, dcfs):

    differences = np.abs(np.subtract.outer(ground_truth, dcfs))
    row_ind, col_ind = linear_sum_assignment(differences)
    selected_differences = differences[row_ind, col_ind]
    mean_difference = np.mean(selected_differences)

    return mean_difference


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
    parser.add_argument("pickle_path", type=str, help="Path to the pickled object")
    parser.add_argument("ground_truth", type=str, help="Path to the ground truth data")
    parser.add_argument("output_path", type=str, help="Path to the output data data")
    parser.add_argument("accuracy_path", type=str, help="Path to the accuracy data")
    parser.add_argument("num_restarts", type=int, help="Number of restarts")
    parser.add_argument("restrict_CNA_trees", type=int, help="Whether or not to enforce CNA Constraint")

    # Parse command line arguments
    args = parser.parse_args()

    # Call the main function with provided arguments
    main(args.pickle_path, args.ground_truth, args.output_path, args.accuracy_path, args.num_restarts, args.restrict_CNA_trees)


     #load in pickled data object 
     #from data import Data, load_from_pickle
     #dat = load_from_pickle("/Users/annahart/CLionProjects/Pharming/s11_m5000_k25_l7/n1000_c0.25_e0/data.pkl")


     #load in ground-truth DCFs

     #initialize object 
     #dec = DCF_Clustering(nrestarts=50,seed=21,verbose=True)
     #all_results = dec.run(dat, k_vals=[i+1 for i in range(8)])
     #results = dec.run(dat, k_vals=[7])
     #print(results)

     #store file of results

     #write to accuracy file

     

