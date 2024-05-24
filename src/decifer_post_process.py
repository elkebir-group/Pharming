
import pandas as pd 
from collections import defaultdict
from genotype_tree import GenotypeTree
import ast 
import numpy as np
from scipy.optimize import minimize_scalar
from sklearn.cluster import SpectralBiclustering
import argparse

def scalar_obj_val(dcf, snv_trees, data):
    obj = 0
    for j in snv_trees:
        tree =  snv_trees[j]
        ell = data.snv_to_seg[j]
        cn_prob = data.thresholded_cn_prop(ell, 0.05)
        alt = data.var_marg([j])
        total = data.total_marg([j])
        obj +=  -1*tree.vectorized_posterior(dcf, alt, total, cn_prob).sum()
    return obj



def get_clusters(labs, vals):
        subclust = {}
        for u in [0,1]:
            indices = np.where(labs==u)[0]
            subclust[u] = [vals[j] for j in indices]
        return subclust

def post_process(old_dcfs, clust_assign, snv_trees, dat, snv_thresh=0.05, cell_thresh=0.05, vaf_thresh=0.1):
    dcfs = old_dcfs.copy()
    #keep only clusters that have at least snv_thresh proportion of SNVs assigned
    nsnvs = sum(len(clust_assign[q]) for q in clust_assign)
    clust = []
    for q in clust_assign:
        if len(clust_assign[q])/nsnvs < snv_thresh:
            print(f"Removing cluster {q} with DCF {dcfs[q]} due to insufficient number of SNVs....")
            del dcfs[q]
        else:
            clust.append(q)


    for q in clust:
        print(f"Checking cluster {q} with DCF {dcfs[q]}....")
        snvs = np.array(clust_assign[q])
        X = dat.var[:,snvs]
        snv_marg = X.sum(axis=0)
        cell_marg = X.sum(axis=1)
 
        snvs = snvs[snv_marg > 0]
        X = X[cell_marg >0, :]
        X = X[:, snv_marg > 0]
        cells = dat.cells[cell_marg > 0]
        try:
            clustering = SpectralBiclustering(n_clusters=(2,2), random_state=5, n_init=25).fit(X)
        except:
            print(f"Biclustering could not be found for cluster {q}")
            continue
        snv_clust = get_clusters(clustering.column_labels_, snvs)
        valid_size = True
        for u,csnvs in  snv_clust.items():
            if len(csnvs)/nsnvs < snv_thresh:
                print(f"{len(csnvs)} SNVs does not meet SNV threshold!")
                valid_size = False
                break
        
        cell_clust = get_clusters(clustering.row_labels_, cells)
        for u, ccells in cell_clust.items():
            if len(ccells)/dat.N < cell_thresh:
                print(f"{len(ccells)} cells does not meet cell threshold!")
                valid_size = False
                break
        
        if not valid_size:
            continue

       
   
        vaf = np.zeros((2,2))
        for u in [0,1]:
            for v in [0,1]:
               var = dat.var[np.ix_(cell_clust[u], snv_clust[v])].sum()
               total = dat.total[np.ix_(cell_clust[u], snv_clust[v])].sum()
               vaf[u,v] = var/total
        
        diag1 = np.all(np.diag(vaf) < vaf_thresh) and np.all(np.fliplr(vaf).diagonal() > vaf_thresh)
        diag2 = np.all(np.fliplr(vaf).diagonal() < vaf_thresh) and np.all(np.diag(vaf) > vaf_thresh)
        if diag1 or diag2:

            print(f"Recommend splitting cluster {q}, finding new cluster centers....")
            print(vaf)
            old_dcf = dcfs[q]
            del dcfs[q]
            for u, split_snvs in snv_clust.items():
                trees = {j: snv_trees[j] for j in split_snvs}
                new= minimize_scalar(scalar_obj_val, args=(trees, dat), 
                                          method='bounded', bounds=[0,1]).x
                # new_dcfs.append(new)
                print(f"Old DCF: {old_dcf} New DCF: {new}")
                if u == 0:
                    dcfs[q] = new 
                else:
                    new_key = max(dcfs.keys()) + 1
                    dcfs[new_key] = new

    return dcfs
   

               


      

    
            

 


    # groups = []
    # group_snvs = defaultdict(list)
    #     #find groups of SNVs incompatible 
    #     for j in clust_assign[q]:
    #         t = snv_trees[j]
    #         consist = False 
    #         if len(groups) ==0:
    #             groups.append(t)
    #             group_snvs[j] =0
    #             continue
    #         while not consist:
    #             for idx,gtree in enumerate(groups):
    #                 if gtree.is_consistent(t):
    #                     consist = True 
    #                     group_snvs[idx].append(j)
    #                     break
    #             if consist:
    #                 continue
    #             groups.append(t)
    #             group_snvs[len(groups)].append(j)
    #             consist =True

    #     if len(groups) > 1:
    #         print(f"Found cluster {q} with incompatible SNV trees")
    #         for idx in len(groups):
    #             print(f"{idx}: {len(group_snvs[idx])}")
            



def read_input(fname):
    #get cluster dcfs 
    #get snvs assignments to clusters 
    #get snv trees

 
    dec = pd.read_table(fname)
    dec= dec[['mut_index', 'cluster', 'state_tree', 'point_estimate_DCF0', 'true_cluster_DCF0']]
    dec[['clust_dcf', 'CI']] = dec['true_cluster_DCF0'].str.split(';', expand=True)
    dec['clust_dcf'] = dec['clust_dcf'].astype(float)
    dcfs = dec[['clust_dcf', 'cluster']].drop_duplicates()
    dcfs = dict(zip(dcfs['cluster'], dcfs['clust_dcf']))
    snv_assign = dict(zip(dec['mut_index'], dec['cluster']))

    snv_trees = {}
    for j, tstr in zip(dec['mut_index'],dec["state_tree"]):
        def add_allele(tup):
            e = list(tup)
            mut_copies = e[2]
            x = e[0]
            y = e[1]

            if mut_copies == 0:
                w,z = 0, 0
            elif mut_copies == x:
                w = mut_copies 
                z =0 
            else:
                z = mut_copies
                w = 0
            return (x,y,w,z)
        clean = (lambda e: tuple([add_allele(ast.literal_eval(u)) for u in  e.split("->")]))
     
        tstr_edges = tstr.split(";")
        t = [clean(e) for e in tstr_edges]
    
        snv_trees[j] = GenotypeTree(t)

    clust = defaultdict(list)
    for j,q in snv_assign.items():
        clust[q].append(j)
    
    return dcfs, clust, snv_trees



def main(args):
    dat = pd.read_pickle(args.pickle_path)
    dcfs, clust, snv_trees = read_input(args.decifer)
    print("\nPost-processing inferred DCFs...")
    new_dcfs = post_process(dcfs, clust, snv_trees,dat, args.snv_thresh, args.cell_thresh, args.vaf_thresh)
    print(new_dcfs)
    with open(args.out, "w+") as file:
        for q, d in new_dcfs.items():
            file.write(f"{d}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform clustering analysis on data from a pickled object")
    parser.add_argument("-d","--pickle_path", type=str, help="Path to the pickled object")
    parser.add_argument("-f", "--decifer", type=str, help="Path to output file of Decifer")
    parser.add_argument("-o", "--out", type=str, help="Path to write post-processed DCFs")  
    parser.add_argument( "--vaf-thresh", type=float, default=0.1, help="VAF threshold")  
    parser.add_argument( "--snv-thresh", type=float, default=0.05, help="SNV threshold")  
    parser.add_argument( "--cell-thresh", type=float, default=0.05, help="cell threshold") 
    
    # cov = 0.25
    # seed = 12
    # fname = f"/scratch/leah/Pharming/simulation_study/decifer/s{seed}_m5000_k25_l5_d2/n1000_c{cov}_e0/decifer_output.tsv" 
    # datfile = f"/scratch/leah/Pharming/simulation_study/sims/s{seed}_m5000_k25_l5_d2/n1000_c{cov}_e0/data.pkl"
    args = parser.parse_args()

    # args = parser.parse_args([
    #     "-d", datfile,
    #     "-f", fname,  
    #     "-o", "test/post_dcfs.txt"  
    #     ])
    main(args)

# gt = pd.read_pickle(gtfile)
# psi = gt.get_psi()




  