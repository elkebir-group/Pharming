import numpy as np
import pandas as pd
import argparse
from sklearn.cluster import OPTICS
import matplotlib.pyplot as plt
import seaborn as sns
from prep_phertilizer import umap_projection

"""
First conduct a umap projection of the data. Then cluster the cells using OPTICS and 
the umap embedding. For each cluster, determine the genotype of the cells in the cluster.

"""


def cluster_cells(X, min_samples=10, eps=3):
    clustering =OPTICS(eps=eps, min_samples=min_samples).fit(X)
    highest_label = np.max(clustering.labels_)
    highest_label += 1
    labs = clustering.labels_
    # labs[labs==-1] = highest_label
    

    return  labs


def genotype_cluster(cells, data, threshold = 0.02):
    #cells are the indexes of the cells and not the labels
    var_counts = data.var[cells,  ].sum(axis=0)
    total_counts = data.total[cells,  ].sum(axis=0)
    vaf = var_counts/total_counts

 
   
    geno = np.zeros(data.m, dtype=int)
    geno[vaf >= threshold] = 1
  
    return geno

  



def write_scite_reduced(fname, reduced_matrix):
    transpose_mat = reduced_matrix.T
    np.savetxt(fname, transpose_mat.astype(int), fmt="%d")
            
def write_scite_genes(fname, m):
    with open(fname, 'w+') as file:
        for i in range(m):
            file.write(f"SNV{i}\n")

def cluster_muts(full_mat, unique_mat):
    column_vecs = full_mat.T.tolist()

    mclusters = []
    for i,c in enumerate(column_vecs):
        for j in range(unique_mat.shape[1]):
            seq = unique_mat[:,j].tolist()
            if c == seq:
                mclusters.append(str(j))
    return(mclusters)


def write_clusters(mlabels, clabels, fname):
    # print(clabels.head())
    clabels = clabels.values
    clab_string =",".join(clabels.astype(str)) + "\n"
    mlab_string = ",".join(mlabels) + "\n"
    with open(fname, 'w') as out:
        out.write(clab_string)
        out.write(mlab_string)
        
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Baseline comparison method for ultra-low coverage inference')
    parser.add_argument("-d", "--data", required=True,
                        help="input data object")
    parser.add_argument("-t" ,"--geno_threshold", type=float, default=0.02,
                        help="percent of cells with a variant read ")
    parser.add_argument("--min_frac", type=float, default=0.05,
                        help="parameter to decide min sample size for OPTIC clustering")
    parser.add_argument("-s", "--seed", type=int,default=1026,
                        help="random number seed")             
    parser.add_argument("--scite", type=str, 
                        help="name of SCITE input file")
    parser.add_argument("--genes", type=str,
                         help= "filename for SCITE geneNames")
    parser.add_argument("--umap", type=str,
                         help=  "filename for png of umap embedding")
    parser.add_argument("-c", "--cell-clusters", type=str,
                         help=  "mapping of cell labels to clusters")
    parser.add_argument("-m", "--mut-clusters", type=str,
                         help=  "mapping of mut labels to clusters")
            
    args = parser.parse_args()


 
    # pth = "simulation_study/sims4"
    # inpath = f"{pth}/s10_m5000_k25_l5_d2/n1000_c0.25_e0"
    # bins = "reads_per_bin_relabeled.csv"
    # outpath = f"test"

    # args=parser.parse_args([
    #     "-d",  f"{inpath}/data.pkl",
    #     "--umap", f'{outpath}/umap.png',
    #     "-s", "10",
    #     "-t", "0.05",
    #     "--cell-clusters", f"{outpath}/cell_clusters.csv",
    #     "--mut-clusters", f"{outpath}/mut_clusters.csv",
    #     "--scite", f"{outpath}/scite.in",
    #     "--genes", f"{outpath}/scite_genes.txt"
    # ])

    dat = pd.read_pickle(args.data)

    embedding = umap_projection(dat.copy_x, dat.copy_y) 
    min_samples = int(args.min_frac*dat.N)
    labels = cluster_cells(embedding, min_samples=min_samples)
    if args.umap is not None:
        plt.scatter(embedding[:,0], 
                    embedding[:,1],
                    c = [sns.color_palette()[x] for x in labels])
        plt.title('UMAP embedding')
        plt.show() #for control
        plt.savefig(args.umap)


    labs = np.sort(np.unique(labels))
    num_clusters= len(labs)
    
    genotypes = {}
    ncells = {}
    
  
    geno_list = []
 


    for l in labs:
        cells_in_cluster = np.where(labels == l)[0]            
        genotypes[l] = genotype_cluster(cells_in_cluster, dat, threshold=args.geno_threshold)
        geno_list.append(genotypes[l])
        ncells[l] = cells_in_cluster.shape[0]
    matrix = np.vstack(geno_list)
    reduced_matrix =np.unique(matrix, axis=1)
 
    mut_clusters = cluster_muts(matrix, reduced_matrix)
    mut_lookup = dat.mut_lookup
    mut_labels = mut_lookup.values

    mut_clust = pd.Series(data=mut_clusters)
    mut_df = mut_clust.reset_index().rename(columns={'index': 'index', 0: 'cluster'})
    cell_clust = pd.Series(data=labels)
    cell_clust_df = cell_clust.reset_index().rename(columns={'index': 'index', 0: 'cluster'})

    if args.cell_clusters is not None:
        cell_clust_df.to_csv(args.cell_clusters, index=False)
    if args.mut_clusters is not None:
        mut_df.to_csv(args.mut_clusters, index= False)


    if args.scite is not None:
        write_scite_reduced(args.scite, reduced_matrix)
        write_scite_genes(args.genes, reduced_matrix.shape[1])
    