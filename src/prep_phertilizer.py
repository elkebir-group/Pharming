# Created by: L.L. Weber
# Created on: 2024-03-13 18:02:44

from data import Data
from utils  import load_pickled_object
import argparse
import umap
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns

def umap_projection(x, y):
    cn_profiles = np.hstack([x,y])
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(cn_profiles)
    return embedding
    

def main(args):
    dat = load_pickled_object(args.data)
    embedding = umap_projection(dat.copy_x, dat.copy_y)
    # plt.scatter(
    # embedding[:, 0],
    # embedding[:, 1])
    # plt.gca().set_aspect('equal', 'datalim')
    # plt.title('UMAP projection of Allele-specific copy numbers', fontsize=24)

    if args.counts is not None:
        with open(args.counts, "w+") as file:
            for i in range(dat.N):
                for j in range(dat.m):
                    if dat.total[i,j] > 0:
                        #|chr | snv | cell |variant base | variant_reads  total_reads |
                        file.write(f"{dat.snv_to_seg[j]}\t{j}\t{i}\tA\t{dat.var[i,j]}\t{dat.total[i,j]}\n")
    

    if args.umap is not None:
        with open(args.umap, "w+") as file:
            file.write("cell,bin1,bin2\n")
            for i in range(dat.N):
                file.write(f"{i},{embedding[i,0]},{embedding[i,1]}\n")





        
    





if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-f", "--counts", required=False,
                        help="output file of read count input for phertilizer")
    parser.add_argument("-u", "--umap", required=False,
                        help="output file of the umap coordinates to give phertilizer as input")
  


    
    args = parser.parse_args()


    # instance = "s11_m5000_k25_l7"
    # folder = "n1000_c0.05_e0" 
    # pth = f"simulation_study/input"

    # args = parser.parse_args([

    #     "-d", f"{pth}/{instance}/{folder}/data.pkl",
    #     "-f", "test/phert_counts.tsv",
    #     "-u", "test/phert_umap.csv"


    # ])

    main(args)


