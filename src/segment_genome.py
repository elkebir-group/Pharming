import numpy as np
import pandas as pd
import argparse

class Segment:
    def __init__(self, tol=0.0, pseudo_counts=1e-6) -> None:
 
    

        self.tol= tol
        self.pseudo_counts = pseudo_counts
 
    
    def compute_cn_dist(self,arr):

        dist = np.bincount(arr, minlength=self.max_cn+1)/arr.shape[0]
    
        return dist


    def kl_div(self, p, q, eps=1e-6):
    
        p[(p == 0) | np.isnan(p)] = self.pseudo_counts
        q[(q == 0) | np.isnan(q)] = self.pseudo_counts

        return np.sum(p * np.log2(p / q))


    def scan(self, arr):


            
        start = 0
        segments = []
        for j in range(arr.shape[0]):
            # print(arr[j])
            # if j==90:
            #     print(j)
            if arr[j] > self.tol:

                end = j
                segments.append(([start, end]))
                start = j +1

        segments.append([start, self.nbins-1])
        segments_df = pd.DataFrame(segments, columns=["start", "end"])
        segments_df.index.names = ["segment_id"]
        return segments_df
    
       

    def fit(self, cn_profiles):
        self.cn_profiles = cn_profiles
        self.nbins = self.cn_profiles.shape[1]
        self.max_cn = np.max(self.cn_profiles)

        self.cn_dist_by_bin =np.apply_along_axis(self.compute_cn_dist, axis=0, arr=self.cn_profiles)
        result = []
        #compute pairwise kl divergence
        for i in range(self.nbins - 1):
          
            column_pair_result = self.kl_div(self.cn_dist_by_bin[:, i], self.cn_dist_by_bin[:, i + 1])
            result.append(column_pair_result)

        segs  =self.scan(np.array(result))
        return segs

      


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-i", "--index", type=str, default="genotype",
        help="index of copy number profiles")
    parser.add_argument("-t", "--tolerance",  type=float, default=0.0,
        help="lower bound KL-divergence to establish a new segment")
    parser.add_argument("-p", "--pseudo", required=False, type=float, default=1e-6,
        help="pseudo counts to use for computing KL- divergence")
    parser.add_argument("-o", "--out",  type=str, 
        help="filename of segmentation  to output")
    
    args = parser.parse_args()

    cn_profile_df = pd.read_csv(args.profiles)
# fname = "/scratch/data/leah/phertilizer2.0/sim_study/input//s13_n1000_m15000_c5_p0.25_l0/copy_number_profiles.csv"
# fname ="/scratch/data/leah/phertilizer2.0/DLP/cn_profiles.csv"

    cn_profile_df = cn_profile_df.set_index(args.index)
    cn_prof = cn_profile_df.values

# cn_prof = np.array([[2,2,3,4,2,2],[2,2,3,2,2,2], [2,2,3,2,2,2]])

    segments = Segment().fit(cn_prof)
    segments.to_csv(args.out)




