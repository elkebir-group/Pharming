import numpy as np
import pandas as pd 
from scipy.spatial import distance
import logging
import seaborn as sns
import matplotlib.pyplot as plt

log = logging.getLogger()
logging.basicConfig(level=logging.WARN)
fname = "dlp/chr1.csv"
dat = pd.read_csv(fname)
dat  = dat.set_index("cn")
# sns.clustermap(dat, cmap="viridis", linewidths=.5)

# Display the plot
# plt.show()
print(dat.head())

matrix= dat.values 

print(matrix.shape)



def prob_distribution(mat, pseudo=0):
    mat_counts = np.sum(mat, axis=1)
    mat_counts =  mat_counts.astype("float")
    mat_counts[mat_counts==0] = pseudo
    probabilities = mat_counts / np.sum(mat_counts)
    return probabilities

def sj_divergence(matrix, start_index, end_index):
    # Extract the segment from the matrix
    segment = matrix[:, start_index:end_index]
    shape1 = segment.shape[1]
    segment = segment.sum(axis=1)
    
    # Calculate the probability distribution of states within the segment
    # segment_probabilities =prob_distribution(segment)
    
    
    # Calculate the probability distribution of states for the rest of the matrix
    rest_matrix = np.hstack((matrix[:, :start_index], matrix[:, end_index:]))
    shape2 = rest_matrix.shape[1]
    rest_matrix = rest_matrix.sum(axis=1)

    assert matrix.shape[1] == (shape1 + shape2)
    # rest_probabilities = prob_distribution(rest_matrix)

    jsd = distance.jensenshannon(segment, rest_matrix, base=2)

    
    return jsd

def cbs(x, ent_thresh, bin_thresh=0):
    max_entropy = np.NINF
    for i in range(x.shape[1]):
        seg_entropy, j = find_max_entropy_segment(i, x, bin_thresh)
        if seg_entropy >= max_entropy:
            s = i 
            e = j
            max_entropy = seg_entropy
    # print(f"{s} to {e}: {max_entropy}")
    thresh = ent_thresh  <= max_entropy

    return thresh, max_entropy, s, e,


def find_max_entropy_segment(start, x, bin_thresh):


    max_entropy = np.NINF
    best = x.shape[1]
    for j in range(start + bin_thresh, x.shape[1]):
        div = sj_divergence(x, start, j)
        if div > max_entropy:
            best = j 
            max_entropy = div 
    return max_entropy, best 
    


   


def rsegment(x, start, end, L, ent_thresh, bin_thresh):
    '''Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''

    threshold, ent, s, e = cbs(x[:,start:end], ent_thresh=ent_thresh, bin_thresh=bin_thresh)
    #candidate
    #[s,e)
    log.info('Proposed partition of {} to {} from {} to {} with entropy {} is {}'.format(start, end, start+s, start+e, ent, threshold))
    if (not threshold)  | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            #recurse on [0,i)
            rsegment(x, start, start+s, L,ent_thresh, bin_thresh)
        if e-s > 0:
            #recurse [i,j)
            rsegment(x, start+s, start+e, L,ent_thresh, bin_thresh)
        if start+e < end:
            #recurse [j,nbins)
            rsegment(x, start+e, end, L,ent_thresh, bin_thresh)
    return L


def scan(arr):
    entropies = []            
    for i in range(arr.shape[1]-1):
        p = arr[:,i]
        q = arr[:,i+1]
        entropies.append(distance.jensenshannon(p, q, base=2))
    return np.array(entropies)
         
      

def segment(x, ent_thresh=0.05, bin_thresh=1):
    '''Segment the matrix by looking for large dviergence in prob distributions between segments
    '''
    entropies  =scan(x)
    ent_thresh =np.quantile(entropies, q=0.75)
    start = 0
    end = x.shape[1]
    L = []
    rsegment(x, start, end, L, ent_thresh=ent_thresh, bin_thresh=1)
    return L


def segment_scan(matrix, start_bin, quant=0.75, bin_thresh=1):
        x =  scan(matrix)
        ent_thresh =np.quantile(x, q=quant)
        start = 0
        segments = []
        for j in range(x.shape[0]):

            if x[j] > ent_thresh:

                end = j
                segments.append([start+start_bin, start_bin+ end+1])
                start = end + 1

        segments.append([start +start_bin, start_bin+ len(x)+1])
        segments_df = pd.DataFrame(segments, columns=["start", "end"])
    
        return segments_df


fname = "dlp/counts.csv"
df = pd.read_csv(fname)
print(df.head())
chrom = df["chr"].unique()
chrom_segs = []
for c in chrom:
    dat = df[df["chr"]==c] 
    dat =dat.drop(["chr"], axis=1)
    dat = dat.set_index("cn") 
    start_bin = dat["bin"].min()
    dat = dat.pivot(columns="bin", values="counts")
    mat = dat.values 
    seg_df = segment_scan(mat, start_bin, quant=0.90)

    seg_df["chr"] = c 
    # print(seg_df)
    chrom_segs.append(seg_df)

segmentation = pd.concat(chrom_segs)
segmentation["end"] = segmentation["end"] -1 

segmentation.to_csv("dlp/segmentation.csv", index=False)








# log.setLevel(logging.INFO)    

# L_prime = segment_scan(entropies, ent_thresh=ent_thresh)
# L =segment(matrix, ent_thresh=0.05)

# print(L_prime)
# print(L)

