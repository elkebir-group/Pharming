from dataclasses import dataclass
import numpy as np
import pandas as pd



'''
N = number of cells
M = number of SNVs
G = number of segments
'''

@dataclass
class Data:
    var : np.array  # N x M matrix with the number of variant read counts
    total : np.array #  N x M matrix with the number of total read counts
    copy_numbers: np.array # N x G matrix with the copy number profile of each cell in each segment
    snv_to_seg: np.array  # an M length array which holds the segment of an SNV 
    cell_lookup : pd.Series # an N length series mapping internal cell index to input cell label
    mut_lookup : pd.Series #an M length series mapping interal SNV index to input SNV label