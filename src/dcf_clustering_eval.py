import numpy as np
from scipy.optimize import linear_sum_assignment


def compute_mean_difference(ground_truth, dcfs):

    differences = np.abs(np.subtract.outer(ground_truth, dcfs))
    row_ind, col_ind = linear_sum_assignment(differences)
    selected_differences = differences[row_ind, col_ind]
    mean_difference = np.mean(selected_differences)

    return mean_difference
