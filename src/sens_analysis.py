import clonal_tree_new
import data 
import numpy as np 
import argparse 
import pandas as pd 
from copy import deepcopy




def introduce_cell_errors(ct,rate, rng):
    phi = ct.phi 
    clones = set(ct.clones())

    for i, v in phi.items():
        if rng.random() < rate:
            candidates = clones - {v}
            clone_choices = np.array(list(candidates))
            phi[i] = rng.choice(clone_choices)


    ct.phi_to_cell_mapping(phi)

def introduce_geno_errors(ct, rate, rng):
    psi = ct.psi
    clones = set(ct.clones())
    for j, v in psi.items():
        j_seg = ct.mut2seg[j]
        # geno = ct.genotypes[v][ct.mut2seg[j][j]]

        if rng.random() < rate:
            candidates = clones - {v}
            clone_choices = np.array(list(candidates))
            clone_choices = [c for c in clone_choices  if len(ct.mut_mapping[c]) > 0]
            psi[j] = rng.choice(clone_choices)
            cand_muts = ct.mut_mapping[psi[j]]

            j_prime = rng.choice(cand_muts) 
            j_prime_seg = ct.mut2seg[j_prime]
            for u in ct.tree:
                new_geno = ct.genotypes[u][j_prime_seg][j_prime]
                ct.genotypes[u][j_seg][j] = new_geno

            



def prep_error_rates(inputs):
    error_rates  = inputs.split(":")
    vals =  tuple([int(e) for e in error_rates])

    error_rates = np.arange(*vals)/100
    if 0 not in error_rates:
        error_rates = np.insert(error_rates, 0, 0)
    return error_rates




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--data",  type=str, required=True,
        help="output filename of pickled data object")
    parser.add_argument("-T", "--clonal-tree",  type=str,  required=True,
        help="output filename of pickled clonal tree object")
    parser.add_argument("-c", "--cell-errors",  type=str, default='5:10',
        help="x:y:z:... where each value specifies a cells error rate to apply")
    parser.add_argument("-m", "--mut-errors",  type=str, default='5:10',
        help="x:y:z where x is start, y is end (not inclusive) amd z is step, works like np.arange")
    parser.add_argument("-s", "--seed",  type=int, default=1026,
        help="random number seed")
    parser.add_argument("-r", "--reps",  type=int, default=10,
        help="number of repetitions of each error rate combo")
    parser.add_argument("-o", "--out-file",  type=str,
        help="output filename of daframe results")
    
    args = parser.parse_args()
    # pth = "/scratch/data/leah/pharming/simulation_study/input/s10_n500_m1000_k50_c1_l5"
    # args = parser.parse_args([
    #      "-T", f"{pth}/clonal_tree.pickle",
    #      "-D", f"{pth}/data.pickle",
    #      "-c", "5:105:5",
    #      "-m", "5:105:5",
    #      "-s",  "10",
    #     "-r", "5",
    #     "-o", "test/res.csv"
    #     ]
    # )

    

    dat = data.load_from_pickle(args.data)
    ct = clonal_tree_new.load(args.clonal_tree)
    gt_cell_mapping = ct.cell_mapping.copy()

    rng = np.random.default_rng(args.seed)

    c_error_rates = prep_error_rates(args.cell_errors)
    m_error_rates = prep_error_rates(args.mut_errors)

    data_rows = []
    base1 = ct.compute_costs(dat)
    base2 = ct.compute_pooled_costs(dat)
    # data_rows.append([-1, 0, base1, base2])
    ct_orig = deepcopy(ct)
    for c in c_error_rates:
        for m in m_error_rates:
        
            for n in range(args.reps):

                ct = ct_orig
                introduce_cell_errors(ct, c, rng)
                introduce_geno_errors(ct, m, rng)
                c1 = ct.compute_costs(dat)
                c2 = ct.compute_pooled_costs(dat)
                row = [n, c, m, c1, c2]
                data_rows.append(row)

                if c==0 and m==0:
                    break 
        # print(f"Cell error rate: {c} cost1: {c1} cost2: {c2}")


    res = pd.DataFrame(data_rows, columns=['trial', 'cell_error',  'mut_error', 'cost_indiv', 'cost_pooled'])

    res['ncells'] = len(ct.get_all_cells())
    res['nmuts']  = len(ct.get_all_muts())
    res['clones'] = len(ct.clones())
    res['base_indiv'] = base1
    res['base_pooled'] = base2
    # print(res.head())


    if args.out_file is not None:
        res.to_csv(args.out_file, index=False)






