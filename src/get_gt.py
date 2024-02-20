from data import Data, load_from_pickle
from clonal_tree import ClonalTree, load 
from cell_mapping import CellAssign
import argparse
import pandas as pd 

def dict_to_df(mydict, colnames):
    df = pd.DataFrame.from_dict(mydict, orient='index', columns=[colnames[1]])
    df.reset_index(inplace=True)
    df.columns = colnames
    return df 






if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-t", "--tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-p" ,"--psi", required=True, type=str,
                        help="where the output file should be written")
    parser.add_argument("-f" ,"--delta", required=True, type=str,
                        help="where the output file should be written")
    
    args = parser.parse_args()
    dat = load_from_pickle(args.data)
    phi = load(args.cell_assign)
    gt = load_from_pickle(args.tree)
    gt.reindex_snvs(dat.mut_lookup)
    gt_psi = gt.get_psi()
    psi_df  = dict_to_df(gt_psi, ["snv", "node"])
    psi_df.to_csv(args.psi, index=False)
    phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    df = dict_to_df(gt_dcfs, ["node", "dcf"])
    df.to_csv(args.delta, index=False)