import argparse
from data import load_from_pickle
from clonal_tree import ClonalTree, load
from cell_mapping import CellAssign
import networkx as nx 




    
if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--tree", required=False,
                        help="input file of ground truth tree")
    parser.add_argument("-c", "--cell-assign", required=False,
                        help="input file of ground truth  cell assignment")
    parser.add_argument("-d" ,"--dcfs", required=False, type=str)
    parser.add_argument("-T" ,"--mut_cluster_tree", required=False, type=str)

    args = parser.parse_args()


    gt = load_from_pickle(args.tree)

 
    phi = load(args.cell_assign)

    
   
    # phi.relabel(dat.cell_lookup)
    gt_dcfs = gt.compute_dcfs(phi)
    root = gt.root


    gt_T_m = gt.mutation_cluster_tree()
    gt_T_m.remove_node(root)




    mapping = {}
    nodes = list(gt_T_m.nodes)
    nodes.sort()
    for i,n in enumerate(nodes):
          mapping[n] =i 
    
    T_m = nx.relabel_nodes(gt_T_m, mapping)
    if args.mut_cluster_tree is not None:
        with open(args.mut_cluster_tree, "w+") as file:
            for u,v in T_m.edges:
                file.write(f"{u}\t{v}\n")

    # utils.pickle_object(T_m, "test/T_m.pkl")
    # utils.draw(T_m, "test/input_0.05/T_m.png")

    gt_delta = {mapping[n]: gt_dcfs[n] for n in mapping if n != root}

    if args.dcfs is not None:
        with open(args.dcfs, "w+") as file:
            for key,val in gt_delta.items():

                file.write(f"{val}\n")

   