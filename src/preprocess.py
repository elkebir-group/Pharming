from clonal_tree_new import ClonalTreeNew
from data import Data
import argparse 
import pandas as pd 
import numpy as np 
import networkx as nx
from genotype import genotype

def data_prep(var_fname, copy_fname):
    col_names = ['segment', 'mutation_label', 'cell_label','var', 'total']



    # read_counts = pd.read_table(
    #     var_fname, sep="\t")
    with open(var_fname, "r+") as file:
        firstline =True
        sparse = []
   
        for line in file:
            if firstline:
                firstline = False
                continue   
            else:
               vals = line.strip().split("\t")
               vals = [int(v) for v in vals]
               sparse.append(vals)
        
    read_counts  = pd.DataFrame(sparse, columns=col_names)
   
    # read_counts = pd.read_table(
    #     var_fname, header=None, names=col_names, skiprows=[0])
    
    # print(read_counts.head())

    copy_numbers = pd.read_csv(copy_fname, header=None,names=["segment", "cell_label", "x", "y"], skiprows=[0])

 
    
    # read_counts['chr_mutation'] = read_counts['mutation_label']

    cell_labels = np.sort(read_counts['cell_label'].unique())
    mut_labels = np.sort(read_counts['mutation_label'].unique())



    #create indexed series of mapping of cell index to label
    cell_lookup = pd.Series(data=cell_labels, name="cell_label").rename_axis("cell")     
    mut_lookup = pd.Series(data=mut_labels, name="mutation_label").rename_axis("mut")

    read_counts = pd.merge(read_counts, cell_lookup.reset_index(), on='cell_label', how='left')
    read_counts = pd.merge(read_counts, mut_lookup.reset_index(), on='mutation_label', how='left')
    print(read_counts.head())
    #in long format
    segs = copy_numbers.loc[:, ["segment"]].drop_duplicates()
    seg_labels = np.sort(segs['segment'].unique())
    seg_lookup = pd.Series(data=seg_labels, name="segment").rename_axis("seg_id")

    copy_numbers = pd.merge(copy_numbers, cell_lookup.reset_index(), on='cell_label', how='left').drop("cell_label", axis=1)
    copy_numbers= pd.merge(copy_numbers, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)


    read_counts = pd.merge(read_counts, seg_lookup.reset_index(), on='segment', how='left').drop("segment", axis=1)
    seg_to_mut_mapping = read_counts.loc[:, ["seg_id", "mut"]].drop_duplicates()
    snv_to_seg = seg_to_mut_mapping.set_index("mut")["seg_id"].to_dict()
    seg_to_snvs =  {value: [k for k, v in snv_to_seg.items() if v == value] for value in set(snv_to_seg.values())}


    read_counts= read_counts.set_index(["cell", "mut"])
    var = read_counts["var"].unstack(level="mut", fill_value=0).to_numpy()
    total = read_counts["total"].unstack(level="mut", fill_value=0).to_numpy()
    copy_numbers = copy_numbers.set_index(["seg_id", "cell"])
    copy_numbers= copy_numbers.unstack(level="seg_id", fill_value=0).to_numpy()

    return Data(var, total, copy_numbers, snv_to_seg, seg_to_snvs, cell_lookup, mut_lookup)


def genotypes_prep(genotypes_fname, genotypes):
    firstline = True

    with open(genotypes_fname, "r+") as file:
        for line in file:
            if firstline:
                firstline= False
                continue
            else:
                row = line.strip().split("\t")
                row = [int(r) for r in row]
                m= row[4]
                g= genotype(row[2], row[3], row[5], row[6])
                genotypes[row[0]][row[1]][m] = g
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c", "--profiles", type=str,
        help="filename of input copy number profiles")
    parser.add_argument("-t", "--tree", type=str,
        help="filename of input tree")
    parser.add_argument("--phi", type=str, help="ground truth cell assignments")
    parser.add_argument("-g" ,"--genotypes", type=str, help="ground truth genotypes")

    parser.add_argument("-D", "--data",  type=str, 
        help="output filename of pickled data object")
    parser.add_argument("-T", "--clonal-tree",  type=str, 
        help="output filename of pickled clonal tree object")
    parser.add_argument( "--draw",  type=str, 
        help="png or pdf of output tree")


    args = parser.parse_args()
    # pth = "/scratch/data/leah/pharming/simulation_study/input/s10_n500_m1000_k50_c1_l5"
    # args = parser.parse_args(
    #     ["-f", f"{pth}/sparse.p0",
    #      "-c" ,f"{pth}/cells.p0",
    #      "-g", f"{pth}/gentoypes.txt",
    #      "--phi", f"{pth}/cellAssignments.p0",
    #      "-t", f"{pth}/tree.txt",
    #      "-T", f"{pth}/ground_truth.pickle",
    #      "-D", f"{pth}/pharming_data.pickle",
    #      "--draw", f"{pth}/tree.png"
    #     ]
    # )

    data = data_prep(args.file, args.profiles)
    print(data)
    if args.data is not None:
        data.save(args.data)

    tree = nx.DiGraph()
    with open(args.tree, "r+") as file:
        firstline =True
   
        for line in file:
            if firstline:
                firstline = False
                continue
                
            else:
               vals = line.strip().split("\t")
               vals = [int(v) for v in vals]
               tree.add_edge(vals[0], vals[1])
    genotypes = {v:  {} for v in tree}
    for v in tree:
        for s in data.segments:
            genotypes[v][s] = {}
    genotypes_prep(args.genotypes, genotypes)
    cell_mapping = {v: [] for v in tree}
    if args.phi is not None:
        cell_assignment = pd.read_csv(args.phi)
        # print(cell_assignment.head())
        for index, row in cell_assignment.iterrows():
            i = row['Cell']
            v = row['Cluster']
            cell_mapping[v].append(i)


    ct = ClonalTreeNew(tree, genotypes, cell_mapping=cell_mapping )
    missing_muts = set(data.muts) - set(ct.get_all_muts())

    if args.phi is not None:
        cost2 = ct.compute_pooled_costs(data)
        costs = ct.compute_costs(data)
        print(f"cost 1: {costs} cost 2: {cost2}")
        



    if args.draw is not None:
        ct.draw(args.draw)
    if args.clonal_tree is not None:
        ct.save(args.clonal_tree)








