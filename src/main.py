
# Created by: L.L. Weber
# Created on: 2024-02-29 17:30:46

import argparse
from data import Data
from pharming import Pharming
from utils import load_pickled_object, pickle_object
import networkx as nx

def main(args):
    print("\nWelcome to the Pharm! Let's start pharming.....\n")
    
    if args.data is not None:
        dat = load_pickled_object(args.data)
    elif args.file and args.copy_numbers is not None:
        dat = load_pickled_object(args.file, args.copy_numbers )
    else:
        IOError("Either both read counts and copy number files \
                must be specified or alternatively, the path to the \
                preprocessed pharming data object!")
    


    if args.segments is None:
        segments = dat.segments
    
    else:
        segments = args.segments

    print("\t".join(["seg", "snvs", "cn states"] ))
    print("-------------------------------------------")
    for ell in segments:
        print(f"{ell}\t{dat.num_snvs(ell)}\t{dat.num_cn_states(ell)}")
    
    # segments = [ell for ell in segments if dat.num_snvs(ell) > 75]


    if args.delta is not None:
        dlist = {}
        delta = {i:  dlist[i] for i in range(len(dlist))}
    else:
        delta =None 
    if args.dcfs is not None:
        delta = {}
        with open(args.dcfs, "r+") as file:
            for i, line in enumerate(file):
                try:
                    delta[i] = float(line.strip())
                except:
     
                     raise ValueError("DCF file is not properly formatted.")

    
    if args.cnatrees is not None:
        cnatrees = load_pickled_object(args.cnatrees)
    else:
         cnatrees = None 
    
    if delta is None:
        k = args.k 
    else:
         k = len(delta)
    

    
    if args.Tm is not None:
        T_m =nx.DiGraph()
        with open(args.Tm, "r+") as file:
            for line in file:
                edges = line.strip().split("\t")
                edges = [int(item) for item in edges]
                T_m.add_edge(*edges)
    else:
        T_m = None 

    ph = Pharming(dcfs = delta,
                cnatrees = cnatrees,
                k= k, 
                T_m = T_m,
                start_state=(args.root_x, args.root_y), 
                seed = args.seed,
                top_n=  args.top_n,
                collapse= args.collapse
                )

  
    solutions = ph.fit(dat,args.lamb, segments, cores=args.cores, ninit_segs=args.ninit_segs)


    if args.scores is not None:
            print("Saving cost values...")
            with open(args.scores, "w+") as file:
                for sol in solutions:
                    file.write(f"{sol.cost}\n")
    
    if args.pickle is not None:
         print("Pickling solutions...")
         pickle_object(solutions, args.pickle)
    
    if args.out is not None:
        print("Drawing clonal trees...")
        for i,sol in enumerate(solutions):
            sol.png(f"{args.out}/ct{i}.png")


     

    if args.all_sol is not None:
        pickle_object(ph.clonal_trees, args.all_sol)


    print("\nPharming complete...let's go sell our trees at the local Pharmer's Market!")

import cProfile

if __name__ == "__main__":

  
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--data", required=False,
                        help="input file of preprocessed data pickle")
    parser.add_argument("-f", "--file", required=False,
                        help="input file for variant and total read counts with unlabled columns: [chr segment snv cell var total]")
    parser.add_argument("-c" ,"--copy-numbers", required=False,
                        help="input files of copy numbers by segment with unlabeled columns [segment cell totalCN]")
    parser.add_argument("-s" ,"--seed", required=False, type=int, default=1026,
                        help="random number seed (default: 1026)")
    parser.add_argument("-j" ,"--cores", required=False, type=int,default=1,
                        help="Max number of cores to use for inferring segment trees")
    parser.add_argument("-l" ,"--lamb", required=False, type=float, default=1e3,
                        help="lambda value, default=1e3")
    parser.add_argument("-n" ,"--top_n", required=False, type=int,default=3,
                        help="number of trees to retain in each step")
    parser.add_argument("-T", "--Tm", type=str,  
                        help="optional filename of  mutation cluster tree")
    parser.add_argument("-k", "--snv-clusters", type=int, default=5,
                        help="number of SNV clusters, if dcfs are also specified, k defaults to the \
                        of specified DCFs")
    parser.add_argument("-D", "--dcfs", type=str, 
                        help="optional filename of dcfs to use")
    parser.add_argument("--delta", type=float, nargs='+',
                        help="list of DCFs to use, ignored if dcf file is provided")
    parser.add_argument("--ninit-segs", type=int, default=3,
                        help="default number of segments for initialization of mutation cluster tree")
    parser.add_argument("-S", "--cnatrees",type=str,
                        help="optional filename of a pickled dictionary of CNA tree (nx.digraph) for each segment")
    parser.add_argument("--root_x", type=int, default=1,
                        help="starting state for maternal (x) allele")
    parser.add_argument("--root_y", type=int, default=1,
                        help="starting state for paternal (y) allele")
    parser.add_argument("--collapse", action="store_true",
                        help="whether linear chains of copy number events should be collapsed prior to integration")
    parser.add_argument("-L", "--segments", required=False, type=int, nargs='+',
                    help="segment ids of trees to build")
    parser.add_argument("-P" ,"--pickle", required=False, type=str,
                        help="directory where the pickled solution list of top n trees should be saved")
    parser.add_argument("-O" ,"--out", required=False, type=str,
                        help="directory where output files should be written")
    parser.add_argument("-J", "--scores", type=str,
        help = "filename of tree scores")
    parser.add_argument( "--all-sol", type=str,
        help = "filename of object to pickle all top clonal trees inferred from each mutation cluster tree")
    parser.add_argument("--profile", type=str)


    args = parser.parse_args()
    


    instance = "s11_m5000_k25_l5"
    # instance = "s12_m5000_k25_l7"
    folder = "n1000_c0.05_e0" 
    pth = f"simulation_study/input"

    gtpth = "test"



    args = parser.parse_args([

        "-d", f"{pth}/{instance}/{folder}/data.pkl",
        "-j", "4",
        "-D", f"{pth}/{instance}/{folder}/dcfs.txt",
        "-T", f"{pth}/{instance}/{folder}/Tm.txt",
        "-n", "3",
        # "-L",  "0", "14", "22",
        "--ninit-segs", "2",
        "-s", "11",
        # "--segment", "0",
        # "--out", f"/Users/leah/Documents/Research/projects/Pharming/test",
        "-J", f"{gtpth}/scores_test.csv",
        "-P", f"{gtpth}/solution_test.pkl",
        "--all-sol", f"{gtpth}/clonal_trees.pkl",
        "--profile", "test/profile.prof",
        "--collapse",
        "-O", f"{gtpth}"

    ])

    profiler = cProfile.Profile()
    profiler.enable()

    # Call your main function here with the parsed arguments
    main(args)

    profiler.disable()
    if args.profile is not None:
        profiler.dump_stats(args.profile)

