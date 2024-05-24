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
    

    if args.segfile is not None:
        with open(args.segfile, "r+") as file:
            segments = [int(line.strip()) for line in file]

    else:
        if args.segments is None:
            segments = dat.segments
        
        else:
            segments = args.segments

    print("\t".join(["seg", "snvs", "cn states", "thresh" ] ))
    print("-------------------------------------------")
    for ell in segments:
        print(f"{ell}\t{dat.num_snvs(ell)}\t{dat.num_cn_states(ell)}\t{len(dat.thresholded_cn_prop(ell, thresh=args.thresh_prop))}")
    

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
        k = args.snv_clusters
    else:
         k = len(delta)
    

    
    if args.Tm is not None:
        T_m =nx.DiGraph()
        # edges = [(1,0), (1,2), (0,3)]
        # T_m.add_edges_from(edges)
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
                start_state=(args.root_x, args.root_y), 
                seed = args.seed,
                top_n=  args.top_n,
                collapse= args.collapse,
                order = args.order,
                ninit_segs = args.ninit_segs,
                ninit_Tm = args.ninit_tm,
                cell_threshold= args.cell_threshold,
                max_loops= args.ntree_iter,
                thresh_prop = args.thresh_prop,
                sum_condition = args.sum_condition,
                )

  
    solutions = ph.fit(dat,args.lamb, segments, cores=args.cores, Tm=T_m)



    if args.pickle is not None:
         print("Pickling solutions...")
         pickle_object(solutions, args.pickle)
    
    if args.out is not None:
        print("Drawing clonal trees...")
        for i,sol in enumerate(solutions):
            sol.compute_likelihood(dat, args.lamb)
            sol.png(f"{args.out}/ct{i}.png")
        sol = solutions[0]
        sol.write_flat_files(args.out, data=dat, lamb=args.lamb)


     

    if args.all_sol is not None:
        pickle_object(ph.clonal_trees, args.all_sol)

    if len(solutions) > 0:
        print(f"\nPharming complete with objective value: {solutions[0].cost}!\n")
        print("Let's go sell our trees at the local Pharmer's Market!")
    else:
        print("Error: no trees inferred, check input data and try again!")

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
                        help="number of SNV clusters, if dcfs are also specified, k defaults to the number of specified DCFs")
    parser.add_argument("-D", "--dcfs", type=str, 
                        help="optional filename of dcfs to use")
    parser.add_argument("--delta", type=float, nargs='+',
                        help="list of DCFs to use, ignored if dcf file is provided")
    parser.add_argument("--ninit-segs", type=int,
                        help="number of segments for initialization of mutation cluster tree")
    parser.add_argument("--ninit-tm", type=int,
                        help="number of mutation cluster trees to consider after pruning with initial segs")
    parser.add_argument("--ntree-iter", type=int, default=1,
                        help="number of iterations to check for new mutation cluster tree to use for inference")
    parser.add_argument("--thresh-prop", type=float, default=0.0,
                        help="proportion threshold for determining CN states")
    parser.add_argument("--order", choices=[ 'random','weighted-random', 'nsnvs', 'in-place', 'cost'], default="weighted-random",
                        help="ordering strategy for progressive integration, choose one of 'random', 'weighted-random', 'nsnvs', 'in-place'")
    parser.add_argument("-S", "--cnatrees",type=str,
                        help="optional filename of a pickled dictionary of CNA tree (nx.digraph) for each segment")
    parser.add_argument("--root_x", type=int, default=1,
                        help="starting state for maternal (x) allele")
    parser.add_argument("--root_y", type=int, default=1,
                        help="starting state for paternal (y) allele")
    parser.add_argument("--collapse", action="store_true",
                        help="whether linear chains of copy number events should be collapsed prior to integration")
    parser.add_argument("--sum-condition", action="store_true",
                        help="use the sum condition to filter mutation cluster trees")
    parser.add_argument("--cell-threshold", type=int,
                        help="if collapsing is used, the minimum number of cells a CNA only clone requires to avoid collapsing, NA if not collapsing.")
    parser.add_argument("-L", "--segments", required=False, type=int, nargs='+',
                    help="segment ids of trees to build")
    parser.add_argument("--segfile", required=False, type=str,
                    help="filename with list of segments")
    parser.add_argument("-P" ,"--pickle", required=False, type=str,
                        help="directory where the pickled solution list of top n trees should be saved")
    parser.add_argument("-O" ,"--out", required=False, type=str,
                        help="directory where output files should be written")
    parser.add_argument( "--all-sol", type=str,
        help = "filename of object to pickle all top clonal trees inferred from each mutation cluster tree")
    parser.add_argument("--profile", type=str)


    args = parser.parse_args()
    



    # gtpth = "test"
    # seed = 14
    # cov = 0.25
    # instance = f"s{seed}_m5000_k25_l5_d2"
    # folder = f"n1000_c{cov}_e0" 
    # pth = f"simulation_study/sims"



    # args = parser.parse_args([

    #     "-d", f"{pth}/{instance}/{folder}/data.pkl",
    #     "-j", "10",
    #     # "-D", "test/test_dcfs_s12.txt",
    #     # "-D", f"simulation_study/sims/{instance}/{folder}/dcfs.txt",
    #                   # "-D", f"test/s14_m5000_c0.25_l5/dcfs.txt",
    #      "-D", f"simulation_study/dcf_clustering_gtk/clustsegs10_r100/{instance}/{folder}/dcfs.txt",
    #     # "-T", f"simulation_study/dcf_clustering_gtk/clustsegs10_r100/{instance}/{folder}/Tm.txt",
    #     "-n", "5",
    #     # "-k", "5",
    #     #  "-L", "1",
    #     "-l", "100",
    #     "--ninit-segs", "10",
    #     "--ninit-tm", "10",
    #     "--cell-threshold", "25",
    #     "--order", "cost",
    #     "-s", f"{seed}",
    #     "-P", f"test/s14_m5000_c0.25_l5/solutions.pkl",
    #     "--profile", "test/profile.prof",
    #     "--collapse",
    #     "--ntree-iter", "1",
    #     "--sum-condition",
    #     "-O", "test/s14_m5000_c0.25_l5" #"{gtpth}"

    # ])

    profiler = cProfile.Profile()
    profiler.enable()

    main(args)

    profiler.disable()
    if args.profile is not None:
        profiler.dump_stats(args.profile)

