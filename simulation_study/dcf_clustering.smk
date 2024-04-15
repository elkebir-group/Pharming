configfile: "simulate.yml"

seeds = [i+10 for i in range(config["nseeds"])]
import sys 
sys.path.append("../src")

rule all:
    input:
        expand("dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
               s=seeds,
               cells=config["cells"],
               snvs=config["snvs"],
               nsegs=config["nsegs"],
               cov=config["cov"],
               mclust=config['mclust'],
               err=config["cerror"],
               dirch = config["dirch"]
        ),
        # expand("dcf_clustering/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/output_unconstrained.pkl",
        #        s=seeds,
        #        cells=config["cells"],
        #        snvs=config["snvs"],
        #        nsegs=config["nsegs"],
        #        cov=config["cov"],
        #        mclust=config['mclust'],
        #        err=config["cerror"]
        # )

rule run_dcf_clustering_constrained:
    input: 
        data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
        dcfs= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
    output: 
        output_data= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/output.pkl",
        dcfs= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt"
    params: 
        restarts=50
    threads: 3
    log: 
        std= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.log",
        err= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/err.log"
    benchmark: "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    shell: 
        "python ../src/dcf_clustering_v2.py -d {input.data} -g {input.dcfs} -o {output.output_data} -D {output.dcfs} "
         "-r {params.restarts} -c -j {threads} > {log.std} 2> {log.err}"    


rule clustering_eval:
    input:         
        gt_dcfs = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        inf_dcfs="dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        result = "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/output.pkl",
        gt_tree= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl"  
    output: 
        scores = "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
        dcfs = "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf_dcfs.txt"
    run:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        import pandas as pd 

        def read_dcfs(fname):
            dcfs = []
            with open(fname, "r+") as file:
                for line in file:
                    line = line.strip()
                    if line[0]=="#":
                        continue
                    else:
                        dcfs.append(float(line))
            return np.array(dcfs)
        
        def get_likelihood(fname):
            #likelihood: -19364.18433407771
            #gt likelihood: -19296.746673797636
            linecount = 0
            with open(fname, "r+") as file:
                for line in file:
                    if linecount > 1:
                        break
                    line = line.strip().split(" ")
                    if linecount ==0:
                        likelihood = float(line[1])
                    else:
                        gt_like = float(line[2])
                    linecount += 1
            return likelihood, gt_like
          
        def compute_mean_difference(ground_truth, dcfs):

            differences = np.abs(np.subtract.outer(ground_truth, dcfs))
            row_ind, col_ind = linear_sum_assignment(differences)
            selected_differences = differences[row_ind, col_ind]
        
         
            return selected_differences
        
        gt = read_dcfs(input['gt_dcfs'])
        inf = read_dcfs(input['inf_dcfs'])
        selected_differences = compute_mean_difference(gt, inf)
        max_diff = np.max(selected_differences)

        mean_diff = np.mean(selected_differences)
        min_diff = np.min(selected_differences)


        inf_like, gt_like = get_likelihood(input["inf_dcfs"])
 
        result = pd.read_pickle(input["result"])
        gt_tree = pd.read_pickle(input["gt_tree"])

    
        def compare_CNA_trees(gt, S, segment):
    

            S_gt = gt.get_cna_tree(segment)
            if len(S_gt) == 1:
            
                return None
            return set(S_gt.edges) == set(S)
        cna_trees = result[2]
        cna_trees_correct = [compare_CNA_trees(gt_tree,cna_trees[ell], ell ) for ell in cna_trees]
        cna_trees_correct= [x for x in cna_trees_correct if x is not None]
        perc_correct = sum(cna_trees_correct)/len(cna_trees_correct)

        res = [[mean_diff,max_diff, min_diff, inf_like, gt_like, perc_correct]]
        print(res)

        df = pd.DataFrame(res, columns = ["mad","max_diff", "min_diff", "inf_likelihood", "gt_likelihood", "perc_cna_tree_correct" ])
        df.to_csv(output["scores"], index=False)
        
        with open(output["dcfs"], "w+") as file:
            for d in inf:
                file.write(f"{d}\n")





# rule run_dcf_clustering_unconstrained:
#     input: 
#         data= "input/input/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#         dcfs= "input/input/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
#         accuracy= "dcf_clustering_v2/accuracy.csv"

#     output: 
#         output_data= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/output_unconstrained.pkl"

#     params: 
#         restarts=1

#     log: 
#         std= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run_unconstrained.log",
#         err= "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/err_unconstrained.log"

#     shell: 
#         "python ../src/dcf_clustering.py {input.data} {input.dcfs} {output.output_data} {input.accuracy} {params.restarts} 0 > {log.std} > {log.err}" 