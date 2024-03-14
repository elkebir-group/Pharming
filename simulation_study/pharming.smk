configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")


rule all:
    # noinspection PyInterpreter
    input:
        expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),


rule pharming:
    input:
        dcfs = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/mut_cluster_tree.txt",
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",

    params:
        lamb = 1e3,
        topn = 3,
        opath = "./pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}"
    threads: 5
    output:
        sol = "pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/solution.pkl",
        profile = "pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/run.prof",
    benchmark:"pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.log",
        err= "pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.err.log"
    shell:
        "python ../src/main.py -d {input.data} -T {input.tm}  -D {input.dcfs} "
        "-s {wildcards.s} "
        "-l {params.lamb} "
        "-n {params.topn} "
        "-j {threads} "
        "-O {params.opath} "
        "--profile {output.profile} "
        "-P {output.sol} > {log.std} 2> {log.err} "

# rule pharming_ilp:
#     input:
#         dcfs = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
#         tm = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/mut_cluster_tree.txt",
#         data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",        
#     params:
#         lamb = 1e3,
#         topn = 3,
#         opath = "./pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}"
#     threads: 7
#     output:
#         sol = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/solution.pkl",
#         profile = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/run.prof",
#     benchmark:"pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
#     log:
#         std= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.log",
#         err= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.err.log"
#     shell:
#         "python ../src/main.py -d {input.data} -T {input.tm}  -D {input.dcfs} "
#         "-s {wildcards.s} "
#         "-l {params.lamb} "
#         "-n {params.topn} "
#         "-j {threads} "
#         "-O {params.opath} "
#         "--profile {output.profile} "
#         "-P {output.sol} > {log.std} 2> {log.err} "
    
rule eval_solutions:
    input:
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
        cellassign = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
        sol = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/solution.pkl",
    params:
        lamb = 1e3,
    output:
        scores = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
    log:
        std= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.log",
        err= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.err.log"
    shell:
        "python ../src/score_tree.py -d {input.data} -t {input.gt} -c {input.cellassign} "
        "-S {input.sol} "
        "-l {params.lamb} "
        "-o {output.scores} > {log.std} 2> {log.err} "




# rule score_tree:
#     input:
#         gt_cell ="input/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/cellclust_gt.csv",
#         gt_mut= "input/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/gt_mut_g{g}.csv",
#         gt_tree = "input/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/true_tree.txt",
#         pred_cell = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/pred_mut_g{g}.csv",
#         pred_mut = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/pred_mut_g{g}.csv",
#         pred_tree = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/tree__g{g}.txt",
    
#     output:"pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/metrics_g{g}.csv"
#     shell:
#         "/scratch/data/leah/phertilizer/simulation_study/scripts/cpp/metrics {input.gt_tree} {input.gt_cell} {input.gt_mut} {input.pred_tree} {input.pred_cell} {input.pred_mut} > {output} "

