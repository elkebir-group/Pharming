configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")


rule all:
    input:
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
            inpath = config["inpath"],
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
        dcfs = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/mut_cluster_tree.txt",
        data= "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",

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


    
rule eval_solutions:
    input:
        data= "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        gt = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
        cellassign = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
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
#         gt_cell ="{inpath}/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/cellclust_gt.csv",
#         gt_mut= "{inpath}/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/gt_mut_g{g}.csv",
#         gt_tree = "{inpath}/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/true_tree.txt",
#         pred_cell = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/pred_mut_g{g}.csv",
#         pred_mut = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/pred_mut_g{g}.csv",
#         pred_tree = "pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/tree__g{g}.txt",
    
#     output:"pharming/s{seed}_n{cells}_m{mutations}_c{clones}_p{prob}_l{loss}/SegTrees/metrics_g{g}.csv"
#     shell:
#         "/scratch/data/leah/phertilizer/simulation_study/scripts/cpp/metrics {input.gt_tree} {input.gt_cell} {input.gt_mut} {input.pred_tree} {input.pred_cell} {input.pred_mut} > {output} "

