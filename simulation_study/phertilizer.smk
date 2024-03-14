
configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")

# #TODO: add rule generate cna trees 
# ruleorder:  simulate > generatesinglecells

rule all:
    # noinspection PyInterpreter
    input:
        expand("phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),


    
rule make_data:
    input: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    output: 
        readcounts = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    shell:
        "python ../src/prep_phertilizer.py -d {input} -f {output.readcounts} -u {output.umap}




rule phertilizer:
    input:
        readcounts = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    output:
        pickle = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.pkl",
        tree = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.txt",
        predcell = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        predmut = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        png = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.png",
        like = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/likelihood.txt",
    benchmark:"phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.log",
        err= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.err.log"
    shell:
        "phertilizer -f {input.readcounts} --bin_count_data {input.umap} --no-umap "
        "--tree_pickle {output.pickle} --tree {output.png} --tree_text --likelihood {output.like} "
        " -n {output.predcell} - m {output.predmut}  --post_process -s {wildcards.seed} "
        " > {log.std} 2> {log.err} "


    
# rule eval_solutions:
#     input:
#         data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
#         gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
#         cellassign = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
#         sol = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/solution.pkl",
#     params:
#         lamb = 1e3,
#     output:
#         scores = "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
#     log:
#         std= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.log",
#         err= "pharming_ilp/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.err.log"
#     shell:
#         "python ../src/score_tree.py -d {input.data} -t {input.gt} -c {input.cellassign} "
#         "-S {input.sol} "
#         "-l {params.lamb} "
#         "-o {output.scores} > {log.std} 2> {log.err} "







        






