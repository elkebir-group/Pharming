configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")

# #TODO: add rule generate cna trees 
ruleorder:  simulate > generatesinglecells

rule all:
    input:
        expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
        expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),



rule simulate:
    input: "cnatrees_nocompleteloss.txt"
    output:
        tree =  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        prop=  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/proportions.tsv",
        genotypes ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
    params:
        pth = config["simpath"],
        cp_thresh = 0.05,
        purity = 0.99,
        sample = 1,
        alpha= 0.001,
        truncalSegs = 2,
        simout_dir = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}",
    log:
        std ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/run.log", 
        err ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/err.log" 
    shell:
        "{params.pth} -r -S {input} "
        " -purity {params.purity} -minProp {params.cp_thresh} "
        " -kk {params.truncalSegs} -f "
        "-s {wildcards.s}  -l {wildcards.mclust} "
        "-k {wildcards.nsegs} -n {wildcards.snvs} -m {params.sample} "
        "-output_file_dir {params.simout_dir}  > {log.std} 2> {log.err}  "


rule generatesinglecells:
    input:
        tree =  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        prop=  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/proportions.tsv",
        genotypes ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
    output:
        phi ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        sparse ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cells.p0",
    params:
        pth = config["genpath"],
        purity = 0.99,
        sample = 1,
        alpha= 0.001,
        simout_dir = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}",
        scout_dir = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}",
    log:
        std ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/run.log", 
        err ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/err.log" 
    shell:
     "{params.pth} -num_cells {wildcards.cells} -read_depth {wildcards.cov} -k {wildcards.nsegs} "
     "-alpha_fp {params.alpha} -out_dir {params.scout_dir} -in_dir {params.simout_dir} -e {wildcards.err} "
     " -m {params.sample} > {log.std} 2> {log.err} "

rule make_data:
    input: 
        sparse ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cells.p0",
    output: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    shell:
        "python ../src/make_data.py -f {input.sparse} -c {input.copy_profiles} -D {output} " 

rule make_gt:
    input: 
        phi ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        tree =  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        genotypes ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
        data = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    params: 
        lamb = 1e3
    output:
        png = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.png",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
        phi = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/phi.pkl",
        dcfs ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/Tm.txt",
    log:
        std ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.log",
        err ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.err.log",
    shell:
        "python ../src/ground_truth_clonal_tree.py -d {input.data}  -t {input.tree} -p {input.phi} "
        "-g {input.genotypes} -T {output.gt} --draw {output.png}  -D {output.dcfs} -P {output.phi} "
        "-l {params.lamb} --mut-cluster-tree {output.tm} > {log.std} 2> {log.err} "



    










        






