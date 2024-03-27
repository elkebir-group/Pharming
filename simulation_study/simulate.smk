configfile: "test.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")

# #TODO: add rule generate cna trees 
ruleorder:  simulate > generatesinglecells

rule all:
    input:
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
            inpath = config["inpath"],
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
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        prop=  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/proportions.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
    params:
        pth = config["simpath"],
        cp_thresh = 0.05,
        purity = 0.99,
        sample = 1,
        dirch = config["dirch"],
        alpha= config["alpha"],
        truncalSegs = 2,
        simout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/run.log", 
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/err.log" 
    shell:
        "{params.pth} -r -S {input} "
        " -purity {params.purity} -minProp {params.cp_thresh} "
        " -kk {params.truncalSegs} -f "
        "-dirich_param {params.dirch} "
        "-s {wildcards.s}  -l {wildcards.mclust} "
        "-k {wildcards.nsegs} -n {wildcards.snvs} -m {params.sample} "
        "-output_file_dir {params.simout_dir}  > {log.std} 2> {log.err}  "


rule generatesinglecells:
    input:
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        prop=  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/proportions.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
    output:
        phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        sparse ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cells.p0",
    params:
        pth = config["genpath"],
        purity = 0.99,
        sample = 1,
        alpha= config["alpha"],
        simout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}",
        scout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/run.log", 
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/err.log" 
    shell:
     "{params.pth} -num_cells {wildcards.cells} -read_depth {wildcards.cov} -k {wildcards.nsegs}  "
     "-alpha_fp {params.alpha} -out_dir {params.scout_dir} -in_dir {params.simout_dir} -e {wildcards.err} "
     " -m {params.sample} > {log.std} 2> {log.err} "

rule make_data:
    input: 
        sparse ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cells.p0",
    output: "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    shell:
        "python ../src/make_data.py -f {input.sparse} -c {input.copy_profiles} -D {output} " 

rule make_gt:
    input: 
        phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
        data = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    params: 
        lamb = 1e3
    output:
        png = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.png",
        gt = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
        phi = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/phi.pkl",
        dcfs ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/Tm.txt",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.log",
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.err.log",
    shell:
        "python ../src/ground_truth_clonal_tree.py -d {input.data}  -t {input.tree} -p {input.phi} "
        "-g {input.genotypes} -T {output.gt} --draw {output.png}  -D {output.dcfs} -P {output.phi} "
        "-l {params.lamb} --mut-cluster-tree {output.tm} > {log.std} 2> {log.err} "



    










        






