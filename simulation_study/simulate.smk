
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
        # expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     err = config["cerror"]
        # ),
        # expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/snv_assignment.csv",
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     err = config["cerror"]
        # ),
        expand("pharming/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/solution.pkl",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
        # expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/lamb{lamb}/obj1_cell_accuracy.csv",
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     lamb = config["lamb"],
        #     err = config["cerror"]
        # ),
        # expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/lamb0/cna_only_cell_accuracy.csv",
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     err = config["cerror"]
        # ),

        

        

        # expand("input/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        #     s = seeds,
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     mclust = config['mclust']
        # ),
        


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

rule make_gt:
    input: 
        tree =  "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/tree.tsv",
        genotypes ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/node.tsv",
    output:
        png = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.png",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
    shell:
        "python ../src/ground_truth_clonal_tree.py  -t {input.tree} "
        "-g {input.genotypes} -T {output.gt} --draw {output.png} "

rule make_phi:
    input:
        phi ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
    output:
        cellassign = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl"
    shell:
        "python ../src/cell_mapping.py -p {input.phi} -t {input.gt} -P {output.cellassign} "
    
rule make_data:
    input: 
        sparse ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cells.p0",
    output: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    shell:
        "python ../src/data.py -f {input.sparse} -c {input.copy_profiles} -D {output} " 



# rule eval_snv_assignment:
#     input:
#         tree = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
#         data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
#     output:"input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/snv_assignment.csv"
#     shell:
#         "python ../src/snv_assignment2.py -d {input.data} -t {input.tree} -o {output} "

rule seg_tree_inference:
    input:
        tree = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        phi = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
    # params:
    #     cores = 4,
    #     lamb = 3
    output:"scores/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv"
    benchmark:"scores/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "scores/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/seg_inf.log",
        err= "scores/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/seg_inf.err.log"
    shell:
        "python ../src/segment_main2.0.py -d {input.data} -t {input.tree}  -c {input.phi} "
        "-S {output} > {log.std} 2> {log.err} "

rule get_gt_dcfs:
    input:
        tree = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
        phi = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",
        # data = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    output:
        dcfs = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/mut_cluster_tree.txt",
    shell:
        "python ../src/get_dcfs.py -t {input.tree} -c {input.phi} -d {output.dcfs} -T {output.tm}"

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





# rule sens_analysis: 
#     input:
#         data ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/data.pickle",
#         ct ="input/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}/clonal_tree.pickle",
#     params:
#         nreps = 10,
#         cell_rates = "10:110:10",
#         mut_rates = "10:110:10",
#     output: "obj_comp/s{s}_n{cells}_m{snvs}_k{nsegs}_c{cov}_l{mclust}_obj.csv"
#     shell:
#         "nice python ../src/sens_analysis.py -D {input.data}  "
#         "-T {input.ct} -c {params.cell_rates} -m {params.mut_rates} "
#         "-s {wildcards.s} -r {params.nreps} "
#         "-o {output} "







        





