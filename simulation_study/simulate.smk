configfile: "simulate.yml"
seeds = [i+10 for i in range(config["nseeds"])]
import sys 
sys.path.append("../src")

# #TODO: add rule generate cna trees 
ruleorder:  generatesinglecells > make_data

rule all:
    input:
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/sparse.p0",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"],
            dirch = config["dirch"]
        ),
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/genotypes.csv",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            dirch = config["dirch"],
            err = config["cerror"]
        ),
        expand("{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/loss_error_report.csv",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            dirch = config["dirch"],
            err = config["cerror"]
        ),


rule simulate:
    input: "cnatrees_nocompleteloss.txt"
    output:
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/tree.tsv",
        prop=  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/proportions.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/node.tsv",
    params:
        pth = config["simpath"],
        cp_thresh =config["cp_thresh"],
        purity = 0.99,
        sample = 1,
        alpha=  config["alpha"],
        truncalSegs = config["truncalsegs"],
        threshold = config["threshold"],
        lossProb = config["lossProb"],
        simout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/run.log", 
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/err.log" 
    shell:
        "{params.pth} -r -S {input} "
        " -purity {params.purity} -minProp {params.cp_thresh} "
        " -kk {params.truncalSegs} -f "
        "-dirich_param {wildcards.dirch}   -threshold {params.threshold} "
        "-s {wildcards.s}  -l {wildcards.mclust} "
        "-lossProb {params.lossProb} "
        "-k {wildcards.nsegs} -n {wildcards.snvs} -m {params.sample} "
        "-output_file_dir {params.simout_dir}  > {log.std} 2> {log.err}  "


rule generatesinglecells:
    input:
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/tree.tsv",
        prop=  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/proportions.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/node.tsv",
    output:
        phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        sparse ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cells.p0",
    params:
        pth = config["genpath"],
        purity = 0.99,
        sample = 1,
        alpha= config["alpha"],
        simout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}",
        scout_dir = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.log", 
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/err.log" 
    shell:
     "{params.pth} -num_cells {wildcards.cells} -read_depth {wildcards.cov} -k {wildcards.nsegs}  "
     "-alpha_fp {params.alpha} -out_dir {params.scout_dir} -in_dir {params.simout_dir} -e {wildcards.err} "
     " -m {params.sample} > {log.std} 2> {log.err} "

rule make_data:
    input: 
        sparse ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/sparse.p0",
        copy_profiles ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cells.p0",
    params: 
        alpha = config["alpha"]
    output: 
        data = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
        cell_lookup = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cells_lookup.csv",
        mut_lookup = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mut_lookup.csv",
        # seg_lookup = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/segment_lookup.csv",

    shell:
        "python ../src/make_data.py -f {input.sparse} -c {input.copy_profiles} -D {output.data} "
        "-a {params.alpha} -l {output.cell_lookup} -m {output.mut_lookup} " 




rule make_gt:
    input: 
        phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellAssignments.p0",
        tree =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/tree.tsv",
        genotypes ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/node.tsv",
        data = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    params: 
        lamb = config["lamblikelihood"]
    output:
        png = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.png",
        gt = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
        phi = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl",
        dcfs ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        tm = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/Tm.txt",
    log:
        std ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.log",
        err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.err.log",
    shell:
        "python ../src/ground_truth_clonal_tree.py -d {input.data}  -t {input.tree} -p {input.phi} "
        "-g {input.genotypes} -T {output.gt} --draw {output.png}  -D {output.dcfs} -P {output.phi} "
        "-l {params.lamb} --mut-cluster-tree {output.tm} > {log.std} 2> {log.err} "


rule write_gt_files:
    input:
        gt = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
        phi = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl",
    params: config["lamblikelihood"]
    output:
        gt_mut= "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mutclust_gt.csv",
        gt_mut_loss=   "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mut_loss_clust_gt.csv",
        gt_phi =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellclust_gt.csv",
        gt_genos =  "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/genotypes.csv",
        tree = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt",
        likelihood = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/likelihood.csv"
    run:
        import utils 
        gt = utils.load_pickled_object(input['gt'])
        phi = utils.load_pickled_object(input['phi'])
        phi.write_phi(output["gt_phi"])
        gt.write_psi(output['gt_mut'])
        gt.write_loss_mapping(output['gt_mut_loss'])
        gt.write_genotypes(output["gt_genos"])
        gt.save_text(output['tree'])
        gt.write_likelihood(output['likelihood'])
        


rule filter_loss_list:
    input: 
        gt = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
        phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl"
    output: 
        filtered_snvs = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/snv_filtered_loss.csv",
        loss_err ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/loss_error_report.csv"
    run:
        import networkx as nx
        import pandas as pd 
        gt = pd.read_pickle(input.gt)
        phi = pd.read_pickle(input.phi)
        gt_loss_snvs = gt.get_lost_snvs()
        filtered_snvs = []
        good_snvs =[]
        for j in gt_loss_snvs:

            for n in gt.preorder():

               if n in gt.mut_mapping and j in gt.mut_mapping[n]:
                    for u, lsnvs in gt.mut_loss_mapping.items():
                        if j in lsnvs:
                            lost_node = u 
                            break
                    lost_path = nx.shortest_path(gt.tree, n, lost_node)
                    has_cells = False
                    for p in lost_path:
                        if p == lost_node:
                            continue
                        if len(phi.get_cells(p)) > 0:
                            has_cells = True
                            filtered_snvs.append([j,lost_node])
                            good_snvs.append(j)
                            break 
  
        with open(output['loss_err'], "w+") as file:
            file.write("num_lost,num_retained,num_filtered\n")
            num_filtered = len(gt_loss_snvs) - len(good_snvs)
            file.write(f"{len(gt_loss_snvs)},{len(good_snvs)},{num_filtered}\n")

        pd.DataFrame(filtered_snvs, columns=["mutation", "cluster"]).to_csv(output['filtered_snvs'], index=False)









        






