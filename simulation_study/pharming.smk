configfile: "simulate.yml"
configfile: "dcf_clustering.yml"
configfile: "pharming.yml"

seeds = [i+10 for i in range(config["nseeds"])]

import sys 
sys.path.append("../src")


rule all:
    input:
        expand("pharming/{dcf_source}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv",
            dcf_source = ["gt"],
            order = config["order"],
            isegs = config["ninit_segs"],
            tm = config["ninit_tm"],
            topn = config["topn"],
            lamb = config["lamb"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            dirch = config["dirch"],
            err = config["cerror"]
        ),
        expand("pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/model_selection2.csv",
            k = [3,4,5,6],
            dcf_source = ["gt"],
            order = config["order"],
            isegs = config["ninit_segs"],
            tm = config["ninit_tm"],
            topn = config["topn"],
            lamb = config["lamb"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            dirch = config["dirch"],
            err = config["cerror"]
        ),
        # expand("pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv",
        #     clust =["dcf_clustering", "dcf_clustering_gtk"],  #"decifer"
        #     prefix = ["dcfs"],
        #     order = config["order"],
        #     isegs = config["ninit_segs"],
        #     tm = config["ninit_tm"],
        #     topn = config["topn"],
        #     lamb = config["lamb"],
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     dirch = config["dirch"],
        #     err = config["cerror"]
        # ),
        # expand("pharming/gt/dcfs/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv",
        #     order = config["order"],
        #     isegs = config["ninit_segs"],
        #     tm = config["ninit_tm"],
        #     topn = config["topn"],
        #     lamb = config["lamb"],
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     dirch = config["dirch"],
        #     err = config["cerror"]
        # ),
        
        # expand("pharming/{inpath}/aggregate_scores.csv",
        #              inpath = config["inpath"],
        # )

rule pharming_gt_dcfs:
    input:
        dcfs = "sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        data= "sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    params:
        cell_thresh = config["cell_threshold"], #increasing threshold 16-May
        root_x = 1,
        root_y = 1,
        max_loops = 1, #config["max_loops"],
        thresh_prop = config["thresh_prop"],
        opath = "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"
    threads: 1
    output:
        sol = "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
        profile = "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
        pred_genos = "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genos.csv",
        pred_mut= "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        pred_mut_loss=   "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut_loss.csv",
        pred_cell =  "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        pred_tree = "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_tree.txt"

    benchmark:"pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
        err= "pharming/gt/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
    shell:
        "python ../src/main.py -d {input.data} -D {input.dcfs} "
        "-s {wildcards.s} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--cell-threshold {params.cell_thresh} "
        "--order {wildcards.order} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--thresh-prop {params.thresh_prop} "
        "--ntree-iter {params.max_loops} "
        "--sum-condition "
        "--profile {output.profile} "
        "-O {params.opath} "
        "-P {output.sol} > {log.std} 2> {log.err} "


rule pharming_k:
    input:
        dcfs = "dcf_clustering/k{k}/clustsegs10_r100_nfull5/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        data= "sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    params:
        cell_thresh = config["cell_threshold"], #increasing threshold 16-May
        root_x = 1,
        root_y = 1,
        max_loops = 1, #config["max_loops"],
        thresh_prop = config["thresh_prop"],
        opath = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"
    threads: 1
    output:
        sol = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
        profile = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
        pred_genos = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genos.csv",
        pred_mut= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        pred_mut_loss=   "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut_loss.csv",
        pred_cell =  "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        pred_tree = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_tree.txt",
        model_selection= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/model_selection.csv"
    benchmark:"pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
        err= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
    shell:
        "python ../src/main.py -d {input.data} -D {input.dcfs} "
        "-s {wildcards.s} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--cell-threshold {params.cell_thresh} "
        "--order {wildcards.order} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--thresh-prop {params.thresh_prop} "
        "--ntree-iter {params.max_loops} "
        "--sum-condition "
        "--profile {output.profile} "
        "-O {params.opath} "
        "--model-selection {output.model_selection} "
        "-P {output.sol} > {log.std} 2> {log.err} "


rule score_tree:
    input:
        gt_phi ="sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellclust_gt.csv",
        gt_mut= "sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mutclust_gt.csv",
        gt_tree ="sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt",
        gt_genos ="sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/genotypes.csv",
        pred_mut= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        pred_cell =  "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        pred_genos = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genos.csv",
        pred_tree = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_tree.txt"
    log:
        err= "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.err.log"
    output:"pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv"
    shell:
        "timeout 20m ./cpp/metrics {input.gt_tree} {input.gt_phi} {input.gt_mut} {input.gt_genos} "
        " {input.pred_tree} {input.pred_cell} {input.pred_mut} {input.pred_genos} > {output} 2> {log.err} "

rule model_selection:
    input:         
        sol = "pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
        data= "sims3/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    output:"pharming/dcf_clustk{k}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/model_selection2.csv"
    run:
       import pandas as pd
       from solution import Solution 
       from clonal_tree import ClonalTree 
       dat = pd.read_pickle(input.data)
       solutions = pd.read_pickle(sol)
       icl, bic, ent = sol[0].ICL(dat, wildcards["lamb"])
       cost, snv, cna = sol[0].compute_likelihood(dat, wildcards["lamb"])
       with open(output[0], "w+") as file:
            file.write("ICL,BIC,ENT,cost,snv_cost,cna_cost\n")
            file.write(f"{icl},{bic},{ent},{cost},{snv},{cna}\n")
# rule pharming:
#     input:
#         dcfs = "{clust}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/{prefix}.txt",
#         data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#     params:
#         cell_thresh = 5,
#         root_x = 1,
#         root_y = 1,
#         max_loops = config["max_loops"],
#         thresh_prop = config["thresh_prop"]
#     threads: 1
#     output:
#         sol = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#         profile = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
#     benchmark:"pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
#     log:
#         std= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
#         err= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
#     shell:
#         "timeout 12h nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
#         "-s {wildcards.s} "
#         "-l {wildcards.lamb} "
#         "-n {wildcards.topn} "
#         "-j {threads} "
#         "--ninit-segs {wildcards.isegs} "
#         "--ninit-tm {wildcards.tm} "
#         "--cell-threshold {params.cell_thresh} "
#         "--root_x {params.root_x} --root_y {params.root_y} " 
#         "--collapse "
#         "--thresh-prop {params.thresh_prop} "
#         "--ntree-iter {params.max_loops} "
#         "--profile {output.profile} "
#         "-P {output.sol} > {log.std} 2> {log.err} "

# rule eval_solutions:
#     input:
#         data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#         gt = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
#         cellassign = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl",
#         sol = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#     params:
#         lamb = config["lamb"],
#     output:
#         scores = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
#     log:
#         std= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.log",
#         err= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.err.log"
#     shell:
#         "nice -n 10 python ../src/score_tree.py -d {input.data} -t {input.gt} -c {input.cellassign} "
#         "-S {input.sol} "
#         "-l {params.lamb} "
#         "-o {output.scores} > {log.std} 2> {log.err} "

# rule write_files:
#     input: "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#     output:
#         pred_mut= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
#         pred_mut_loss=   "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut_loss.csv",
#         pred_cell =  "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
#         pred_genos = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genotypes.csv",
#         pred_tree = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt"
#     run:
#         import utils 
#         import pandas as pd 
#         sol_list = utils.load_pickled_object(input[0])
#         sol = sol_list[0]
#         sol.phi.write_phi(output["pred_cell"])
#         sol.ct.write_psi(output["pred_mut"])
#         sol.ct.write_loss_mapping(output['pred_mut_loss'])
#         sol.ct.write_genotypes(output["pred_genos"])
#         sol.ct.save_text(output['pred_tree'])











# rule aggregate_scores:
#     input:
#         # paths = [path for path in glob.glob("phertilizer/*/s*_m*_k*_l*/n*_c*_e*/")]
#         paths = expand("pharming//{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
#          inpath = config["inpath"],
#             order = config["order"],
#             isegs = config["ninit_segs"],
#             tm = config["ninit_tm"],
#             topn = config["topn"],
#             lamb = config["lamb"],
#             s =seeds,
#             cells = config["cells"],
#             snvs = config["snvs"],
#             nsegs = config["nsegs"],
#             cov = config["cov"],
#             mclust = config['mclust'],
#             dirch = config["dirch"],
#             err = config["cerror"]
#         ),
#     output:
#         "pharming/sims/aggregate_scores.csv"
#     run:
#         # Read each CSV file, skipping missing files
#         dfs = []
#         import pandas as pd
#         for path in input.paths:

#             file = path
#             if os.path.exists(file):
#                 df = pd.read_csv(file)
#                 df["folder"] = path.split("/")[4]
#                 df["instance"] = path.split("/")[5]
#                 df["params"] = path.split("/")[3]
#                 df["prog_order"] =path.split("/")[2]
           
#                 dfs.append(df)
#             else:
#                 print(f"Warning: Input file {file} is missing. Skipping...")
        
#         # Concatenate dataframes and save as CSV
#         if dfs:
#             result = pd.concat(dfs, ignore_index=True)
#             result.to_csv(output[0], index=False)
#         else:
#             print("No input files found. Skipping aggregation.")




# rule pharming_decifer:
#     input:
#         dcfs = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
#         data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#     params:
#         cell_thresh = 5,
#         root_x = 1,
#         root_y = 1,
#         max_loops = config["max_loops"]
#     threads: 1
#     output:
#         sol = "pharming/decifer/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#         profile = "pharming/decifer/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
#     benchmark:"pharming/decifer/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
#     log:
#         std= "pharming/decifer/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
#         err= "pharming/decifer/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
#     shell:
#         "nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
#         "-s {wildcards.s} "
#         "-l {wildcards.lamb} "
#         "-n {wildcards.topn} "
#         "-j {threads} "
#         "--ninit-segs {wildcards.isegs} "
#         "--ninit-tm {wildcards.tm} "
#         "--cell-threshold {params.cell_thresh} "
#         "--root_x {params.root_x} --root_y {params.root_y} " 
#         "--collapse "
#         "--ntree-iter {params.max_loops} "
#         "--profile {output.profile} "
#         "-P {output.sol} > {log.std} 2> {log.err} "


# rule pharming_dcf_clust:
#     input:
#         dcfs = "dcf_clustering_v2/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf_dcfs.txt",
#         data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#     params:
#         cell_thresh = 5,
#         root_x = 1,
#         root_y = 1,
#         max_loops = config["max_loops"]
#     threads: 1
#     output:
#         sol = "pharming/dcf_clust/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#         profile = "pharming/dcf_clust/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
#     benchmark:"pharming/dcf_clust/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
#     log:
#         std= "pharming/dcf_clust/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
#         err= "pharming/dcf_clust/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
#     shell:
#         "nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
#         "-s {wildcards.s} "
#         "-l {wildcards.lamb} "
#         "-n {wildcards.topn} "
#         "-j {threads} "
#         "--ninit-segs {wildcards.isegs} "
#         "--ninit-tm {wildcards.tm} "
#         "--cell-threshold {params.cell_thresh} "
#         "--root_x {params.root_x} --root_y {params.root_y} " 
#         "--collapse "
#         "--ntree-iter {params.max_loops} "
#         "--profile {output.profile} "
#         "-P {output.sol} > {log.std} 2> {log.err} "