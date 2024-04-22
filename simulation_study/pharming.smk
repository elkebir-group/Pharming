configfile: "simulate.yml"
configfile: "pharming.yml"

seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")


rule all:
    input:
        expand("pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
            inpath = ["sims"],
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
        expand("pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
            clust =["dcf_clustering", "decifer"],
            prefix = ["dcfs", "post_dcfs"],
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
        
        # expand("pharming/{inpath}/aggregate_scores.csv",
        #              inpath = config["inpath"],
        # )

rule pharming_gt_dcfs:
    input:
        dcfs = "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt",
        data= "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    params:
        cell_thresh = 5,
        root_x = 1,
        root_y = 1,
    threads: 1
    output:
        sol = "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
        profile = "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
    benchmark:"pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
        err= "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
    shell:
        "nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
        "-s {wildcards.s} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--cell-threshold {params.cell_thresh} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--profile {output.profile} "
        "-P {output.sol} > {log.std} 2> {log.err} "




rule pharming:
    input:
        dcfs = "{clust}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/{prefix}.txt",
        data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    params:
        cell_thresh = 5,
        root_x = 1,
        root_y = 1,
        max_loops = config["max_loops"],
        thresh_prop = config["thresh_prop"]
    threads: 1
    output:
        sol = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
        profile = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.prof",
    benchmark:"pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.log",
        err= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pharm.err.log"
    shell:
        "nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
        "-s {wildcards.s} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--cell-threshold {params.cell_thresh} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--thresh-prop {params.thresh_prop} "
        "--ntree-iter {params.max_loops} "
        "--profile {output.profile} "
        "-P {output.sol} > {log.std} 2> {log.err} "

rule eval_solutions:
    input:
        data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
        gt = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
        cellassign = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl",
        sol = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
    params:
        lamb = config["lamb"],
    output:
        scores = "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
    log:
        std= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.log",
        err= "pharming/{clust}/{prefix}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.err.log"
    shell:
        "nice -n 10 python ../src/score_tree.py -d {input.data} -t {input.gt} -c {input.cellassign} "
        "-S {input.sol} "
        "-l {params.lamb} "
        "-o {output.scores} > {log.std} 2> {log.err} "

# rule write_files:
#     input: "pharming/sims/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solutions.pkl",
#     output:
#         pred_mut= "pharming/sims/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
#         pred_mut_loss=   "pharming/sims/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut_loss.csv",
#         pred_cell =  "pharming/sims/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
#         tree = "pharming/sims/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt"
#     run:
#         import utils 
#         import pandas as pd 
#         sol_list = utils.load_pickled_object(input[0])
#         sol = sol_list[0]
#         sol.phi.write_phi(output["pred_cell"], include_subcluster=True)
#         psi = sol.ct.get_psi()
#         df = pd.DataFrame(list(psi.items()), columns=['mutation', 'cluster'])
#         df["seg"] = df["mutation"].map(sol.ct.mut_to_seg)
#         df["mutation"] = df["seg"].astype("str") + "_" +  df["mutation"].astype("str")
#         df = df[["cluster", "mutation"]]
#         df.to_csv(output["pred_mut"], index=False)
#         sol.ct.write_loss_mapping(output['pred_mut_loss'])
#         sol.ct.save_text(output['tree'])




# rule score_tree:
#     input:
#         gt_phi ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellclust_gt.csv",
#         gt_mut= "{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mutclust_gt.csv",
#         gt_tree ="{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt",
#         pred_cell =  "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
#         pred_mut= "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
#         pred_tree = "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt"
#     log:
#         err= "pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.err.log"
#     output:"pharming/{inpath}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv"
#     shell:
#         "./scripts/cpp/metrics {input.gt_tree} {input.gt_phi} {input.gt_mut} {input.pred_tree} {input.pred_cell} {input.pred_mut} > {output} 2> {log.err} "

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