
configfile: "simulate.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")
import pandas as pd 


rule all:
    # noinspection PyInterpreter
    input:
        expand("phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"],
            dirch  =config["dirch"]
        ),
        # expand("phertilizer/{inpath}/aggregate_scores.csv",
        #        inpath = config["inpath"]
        # )


    
rule make_data:
    input: "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    output: 
        readcounts = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    shell:
        "python ../src/prep_phertilizer.py -d {input} -f {output.readcounts} -u {output.umap}"


rule phertilizer:
    input:
        readcounts = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    output:
        pickle = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.pkl",
        png = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.png",
        like = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/likelihood.txt",
        phi = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.csv",
        psi = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/psi.csv",
        tree = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf_tree.txt",
    benchmark:"phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf.log",
        err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf.err.log"
    shell:
        "nice -n 10 phertilizer -f {input.readcounts} --bin_count_data {input.umap} --no-umap "
        "--tree_pickle {output.pickle} --tree {output.png}  --likelihood {output.like} "
        "-n {output.phi} -m {output.psi} --tree_text {output.tree}  --post_process -s {wildcards.s} "
        " > {log.std} 2> {log.err} "


rule convert_to_pharming_obj:
    input:
        data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
        phi = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.csv",
        psi = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/psi.csv",
        tree = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/inf_tree.txt",
    params:
        lamb = 1e3
    output:
        solution = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solution.pkl",
    log:
        std= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/convert.log",
        err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/convert.err.log"
    shell:
        "python ../src/convert_phert_to_pharm.py "
        "-d {input.data} "
        "-t {input.tree} "
        "-n {input.phi} "
        "-m {input.psi} "
        "-l {params.lamb} "
        "-o {output.solution} > {log.std} 2> {log.err} "



rule write_files:
    input: "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/solution.pkl",
    output:
        pred_mut= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        pred_cell =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        pred_genos = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genotypes.csv",
        pred_tree = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt",
        likelihood ="phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/likelihood.csv"
    run:
        import utils 
        import pandas as pd 
        sol= utils.load_pickled_object(input[0])
        sol.phi.write_phi(output["pred_cell"])
        sol.ct.write_psi(output["pred_mut"])
        # sol.ct.write_loss_mapping(output['pred_mut_loss'])
        sol.ct.write_genotypes(output["pred_genos"])
        sol.ct.save_text(output['pred_tree'])
        sol.ct.write_likelihood(output['likelihood'])


rule score_tree:
    input:
        gt_phi ="sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/cellclust_gt.csv",
        gt_mut= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/mutclust_gt.csv",
        gt_tree ="sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt",
        gt_genos ="sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/genotypes.csv",
        pred_mut= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        pred_cell =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        pred_genos = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/pred_genotypes.csv",
        pred_tree = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.txt"
    log:
        err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.err.log"
    output:"phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/metrics.csv"
    shell:
        "timeout 20m ./cpp/metrics {input.gt_tree} {input.gt_phi} {input.gt_mut} {input.gt_genos} "
        " {input.pred_tree} {input.pred_cell} {input.pred_mut} {input.pred_genos} > {output} 2> {log.err} "


# rule eval_solutions:
#     input:
#         data= "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
#         gt = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/gt.pkl",
#         cellassign = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/phi.pkl",
#         phert =   "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/tree.pkl",
#     params:
#         lamb = 0
#     output:
#         scores = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
#     log:
#         std= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.log",
#         err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/score.err.log"
#     shell:
#         "nice -n 5 python ../src/score_phertilizer.py -d {input.data} -t {input.gt} -c {input.cellassign} "
#         "-T {input.phert} "
#         "-l {params.lamb} "
#         "-o {output.scores} > {log.std} 2> {log.err} "

# rule aggregate_scores:
#     input:
#         # paths = [path for path in glob.glob("phertilizer/*/s*_m*_k*_l*/n*_c*_e*/")]
#         paths = expand("phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/scores.csv",
#             inpath = config["inpath"],
#             s =seeds,
#             cells = config["cells"],
#             snvs = config["snvs"],
#             nsegs = config["nsegs"],
#             cov = config["cov"],
#             mclust = config['mclust'],
#             err = config["cerror"]
#         ),
#     output:
#         "phertilizer/{inpath}/aggregate_scores.csv"
#     run:
#         # Read each CSV file, skipping missing files
#         dfs = []
#         import pandas as pd
#         for path in input.paths:

#             file = path
#             if os.path.exists(file):
#                 df = pd.read_csv(file)
#                 df["folder"] = path.split("/")[2]
#                 df["instance"] = path.split("/")[3]
           
#                 dfs.append(df)
#             else:
#                 print(f"Warning: Input file {file} is missing. Skipping...")
        
#         # Concatenate dataframes and save as CSV
#         if dfs:
#             result = pd.concat(dfs, ignore_index=True)
#             result.to_csv(output[0], index=False)
#         else:
#             print("No input files found. Skipping aggregation.")





        






