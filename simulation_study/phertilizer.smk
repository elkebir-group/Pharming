
configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")
import pandas as pd 



rule all:
    # noinspection PyInterpreter
    input:
        # expand("phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
        #     s =seeds,
        #     cells = config["cells"],
        #     snvs = config["snvs"],
        #     nsegs = config["nsegs"],
        #     cov = config["cov"],
        #     mclust = config['mclust'],
        #     err = config["cerror"]
        # ),
        "phertilizer/aggregate_scores.csv"


    
rule make_data:
    input: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    output: 
        readcounts = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    shell:
        "python ../src/prep_phertilizer.py -d {input} -f {output.readcounts} -u {output.umap}"


# rule phertilizer:
#     input:
#         readcounts = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
#         umap =  "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
#     output:
#         pickle = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.pkl",
#         tree = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.txt",
#         predcell = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_cell.csv",
#         predmut = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_mut.csv",
#         png = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.png",
#         like = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/likelihood.txt",
#     benchmark:"phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
#     log:
#         std= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.log",
#         err= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.err.log"
#     shell:
#         "phertilizer -f {input.readcounts} --bin_count_data {input.umap} --no-umap "
#         "--tree_pickle {output.pickle} --tree {output.png} --tree_text {output.tree} --likelihood {output.like} "
#         " -n {output.predcell} -m {output.predmut}  --post_process -s {wildcards.s} "
#         " > {log.std} 2> {log.err} "


    
rule eval_solutions:
    input:
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
        cellassign = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/phi.pkl",
        phert =   "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.pkl",
    params:
        lamb = 0
    output:
        scores = "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
    log:
        std= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.log",
        err= "phertilizer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.err.log"
    shell:
        "nice -n 10 python ../src/score_phertilizer.py -d {input.data} -t {input.gt} -c {input.cellassign} "
        "-T {input.phert} "
        "-l {params.lamb} "
        "-o {output.scores} > {log.std} 2> {log.err} "



rule aggregate_scores:
    input:
        paths = [path for path in glob.glob("phertilizer/s*_m*_k*_l*/n*_c*_e*/")]
    output:
        "phertilizer/aggregate_scores.csv"
    run:
        # Read each CSV file, skipping missing files
        dfs = []
        import pandas as pd
        for path in input.paths:

            file = f"{path}/scores.csv"
            if os.path.exists(file):
                df = pd.read_csv(file)
                df["folder"] = path.split("/")[1]
                df["instance"] = path.split("/")[2]
           
                dfs.append(df)
            else:
                print(f"Warning: Input file {file} is missing. Skipping...")
        
        # Concatenate dataframes and save as CSV
        if dfs:
            result = pd.concat(dfs, ignore_index=True)
            result.to_csv(output[0], index=False)
        else:
            print("No input files found. Skipping aggregation.")





        






