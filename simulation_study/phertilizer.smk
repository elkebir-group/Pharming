
configfile: "test.yml"
seeds = [i+10 for i in range(config["nseeds"])]
# seeds = [11]
import sys 
sys.path.append("../src")
import pandas as pd 



rule all:
    # noinspection PyInterpreter
    input:
        expand("phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
        expand("phertilizer/{inpath}/aggregate_scores.csv",
               inpath = config["inpath"]
        )


    
rule make_data:
    input: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    output: 
        readcounts = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    shell:
        "python ../src/prep_phertilizer.py -d {input} -f {output.readcounts} -u {output.umap}"


rule phertilizer:
    input:
        readcounts = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/read_counts.tsv",
        umap =  "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cn_umap.csv"
    output:
        pickle = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.pkl",
        tree = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.txt",
        predcell = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_cell.csv",
        predmut = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/pred_mut.csv",
        png = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.png",
        like = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/likelihood.txt",
    benchmark:"phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/benchmark.log"
    log:
        std= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.log",
        err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/inf.err.log"
    shell:
        "nice -n 10 phertilizer -f {input.readcounts} --bin_count_data {input.umap} --no-umap "
        "--tree_pickle {output.pickle} --tree {output.png} --tree_text {output.tree} --likelihood {output.like} "
        " -n {output.predcell} -m {output.predmut}  --post_process -s {wildcards.s} "
        " > {log.std} 2> {log.err} "


    
rule eval_solutions:
    input:
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        gt = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt.pkl",
        cellassign = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/phi.pkl",
        phert =   "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/tree.pkl",
    params:
        lamb = 0
    output:
        scores = "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
    log:
        std= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.log",
        err= "phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/score.err.log"
    shell:
        "nice -n 10 python ../src/score_phertilizer.py -d {input.data} -t {input.gt} -c {input.cellassign} "
        "-T {input.phert} "
        "-l {params.lamb} "
        "-o {output.scores} > {log.std} 2> {log.err} "



rule aggregate_scores:
    input:
        # paths = [path for path in glob.glob("phertilizer/*/s*_m*_k*_l*/n*_c*_e*/")]
        paths = expand("phertilizer/{inpath}/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/scores.csv",
            inpath = config["inpath"],
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
    output:
        "phertilizer/{inpath}/aggregate_scores.csv"
    run:
        # Read each CSV file, skipping missing files
        dfs = []
        import pandas as pd
        for path in input.paths:

            file = path
            if os.path.exists(file):
                df = pd.read_csv(file)
                df["folder"] = path.split("/")[2]
                df["instance"] = path.split("/")[3]
           
                dfs.append(df)
            else:
                print(f"Warning: Input file {file} is missing. Skipping...")
        
        # Concatenate dataframes and save as CSV
        if dfs:
            result = pd.concat(dfs, ignore_index=True)
            result.to_csv(output[0], index=False)
        else:
            print("No input files found. Skipping aggregation.")





        






