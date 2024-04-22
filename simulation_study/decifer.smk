configfile: "simulate.yml"
seeds = [i+10 for i in range(config["nseeds"])]
import sys 
sys.path.append("../src")


rule all:
    # noinspection PyInterpreter
    input:
        expand("decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/post_dcfs.txt",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"],
            dirch = config["dirch"]
        ),

rule prep_decifer:
    input: "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
    output:"decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/input.tsv"
    shell:
        "python ../src/prep_decifer.py -d {input} -o {output} "


rule run_decifer:
    input: 
        vaf= "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/input.tsv",
        purity= "purity.tsv",
    output: "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/decifer_output.tsv"
    params: 
        restarts= 50,
        mink = lambda wildcards: int(wildcards.mclust) - 1,
        maxk = lambda wildcards: int(wildcards.mclust)  + 1,
        prefix ="decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/decifer"
    threads: 1
    log: 
        std = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/run.log",
        err= "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/err.log",
    benchmark:"decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/benchmark.log",
    shell: 
     " decifer {input.vaf} -p {input.purity} -k {params.mink} -K {params.maxk} -j {threads} "
     " -r {params.restarts}  --seed {wildcards.s} -o {params.prefix} > {log.std} 2> {log.err} "

rule get_dcfs:
    input: "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/decifer_output.tsv"
    output: "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/dcfs.txt"
    run:
        import pandas as pd 
        dec_out = pd.read_table(input[0])
        dec_clust = dec_out[['mut_index', 'cluster', 'point_estimate_DCF0', 'true_cluster_DCF0']]
        dec_clust[['clust_dcf', 'CI']] = dec_clust['true_cluster_DCF0'].str.split(';', expand=True)
        dec_clust['clust_dcf'] = pd.to_numeric(dec_clust['clust_dcf'])
        dcfs = dec_clust['clust_dcf'].unique()
        with open(output[0], 'w') as file:
            for value in dcfs:
                file.write(str(value) + '\n')
    
rule post_process:
        input: 
            dec = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/decifer_output.tsv",
            dat = "sims/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/data.pkl",
        output: "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/post_dcfs.txt",
        params:
            vaf = 0.1,
            snv = 0.05,
            cell_thresh = 0.05
        log: 
            std = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/post.log",
            err= "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}/post.err.log",
        shell:
            "python ../src/decifer_post_process.py -d {input.dat} -f {input.dec} "
            "--vaf-thresh {params.vaf} --snv-thresh {params.snv} --cell-thresh {params.cell_thresh} "
            "-o {output}  > {log.std} 2> {log.err} "
        