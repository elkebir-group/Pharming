configfile: "config.yml"
seeds = [i+10 for i in range(config["nseeds"])]
import sys 
sys.path.append("../src")


rule all:
    # noinspection PyInterpreter
    input:
        expand("decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/decifer_output.tsv",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),
        expand("decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt_delta.csv",
            s =seeds,
            cells = config["cells"],
            snvs = config["snvs"],
            nsegs = config["nsegs"],
            cov = config["cov"],
            mclust = config['mclust'],
            err = config["cerror"]
        ),


        

        

        


rule prep_decifer:
    input: "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
    output:"decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/input.tsv"
    shell:
        "python ../src/prep_decifer.py -d {input} -o {output} "


rule run_decifer:
    input: 
        vaf= "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/input.tsv",
        purity= "purity.tsv",
    output: "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/decifer_output.tsv"
    params: 
        restarts= 50,
        prefix ="decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/decifer"
    log: 
        std = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/run.log",
        err= "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/err.log",
    shell: 
     " decifer {input.vaf} -p {input.purity} -k {wildcards.mclust} -K {wildcards.mclust} "
     " -r {params.restarts} --seed {wildcards.s} -o {params.prefix} > {log.std} 2> {log.err} "
    

rule construct_gt:
    input:
        tree = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/gt.pkl",
        data= "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/data.pkl",
        phi = "input/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/cellAssign.pkl",

    output: 
        psi = "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt_psi.csv",
        delta=  "decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/gt_delta.csv",
    shell: 
     " python ../src/get_gt.py -d {input.data} -t {input.tree} -c {input.phi} "
     " -p {output.psi} -f {output.delta} "
