configfile: "bysample.yml"

rule all:
   input:
        expand("decifer/{sample}/dcfs.txt",
                sample = config["sample"]
            )

rule make_data:
    input:
        read_counts = "{sample}/input/read_counts.csv",
        copy_numbers = "{sample}/input/copy_number_profiles.csv"
    output: 
        dat = "{sample}/input/data.pkl",
        cell_lookup = "{sample}/input/cell_lookup.csv",
        mut_lookup = "{sample}/input/mut_lookup.csv"
    log: 
        std=  "{sample}/input/make_data.log",
        err= "{sample}/input/make_data.err.log"
    shell:
     "python ../src/make_data.py -f {input.read_counts} -c {input.copy_numbers} "
     " -D {output.dat} -m {output.mut_lookup} -l {output.cell_lookup} > {log.std} 2> {log.err}"


rule prep_decifer:
    input: "{sample}/input/data.pkl",
    output:"decifer/{sample}/input.tsv"
    shell:
        "python ../src/prep_decifer.py -d {input} -o {output} "


rule run_decifer:
    input: 
        vaf= "decifer/{sample}/input.tsv",
        purity= "purity.tsv",
    output: "decifer/{sample}/decifer_output.tsv"
    params: 
        restarts= 50,
        mink = config["mink"],
        maxk = config["maxk"],
        seed = 21,
        prefix ="decifer/{sample}/decifer"
    threads: 3
    log: 
        std = "decifer/{sample}/run.log",
        err= "decifer/{sample}/err.log",
    benchmark:"decifer/{sample}/benchmark.log",
    shell: 
     "decifer {input.vaf} -p {input.purity} -k {params.mink} -K {params.maxk} -j {threads} "
     " -r {params.restarts}  --seed {params.seed} -o {params.prefix} > {log.std} 2> {log.err} "


rule get_dcfs:
    input: "decifer/{sample}/decifer_output.tsv"
    output: "decifer/{sample}/dcfs.txt"
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
    
