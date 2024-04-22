configfile: "bysample.yml"
# include --use-conda when running pipeline to activate conda r_env for rule prep_input
rule all:
   input:
        expand("{sample}/decifer/dcfs.txt",
                sample = config["sample"]
            ),
        expand("{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/solutions.pkl",
               sample = config["sample"],
               clust = ["decifer"],
               isegs = config["ninit_segs"],
               tm = config["ninit_tm"],
               topn = config["topn"],
               lamb = config["lamb"]
        )
            # expand("dcf_clustering/{sample}/k{k}/dcfs.txt",
            #     sample = config["sample"],
            #     k= range(config['mink'], config['maxk']+1)
            # )

rule prep_input:
    input: 
        copy_numbers = "DeepCopyPrediction.csv",
        var_reads = "variant_data.filt.tsv"
    output:
        copy_numbers = "{sample}/input/copy_number_profiles.csv",
        var_reads = "{sample}/input/read_counts.csv",
        sampled_segs  = "{sample}/input/selected_segments.csv",
    conda: "r_env"
    params: 
        figpath = "./{sample}/input",
        seed = config["seed"],
        nsegs = config["nsegs"]
    script:
        "scripts/prep.R"




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
    output:"{sample}/decifer/input.tsv"
    params: 
        thresh = config["thresh"]
    shell:
        "python ../src/prep_decifer.py -d {input} -o {output} --cn-thresh {params.thresh} "


rule run_decifer:
    input: 
        vaf= "{sample}/decifer/input.tsv",
        purity= "purity.tsv",
    output: "{sample}/decifer/decifer_output.tsv"
    params: 
        restarts= 50,
        mink = config["mink"],
        maxk = config["maxk"],
        seed = config["seed"],
        prefix ="{sample}/decifer/decifer"
    threads: 3
    log: 
        std = "{sample}/decifer/run.log",
        err= "{sample}/decifer/err.log",
    benchmark:"{sample}/decifer/benchmark.log",
    shell: 
     "decifer {input.vaf} -p {input.purity} -k {params.mink} -K {params.maxk} -j {threads} "
     " -r {params.restarts}  --seed {params.seed} -o {params.prefix} > {log.std} 2> {log.err} "


rule get_dcfs:
    input: "{sample}/decifer/decifer_output.tsv"
    output: "{sample}/decifer/dcfs.txt"
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
    
rule pharming:
    input:
        dcfs = "{sample}/{clust}/dcfs.txt",
        data= "{sample}/input/data.pkl",
    params:
        cell_thresh = 5,
        root_x = 1,
        root_y = 1,
        max_loops = config["max_loops"],
        outdir = "./{sample}/pharming/{clust}",
        seed = config["seed"],
        thresh = config["thresh"]
    threads: 1
    output:
        sol = "{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/solutions.pkl",
        profile = "{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/run.prof",
    benchmark:"{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/benchmark.log"
    log:
        std= "{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pharm.log",
        err= "{sample}/pharming/{clust}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pharm.err.log"
    shell:
        "nice -n 10 python ../src/main.py -d {input.data} -D {input.dcfs} "
        "-s {params.seed} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--thresh-prop {params.thresh} "
        "--cell-threshold {params.cell_thresh} "
        "-O {params.outdir} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--ntree-iter {params.max_loops} "
        "--profile {output.profile} "
        "-P {output.sol} > {log.std} 2> {log.err} "

# rule run_dcf_clustering_constrained:
#     input: 
#         data=  "{sample}/input/data.pkl",
#     output: 
#         output_data= "dcf_clustering/{sample}/k{k}/results.pkl",
#         dcfs= "dcf_clustering/{sample}/k{k}/dcfs.txt"
#     params: 
#         restarts=21,
#         seed = 21
#     threads: 3
#     log: 
#         std= "dcf_clustering/{sample}/k{k}/run.log",
#         err= "dcf_clustering/{sample}/k{k}/err.log"
#     benchmark: "dcf_clustering/{sample}/k{k}/benchmark.log"
#     shell: 
#         "python ../src/dcf_clustering_v2.py -d {input.data} -k {wildcards.k} -o {output.output_data} -D {output.dcfs} "
#          "-r {params.restarts} -c -j {threads} -s {params.seed} > {log.std} 2> {log.err} "    