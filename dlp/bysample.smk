configfile: "bysample.yml"
import sys
import pandas as pd 
sys.path.append("../src")



# include --use-conda when running pipeline to activate conda r_env for rule prep_input
rule all:
   input:
        expand("{sample}/input/data.pkl",
            sample = config["sample"]
        ),
        expand( "{sample}/decifer/k{k}/post_dcfs.txt",
                k= range(config['mink'], config['maxk']+1),
                sample = config["sample"]
        ),
        expand( "{sample}/input/s{seed}_r{nsegs}/sampled_segments.csv",
                seed = [20 + i for i in range(config["nseeds"])],
                sample = config["sample"],
                nsegs = config["nsegs"]
        ),
        expand("{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/cmb.csv",
                sample = config["sample"],
                k= range(config['mink'], config['maxk']+1),
                 seed = [20 + i for i in range(config["nseeds"])],
                nsegs = config["nsegs"],
               clust = ["decifer"],
               isegs = config["ninit_segs"],
               tm = config["ninit_tm"],
               topn = config["topn"],
               lamb = config["lamb"]

        ),
        expand("{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pred_cell.csv",
                sample = config["sample"],
                k= range(config['mink'], config['maxk']+1),
                 seed = [20 + i for i in range(config["nseeds"])],
                nsegs = config["nsegs"],
               clust = ["decifer"],
               isegs = config["ninit_segs"],
               tm = config["ninit_tm"],
               topn = config["topn"],
               lamb = config["lamb"]

        )

rule prep_input:
    input: 
        copy_numbers = "raw/DeepCopyPrediction.csv",
        var_reads = "input/read_counts.csv"
    output:
        copy_numbers = "{sample}/input/copy_number_profiles.csv",
        var_reads = "{sample}/input/read_counts.csv",
        sampled_segs  = "{sample}/input/selected_segments.csv",
    conda: "r_env"
    params: 
        figpath = "./{sample}/input",
        thresh = config["thresh"],
        minvar = config["minvar"]
        # nsegs = config["nsegs"]
    script:
        "scripts/prep.R"




rule make_data:
    input:
        read_counts = "{sample}/input/read_counts.csv",
        copy_numbers = "{sample}/input/copy_number_profiles.csv"
    output: 
        dat = "{sample}/input/data.pkl",
        cell_lookup = "{sample}/input/cell_lookup.csv",
        mut_lookup = "{sample}/input/mut_lookup.csv",
        seg_lookup = "{sample}/input/segment_lookup.csv"
    log: 
        std=  "{sample}/input/make_data.log",
        err= "{sample}/input/make_data.err.log"
    shell:
     "python ../src/make_data.py -f {input.read_counts} -c {input.copy_numbers} "
     "-s {output.seg_lookup} "
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
    output: "{sample}/decifer/k{k}/decifer_output.tsv"
    params: 
        restarts= config["restarts"],
        # mink = config["mink"],
        # maxk = config["maxk"] +1,
        seed = config["seed"],
        prefix ="{sample}/decifer/k{k}/decifer"
    threads: 5
    log: 
        std = "{sample}/decifer/k{k}/run.log",
        err= "{sample}/decifer/k{k}/err.log",
    benchmark:"{sample}/decifer/k{k}/benchmark.log",
    shell: 
     "decifer {input.vaf} -p {input.purity} -k {wildcards.k} -K {wildcards.k} -j {threads} "
     " -r {params.restarts}  --seed {params.seed} -o {params.prefix} > {log.std} 2> {log.err} "


rule get_dcfs:
    input: "{sample}/decifer/k{k}/decifer_output.tsv"
    output: "{sample}/decifer/k{k}/dcfs.txt"
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
        dec = "{sample}/decifer/k{k}/decifer_output.tsv",
        dat =  "{sample}/input/data.pkl",
    output: "{sample}/decifer/k{k}/post_dcfs.txt"
    params:
        vaf = 0.1,
        snv = 0.05,
        cell_thresh = 0.05
    log: 
        std = "{sample}/decifer/k{k}/post.log",
        err= "{sample}/decifer/k{k}/post.err.log",
    shell:
        "python ../src/decifer_post_process.py -d {input.dat} -f {input.dec} "
        "--vaf-thresh {params.vaf} --snv-thresh {params.snv} --cell-thresh {params.cell_thresh} "
        "-o {output}  > {log.std} 2> {log.err} "
# rule dcf_clustering_k:
#     input: 
#         data= "{sample}/input/data.pkl",
#     output: 
#         output_data= "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/output.pkl",
#         dcfs= "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/dcfs.txt",
#         post_dcfs = "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/post_dcfs.txt"
#     params: 
#         thresh_prop = config["thresh"]
#     threads: 5
#     log: 
#         std= "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/run.log",
#         err= "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/err.log"
#     benchmark: "{sample}/dcf_clustering/k{k}/clustsegs{clustsegs}_r{restarts}_nfull{nfull}_maxcn{maxcn}/benchmark.log"
#     shell: 
#         "python ../src/dcf_clustering_v2.py -d {input.data}  -o {output.output_data} -D {output.dcfs} "
#         " -k {wildcards.k} --max-cn-states {wildcards.maxcn} -s {wildcards.k} -P {output.post_dcfs} "
#         "--verbose --nsegs {wildcards.clustsegs} -r {wildcards.restarts} --nfull {wildcards.nfull} --thresh-prop {params.thresh_prop} "
#         " -c -j {threads} > {log.std} 2> {log.err}"    

rule sampled_segs:
    input: "{sample}/input/data.pkl",
    output: "{sample}/input/s{seed}_r{nsegs}/sampled_segments.csv"
    params:
        thresh = config["thresh"],
        max_cn_state = config["maxcn"],
        min_snvs = config["min_snvs"]
    run:
        import pandas as pd
        import numpy as np 
        from data import Data
        rng = np.random.default_rng(int(wildcards["seed"]))

        dat = pd.read_pickle(input[0])
        print(dat)
        cand_segs = dat.sample_segments(int(wildcards["nsegs"]), rng=rng, 
                         max_cn_states=params["max_cn_state"], 
                         min_snvs =params["min_snvs"], thresh=params["thresh"])
            
    
        with open(output[0], "w+") as file:
            for ell in cand_segs:
                file.write(f"{ell}\n")
        


rule pharming:
    input:
        dcfs = "{sample}/decifer/k{k}/dcfs.txt",
        data= "{sample}/input/data.pkl",
        segfile= "{sample}/input/s{seed}_r{nsegs}/sampled_segments.csv"
    params:
        root_x = 1,
        root_y = 1,
        max_loops = config["max_loops"],
        opath = "./{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}",
        seed = config["seed"],
        thresh_prop = config["thresh"],
        order = config["order"]
    threads: 1
    output:
        sol = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/solutions.pkl",
        profile = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/run.prof",
        model_selection= "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/model_selection.csv"
    benchmark:"{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/benchmark.log"
    log:
        std= "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pharm.log",
        err= "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pharm.err.log"
    shell:
        "python ../src/main.py -d {input.data} -D {input.dcfs} "
        "--segfile {input.segfile} "
        "-s {wildcards.seed} "
        "-l {wildcards.lamb} "
        "-n {wildcards.topn} "
        "-j {threads} "
        "--ninit-segs {wildcards.isegs} "
        "--ninit-tm {wildcards.tm} "
        "--order {params.order} "
        "--root_x {params.root_x} --root_y {params.root_y} " 
        "--collapse "
        "--thresh-prop {params.thresh_prop} "
        "--sum-condition "
        "--profile {output.profile} "
        "-O {params.opath} "
        "--model-selection {output.model_selection} "
        "-P {output.sol} > {log.std} 2> {log.err} "



rule cmb:
    input:
        data= "{sample}/input/data.pkl",
        sol = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/solutions.pkl",
    params:
        min_cells = config["mincells"]
    output:  
        cmb = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/cmb.csv",
       vaf = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/vaf.csv",
    shell:
        "python ../src/cmb.py -d {input.data} -s {input.sol} -o {output.cmb} --min-cells {params.min_cells} -v {output.vafs} "


rule write_flat_files:
    input:
        data= "{sample}/input/data.pkl",
        sol = "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/solutions.pkl",
    params:
        opath= "./{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}"
    output: "{sample}/pharming/{clust}_k{k}/s{seed}_r{nsegs}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/pred_cell.csv"
    run:
        
        dat = pd.read_pickle(input['data'])
        sols = pd.read_pickle(input['sol'])
        bst = sols[0]
        bst.write_flat_files(params['opath'], dat, float(wildcards["lamb"]) )
        bst.drawPrettyTree(f"{params['opath']}/prettyTree.pdf", f"{params['opath']}/labels.csv")



    # # args = parser.parse_args([
    #     "-d", f"{pth}/input/data.pkl",
    #     "-s", f"{pth}/{pth2}/solutions.pkl",
    #     "-o", f"{pth}/{pth2}/cmb.csv"
    # ])



# rule run_dcf_clustering:
#     input: 
#         data=  "{sample}/input/data.pkl",
#     output: 
#         output_data= "{sample}/dcf_clustering/k{k}/results.pkl",
#         dcfs= "{sample}/dcf_clustering/k{k}/dcfs.txt",
#         post ="{sample}/dcf_clustering/k{k}/post_dcfs.txt"
#     params: 
#         restarts=20,
#         seed = 21
#     threads: 3
#     log: 
#         std= "{sample}/dcf_clustering/k{k}/run.log",
#         err= "{sample}/dcf_clustering/k{k}/err.log"
#     benchmark: "{sample}/dcf_clustering/k{k}/benchmark.log"
#     shell: 
#         "python ../src/dcf_clustering_v2.py -d {input.data} -k {wildcards.k} -o {output.output_data} -D {output.dcfs} "
#          "-P {output.post} -r {params.restarts} -c -j {threads} -s {params.seed} > {log.std} 2> {log.err} "    