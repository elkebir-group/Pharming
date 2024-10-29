# Pharming
Solves the joint clonal tree inference problem and infers a joint CNA and SNV tree from low pass single-cell DNA sequencing.

![Overview](Pharming.png)


## Dependencies
- `networkx`
- `pandas`
- `numpy`
- `pybind11`
- `gurobipy`
- [`clonelib`](https://github.com/elkebir-group/clonesim)

*Pharming is dependent on gurobipy. Follow these [instructions](https://www.gurobi.com/features/academic-named-user-license/) to obtain a free license for academic use.*

*clonelib is packaged with [clonesim](https://github.com/elkebir-group/clonesim). Follow the clonesim instructions to install the python clonelin library.*



## Input
Pharming requires two input files:
1. alternate and total read counts for each SNV
2. allele-specific copy number profiles

## Output
Pharming has two main outputs:
1. A pickled `Solution` object containing the clonal tree, cost and cell assignment. The `Solution` object also functions for analyzing and visualizing the output
2. Visualization of the `Solution` object


## Usage
```
$ python src/main.py --help
usage: main.py [-h] [-d DATA] [-f FILE] [-c COPY_NUMBERS]
               [-s SEED] [-j CORES] [-l LAMB] [-n TOP_N] [-T TM]
               [-k SNV_CLUSTERS] [-D DCFS]
               [--delta DELTA [DELTA ...]]
               [--ninit-segs NINIT_SEGS] [--ninit-tm NINIT_TM]
               [--thresh-prop THRESH_PROP]
               [--order {random,weighted-random,nsnvs,in-place,cost}]
               [--root_x ROOT_X] [--root_y ROOT_Y] [--collapse]
               [--sum-condition]
               [--cell-threshold CELL_THRESHOLD]
               [-L SEGMENTS [SEGMENTS ...]]
               [--excl-segments EXCL_SEGMENTS [EXCL_SEGMENTS ...]]
               [--segfile SEGFILE] [-P PICKLE] [-O OUT]
               [--all-sol ALL_SOL]
               [--model-selection MODEL_SELECTION] [--tree TREE]
               [--labels LABELS]

optional arguments:
  -h, --help            show this help message and exit
  -d DATA, --data DATA  input file of preprocessed data pickle
  -f FILE, --file FILE  input file for variant and total read
                        counts with unlabled columns: [chr
                        segment snv cell var total]
  -c COPY_NUMBERS, --copy-numbers COPY_NUMBERS
                        input files of copy numbers by segment
                        with unlabeled columns [segment cell
                        totalCN]
  -s SEED, --seed SEED  random number seed (default: 1026)
  -j CORES, --cores CORES
                        Max number of cores to use for inferring
                        segment trees
  -l LAMB, --lamb LAMB  lambda value, default=1e3
  -n TOP_N, --top_n TOP_N
                        number of trees to retain in each step
  -T TM, --Tm TM        optional filename of mutation cluster
                        tree
  -k SNV_CLUSTERS, --snv-clusters SNV_CLUSTERS
                        number of SNV clusters, if dcfs are also
                        specified, k defaults to the number of
                        specified DCFs
  -D DCFS, --dcfs DCFS  optional filename of dcfs to use
  --delta DELTA [DELTA ...]
                        list of DCFs to use, ignored if dcf file
                        is provided
  --ninit-segs NINIT_SEGS
                        number of segments for initialization of
                        mutation cluster tree
  --ninit-tm NINIT_TM   number of mutation cluster trees to
                        consider after pruning with initial segs
  --thresh-prop THRESH_PROP
                        proportion threshold for determining CN
                        states
  --order {random,weighted-random,nsnvs,in-place,cost}
                        ordering strategy for progressive
                        integration, choose one of 'random',
                        'weighted-random', 'nsnvs', 'in-place'
  --root_x ROOT_X       starting state for maternal (x) allele
  --root_y ROOT_Y       starting state for paternal (y) allele
  --collapse            whether linear chains of copy number
                        events should be collapsed prior to
                        integration
  --sum-condition       use the sum condition to filter mutation
                        cluster trees
  --cell-threshold CELL_THRESHOLD
                        if collapsing is used, the minimum
                        number of cells a CNA only clone
                        requires to avoid collapsing, NA if not
                        collapsing.
  -L SEGMENTS [SEGMENTS ...], --segments SEGMENTS [SEGMENTS ...]
                        segment ids of trees to build
  --excl-segments EXCL_SEGMENTS [EXCL_SEGMENTS ...]
                        segment ids to exclude
  --segfile SEGFILE     filename with list of segments
  -P PICKLE, --pickle PICKLE
                        directory where the pickled solution
                        list of top n trees should be saved
  -O OUT, --out OUT     directory where output files should be
                        written
  --all-sol ALL_SOL     filename of object to pickle all top
                        clonal trees inferred from each mutation
                        cluster tree
  --model-selection MODEL_SELECTION
                        filename to write model selection
                        information
  --tree TREE           filename to draw a pretty clonal tree
  --labels LABELS       filename to save the encoding for the
                        labels of the pretty tree.

```


