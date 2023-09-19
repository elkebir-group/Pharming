# Pharming
Solves the clonal tree inference with copy numbers problem

Recursively clone the repository to access submodules.
```
git clone --recursive git@github.com:elkebir-group/Pharming.git
```

Create environment: 
```
cd Pharming
conda env create -f pharming.yml -n pharming
conda activate pharming
```

Install cnatrees submodule as a python package:
```
cd cnatrees
python setup.py install
```
*TODO: Fix harcoded paths for for lemon and boost*

To test cnatrees:
```
python test.py
```

Pharming is dependent on gurobi. Follow these [instructions](https://www.gurobi.com/features/academic-named-user-license/) to obtain a free license for academic use.



