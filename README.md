## Requirements
You need installed **Python 3** on your system and have bash alias `python`. 
Also some libraries are required: `numpy, scipy, networkx, pygraphviz`.

## About
This repository contains the implementation of the INFER model and supplementary code.

### Breakpoint Graphs
All graphs are presented in `graphs` package and have methods for counting required statistics and drawing itself such as:
- `p_even()` — number of even paths in graph;
- `p_odd()` — number of odd paths in graph;
- `p_m(m)` — number of `m`-paths in graph;
- `chr()` — estimated number of chromosomes in graph;
- `c()` — number of cycles paths in graph;
- `c_m(m)` — number of `m`-cycles paths in graph;
- `n()` — number of all blocks in graph;
- `d()` — the minimal DCJ distance between;
- `b()` — is half the number of breakpoints;
- `save_pygraphviz(filename)` — save graph to `filename` file, 
you can specify `.svg` exension of file for getting an image or do not specify for getting file in *pygraphviz* format.

#### Cyclic genomes
Cyclic genome graphs are presented with `CyclicGenomeGraph` class in file `cyclic_genome_graph.py` file 
and can be used for simulate genome rearrangements with `do_k2_break()` method. 
You need to pass `n` to constructor 
and you can also pass distribution from `scipy.stats` library and it's params to constructor method for specifying fragalities distribution.

Examples:
```python
from src.graphs.cyclic_genome_graph import CyclicGenomeGraph

g = CyclicGenomeGraph(n=1000) # Cyclic genome graph with `1000` fragile regions in initial state
g.do_k2_break()
```

```python
from src.graphs.cyclic_genome_graph import CyclicGenomeGraph

g = CyclicGenomeGraph(n=1000, distribution="gamma", params=[1/3, 0, 1]) # Cyclic genome graph with `1000` fragile regions in initial state with gamma distribution on edges
g.do_k2_break()
```

#### Linear genomes
Linear genome graphs are presented with `LinearGenomeGraph` class in file `linear_genome_graph.py` file.
It's completely the same as a `CyclicGenomeGraph` in usage, but also need to specify number of chromosomes parameter in constructor.

Examples:
```python
from src.graphs.linear_genome_graph import LinearGenomeGraph

g = LinearGenomeGraph(n=1000, chrs=10) # Cyclic genome graph with `1000` fragile regions and `10` chromosomes in initial state
g.do_k2_break()
```

```python
from src.graphs.cyclic_genome_graph import CyclicGenomeGraph

g = CyclicGenomeGraph(n=1000, chrs=10, distribution="gamma", params=[1/3, 0, 1]) # Cyclic genome graph with `1000` fragile regions and `10` chromosomes in initial state with gamma distribution on edges
g.do_k2_break()
```

#### Real data genomes
Real genomes can be analysed by `RealDataGraph` class from file `real_data_graph.py` file.
Also this class is wrapped with bash script and can be used without writing source code.
So the application of real data will be described below.


### Estimators
All estimators are presented in `estimators` package and have methods: 
- `predict_k(g)` — pass graph from `graphs` package for prediction of `k` parameter;
- `predict_n(g)` — pass graph from `graphs` package for prediction of `n` parameter;
- `predict(g)` — for predicting both previous. 

#### Uniform Estimator
Uniform estimator is presented with with `UniformDBEstimator` class in `uniform_db_estimator.py` file. 
This estimator assumes that probabilities are uniform distributed.

#### Flat Dirichlet Estimator
Flat dirichlet estimator is presented with with `FlatDirichletDBEstimator` class in `flat_dirichlet_estimator.py` file. 
This estimator assumes that probabilities are distributed with flat dirichlet distribution.

#### Original Flat Dirichlet Estimator
Implemented flat dirichlet estimator from [original paper](https://academic.oup.com/gbe/article/8/5/1427/2939585) is presented with with `TannierEstimator` class in `tannier_dbc2_estimator.py` file
and also requires to pass `c(2)` parameter to `predict*()` methods. 
This estimator assumes that probabilities are distributed with flat dirichlet distribution.

#### Non-flat Dirichlet Estimator
Non-flat dirichlet estimator is presented with with `DirichletDBEstimator` class in `dirichlet_db_estimator.py` file. 
This estimator assumes that probabilities are distributed with dirichlet distribution with parameter `alpha`.
You need to **pass `alpha` to constructor**.

#### Corrected Non-flat Dirichlet Estimator
Corrected non-flat dirichlet estimator is presented with with `CorrectedDirichletDBEstimator` class in `dirichlet_db_estimator.py` file. 
This estimator assumes that probabilities are distributed with dirichlet distribution with parameter `alpha` and also applies the correction for linear case.
You need to **pass `alpha` to constructor**.