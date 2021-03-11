# PGM Inference Algorithms

C++ library for variational inference on Bayesian network (BN) and Markov random field (MRF). Implemented Mean Field (MF) and Loopy Belief Propagation (LBP) algorithms. 

Usages:

```c++
network net;
// read BN/MRF graph from file
net.readfromfile(inputfile);
// approximate inference using mean field
network::mfresults mf = net.mf();
// approximate inference using loopy belief propagation
network::lbpresults lbp = net.lbp();
```

Input file format:

```
96 # the number of factors in the graph

2 # the number of variables in the factor
0 6 # variable names
2 2 # the cardinality of each variable
4 # the number of non-zero entries in the factor
0  0.0959631
1    10.4207
2    10.4207
3  0.0959631

...
```

## Comparision between MF and LBP

*Mean Field* and *Loopy Belief Propagation* are two common [variational inference](https://en.wikipedia.org/wiki/Variational_Bayesian_methods) approaches for graphical models (Bayesian network and Markov random field). In many cases, when exact inference is infeasible, variational inference uses a simplified model to run efficient inference and give approximate posterior distributions.

Mean field assumes the variational distribution $Q$ factorize over each variables: $Q(X) = \prod_i Q(x_i)$. Loopy belief propagation is a generalization of belief propagation on cluster graphs contain cycles. However, few research has been done to compare the performance of MF and LBP. 

Pairwise Markov Network is a class of MRF where all factors associated with either one node or one pair, it is commonly used as a benchmark for comparison of inference algorithms. Here I used three types of pairwise MRF structures: 2D grid, regular graph, and complete graph with binary variables.

The factors in 2D grids were generated using [Ising model](https://en.wikipedia.org/wiki/Ising_model). The energy function of Ising model has the following form: $\epsilon_i (x_i )=u_i x_i$ for unary factor, and $\epsilon_{i,j} (x_i,x_j )=w_{i,j} x_i x_j$ for pairwise factor. 

Here, each $u_i$ is uniformly draw from $[-1,1]$, and each $w_{i,j}$ is uniformly draw from $[-C,C]$, where C is a positive constant defined the interaction strength of the system. 

Additionally, the logarithm of each entry of factors in regular and complete graph is draw independently from a normal distribution with mean 0 and standard deviation $\beta$, and called ExpGauss.

The *junctionTree* algorithm from libDAI is used to calculated the exact beliefs as the "ground truth".

To study the impact of degree and interaction strength on performance of the two algorithms, test datasets were generated with the following combinations:

- 2D Ising:  $N_g$ = 15, C = 1, 2, 3, 4, 5
- Complete ExpGauss: $N_c$ = 20, $\beta$ = 1, 2, 3, 4, 5
- Regular ExpGauss: $N_r$ = 20, d = 5, 10, 15, $\beta$ = 1
- Regular Ising: $N_r$ = 20, d = 5, 10, 15, C = 1

The results of experiments can be found at [here](https://tye42.github.io/2018/03/24/Comparison-between-Mean-Field-and-Loopy-Belief-Propagation.html).
