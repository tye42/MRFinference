# MRFinference

[Project site](https://tye42.github.io/2018/03/24/Comparison-between-Mean-Field-and-Loopy-Belief-Propagation.html)

C++ library for variational inference on Bayesian network and Markov random field. Implemented Mean Field and Loopy Belief Propagation algorithms. 

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