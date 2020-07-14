# Codaloss: direct interaction network inference for compositional data.
=====

Description
-----
One essential step in the analysis of microbiome compositional data is inferring the direct interaction network among microbial species. Here, we propose a novel loss function called codaloss for direct microbes interaction network estimation under the sparsity assumptions. We develop an alternating direction optimization algorithm to obtain sparse solution of codaloss as estimator. Compared to other state-of-the-art methods, our model makes less assumptions about the microbial networks and outperforms them in network inference.

Requirement
------
gtools
MASS
dplyr
huge (1.2.7 version)
Rcpp
ROCR
clime
...

R files in codaloss
-----
basical_coda.R -> solutions for some basical optimization ploblems.
convenientfunc.R -> some functions for checking parameters, evaluation and ...
updatefunc.R -> updation functions for solving codaloss.
codalosscv.R -> cross validation for codaloss.
codaloss.R -> main function.
symmatrix.cpp, update.cpp -> some heavy matrix calculation.

Example
-----


