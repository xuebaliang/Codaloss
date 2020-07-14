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

huge(1.2.7 version) 

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
We generate several simulation data to run each method as an example. Specifically, you can download all files in github and focus on the "simulation.R" file. The only thing you need to be careful about is the determination of relative addresses. For example, we set sample number n as 200 and the species number p as 50. At random seed 1, we run the R codes, and you will obtain the following results. 

AUC result

      network       n     p     glasso    clime   SPIEC-EASI(gl)   gCoda    CDtrace     CDtr    codaloss
      
[1,] "random"     "200"  "50"  "0.5658"  "0.5681"  "0.859"        "0.8466"  "0.7924"  "0.8226"  "0.9038"

[2,] "neighbor"   "200"  "50"  "0.5754"  "0.5744"  "0.7634"       "0.7763"  "0.8578"  "0.675"   "0.8639"

[3,] "band"       "200"  "50"  "0.6302"  "0.6065"  "0.7885"       "0.7283"  "0.813"   "0.7022"  "0.8283"

[4,] "hub"        "200"  "50"  "0.6073"  "0.6681"  "0.7378"       "0.713"   "0.7438"  "0.77"    "0.8212"

[5,] "block"      "200"  "50"  "0.6042"  "0.5952"  "0.7694"       "0.6639"  "0.7848"  "0.7387"  "0.8529"

[6,] "scalefree"  "200"  "50"  "0.5794"  "0.6749"  "0.7598"       "0.7078"  "0.7585"  "0.7929"  "0.8164"

Theta bias result

      network       n     p     glasso      clime    SPIEC-EASI(gl)    gCoda     CDtrace    CDtr    codaloss
      
[1,] "random"     "200"  "50"  "7.7486"  "2929.1818"  "39.5574"       "4.6485"  "4.0563"  "3.8406"  "2.772" 

[2,] "neighbor"   "200"  "50"  "9.227"   "2283.9225"  "20.4278"       "5.3223"  "3.5172"  "5.2669"  "3.2591"

[3,] "band"       "200"  "50"  "7.094"   "2011.1383"  "27.2118"       "4.8138"  "3.1996"  "3.9705"  "3.0484"

[4,] "hub"        "200"  "50"  "5.5016"  "1516.1031"  "40.9963"       "4.532"   "3.251"   "3.1306"  "2.8052"

[5,] "block"      "200"  "50"  "8.4885"  "1583.2738"  "31.654"        "4.819"   "3.1582"  "3.509"   "2.597" 

[6,] "scalefree"  "200"  "50"  "5.6554"  "1503.0662"  "31.1912"       "4.2624"  "3.0337"  "2.9165"  "2.7099"

Usage for codaloss
-----
As you can see in the line 296-332 in "simulation.R" file, codaloss need the users to provide the compositional data as the input. Codaloss would first use the cross validation procedure to obtain the optimal lasso penalty parameter "weightoptlambda" and then run the whole model again to get the final estimation of precision matrix "Theta". The "geteval" function could supply some evaluation indexes, like AUC, Thetabias, codabias and invbias. Lastly, we would make the whole codes into a R package, and then it will be released. 

Contributing
-----
Authors email: clandzyy@pku.edu.cn; heshun@pku.edu.cn.
