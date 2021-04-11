# PEPN-GRN
Petri net-based approach for the inference of gene regulatory networks from noisy gene expression data.

The PEPN-GRN repo contains code of three variants of the PEPN-GRN method and 3-bin discretized five 10-gene and five 100-gene DREAM4 time series data sets. The PEPN-GRN folder contains a "PEPN-GRN:tutorial" pdf file with a step-by-step tutorial showing how to run the code. Predicted edges in a ranked order, and AUROC and AUPR plots are generated as result.

**Inside PEPN-GRN folder:**
'sourcefiles' folder contains all the code files required for _**PEPN-GRN_v1 and PEPN-GRN_v2 implementation**_.

'multi-bin-disc-dream4-data-repository' contains EFD, EWD and Kmeans discretized data sets of five 10-gene and five 100-gene DREAM4 time series data sets. It also contains the positive and negative edges in the groundtruth of five 10-gene and five 100-gene networks.

'run_pepn-grn.sh' is the shell script to run the PEPN-GRN_v1 and PEPN-GRN_v2 code for discretized 10-gene and 100-gene DREAM4 data sets. To run PEPN-GRN code on other time series data sets, 'run_pepn-grn.sh' needs to be customized accordingly.

_**FOR PEPN-GRN_v3 implementation:**_
For PEPN-GRN_v3 implementation, first run PEPN-GRN_v1 to obtain edge probabilities. The edge probabilities generated are used as edge features in the logistic regression part of PEPN-GRN_v3. 
'log_regression' folder contains another 'sourcefiles' folder containing code files for logistic regression implementation of PEPN-GRN_v3 method. Edge probabilities from PEPN-GRN_v1 are stored in the 'input-files' folder under 'log_regression' folder.
'run_log-regression.sh' is the shell script to run the logistic regression part of the PEPN-GRN_v3 on discretized 10-gene and 100-gene DREAM4 data sets.


Note: All the source code is written in MATLAB (using MATLAB R2015a). A shell script is used to run the code.

PEPN-GRN is based on Petri nets. All the code is written by Vatsa, D.
