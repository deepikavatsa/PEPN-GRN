# PEPN-GRN
Petri net-based approach for the inference of gene regulatory networks from noisy gene expression data.

The PEPN-GRN repo contains code of three variants of the PEPN-GRN method namely PEPN-GRN_v1, PEPN-GRN_v2, and PEPN-GRN_v3. It also contains 3-bin discretized five 10-gene and five 100-gene DREAM4 time series data sets. Predicted edges in a ranked order, and AUROC and AUPR plots are generated as result.

**Inside PEPN-GRN repo:**

'_sourcefiles_' folder contains all the code files required for _**PEPN-GRN_v1 and PEPN-GRN_v2 implementation**_.

'_multi-bin-disc-dream4-data-repository_' folder contains EFD, EWD and Kmeans discretized data sets of five 10-gene and five 100-gene DREAM4 time series data sets. It also contains the positive and negative edges in the groundtruth of five 10-gene and five 100-gene networks.

'_run_pepn-grn.sh_' is the shell script to run the PEPN-GRN_v1 and PEPN-GRN_v2 code for discretized 10-gene and 100-gene DREAM4 data sets. To run PEPN-GRN code on other time series data sets, '_run_pepn-grn.sh_' needs to be customized accordingly.

_**FOR PEPN-GRN_v3 implementation:**_
For PEPN-GRN_v3 implementation, first run PEPN-GRN_v1 to obtain edge probabilities. The edge probabilities generated are used as edge features in the logistic regression part of PEPN-GRN_v3. 
'_log-regression_' folder contains another '_sourcefiles_' folder containing code files for logistic regression implementation of PEPN-GRN_v3 method. Edge probabilities from PEPN-GRN_v1 are stored in the '_input-files_' folder under '_log-regression_' folder.
'_run_log-regression.sh_' is the shell script to run the logistic regression part of the PEPN-GRN_v3 on discretized 10-gene and 100-gene DREAM4 data sets.


Note: All the source code is written in MATLAB (using MATLAB R2015a). A shell script is used to run the code.

All the code is written by Vatsa, D.
