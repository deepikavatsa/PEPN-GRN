# PEPN-GRN
Petri net-based approach for the inference of gene regulatory networks from noisy gene expression data.

The PEPN-GRN repo contains code of three variants of the PEPN-GRN method namely PEPN-GRN_v1, PEPN-GRN_v2, and PEPN-GRN_v3. It also contains 3-bin discretized five 10-gene and five 100-gene DREAM4 time series data sets. Predicted edges in a ranked order, and AUROC and AUPR plots are generated as result.

**Inside PEPN-GRN repo:**

'_source_v1_' folder contains all the code files required for _**PEPN-GRN_v1 implementation**_.

'_source_v2_' folder contains all the code files required for _**PEPN-GRN_v2 implementation**_.

'_source_v3_' folder contains all the code files required for _**PEPN-GRN_v3 implementation**_.

'_multi-bin-disc-dream4-data-repository_' folder contains EFD, EWD and Kmeans 3-bin discretized data sets of five 10-gene and five 100-gene DREAM4 time series data sets. It also contains the groundtruth of five 10-gene and five 100-gene networks.

'_run_pepn-grn.sh_' is the shell script to run the implementation of three variants of PEPN-GRN method for discretized 10-gene and 100-gene DREAM4 data sets. To run PEPN-GRN code on other time series data sets, '_run_pepn-grn.sh_' needs to be customized accordingly.

Note: All the source code is written in MATLAB (using MATLAB R2015a). A shell script is used to run the code.

All the code is written by Deepika Vatsa
