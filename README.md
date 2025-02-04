This repository have functions related to the published paper Sofer T. Confidence intervals for heritability via Haseman-Elston regression. Stat Appl Genet Mol Biol. 2017 Sep 26;16(4):259-273. doi: 10.1515/sagmb-2016-0076. PMID: 28862991.

File TSofer_varComp_CI_functions.R contains the original functions (perhaps with some bug fixes?) from the paper, used to estimate variance components, heritability, and confidence intervals from a single study. 
File Single_study_multiple_VC.R demonstrates how the functions in the file TSofer_varComp_CI_functions.R can be used. 

File pooled_vs_stratified_one_VCs.R provide code demonstrating heritability and confidence interval estimation meta-analysis, as described in the 2017 paper. This code can be used only when focusing on a single variance component (as is the common scenario). 

File Dev_heritability_CI_when_zero_heritability.R provides a quick implementation of heritability and CI estimation when there is only one correlation matrix used to model relatedness between individual in the data (e.g. just a kinship matrix). It fixes a problem in the old code that occured when heritability estimates where exactly 0.
File Dev_test_new_functions.R has some code for testing the updated functions. 
