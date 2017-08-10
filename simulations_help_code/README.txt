The files in this folder provide code and instructions to run simulations described in the manuscript “Confidence intervals for heritability via Haseman-Elston regression”. 

Notes:
1. The simulations in the main manuscript use information from the Hispanic Community Health Study/Study of Latinos (estimated kinship matrix, household matrix, community block unit). There are NOT provided here. 
2. Instructions provided here are useful to replicate simulations that are provided in the supplementary material of the manuscript. 

File descriptions:
1. prepare_matrices_for_simulations.R: provides code to generate correlation matrices and some additional required data to run the simulations. 
2. run_HE_analysis.R: code for using functions from the Heritability_CIs repository based on the Haseman-Elston (HE) regression to run simulations. 
3. run_REML_GENESIS_analysis.R: code for using functions from the GENESIS R package implementing REML estimates of the variance components to run simulations.
4. calculate_results_measures: calculates measures of simulation results based on output from the functions run_HE_analysis.R and run_REML_GENESIS_analysis.R. 
5. code_for_using_ALBI_package.R: as the file name! code for preparing some input arguments, submission line, and how to prepare results. 
6. code_for_using_heritability_R_package.R: as the file name! code for preparing some input arguments, submission line, and how to prepare results. 

For help contact me at tsofer@uw.edu or at other current email address (it can be always found using google). 

 

 