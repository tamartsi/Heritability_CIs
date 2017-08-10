
### For a set of simualtions in which we used a fixed GRM/kinship matrix, we first calculate the
### eigenvalues and eigenvectors of the matrix (once) and save them in an input folder.

sim.dir <- ""
require(GWASTools)


output.folder <- paste0(sim.dir, "/ALBI/output_files/")
input.folder <- paste0(sim.dir, "/ALBI/input_files/")

covMatList <- getobj(paste0(sim.dir, "/simulated_data/covMatList.RData"))
sqrt.covMat <- getobj((paste0(sim.dir, "simulated_data/sqrt_covMat.RData"))
scanAnnot <- getobj((paste0(sim.dir, "simulated_data/scanAnnot_for_sim.RData"))

eigen.kinship <- eigen(covMatList[[1]], symmetric = TRUE)
lambda <- (eigen.kinship$values)
write.table(t(lambda), file = paste0(input.folder, "sim_eigenvalues.txt"), quote = F, row.names = F, col.names = F)
eigenvec <- eigen.kinship$vectors
write.table(eigenvec, file = paste0(input.folder, "sim_eigenctors.txt"), quote = F, row.names = F, col.names = F)


### now the ALBI calcualte confidence intervals independently of the actual simulated outcome.
## this is the command that calcualtes confidence intervals for 100 possible values of estimated heritability:

albi.command <- paste0("python2.7 ", paste0(sim.dir, "/ALBI/albi-master-2/albi.py"),  "--kinship_eigenvalues ", paste0(sim.dir, "/ALBI/input_files/sim_eigenvalues.txt "),  "--kinship_eigenvectors ", paste0(sim.dir, "/ALBI/input_files/sim_eigenctors.txt "),   "--covariates ", paste0(input.folder, "covars_",seed,  ".txt "),   " --estimate_grid ", 100, " --confidence 0.95 --output_filename ", paste0(output.folder, "CIs.txt"))

# We can submit this command through R:

system(albi.command)

### read the estimated CIs:
CIs <- read.table(paste0(output.folder,  "CIs.txt", header = TRUE)

### to obtain confidence intervals from heritability that was estimated using REML (e.g. using the GENESIS package), suppse 
### that the heritability value is her. 
## Then: 

## prepare a data frame for results from heritability estimates in outfiles (or another set up that you made). 
res <- data.frame(type = rep("Bootstrap", length = n.sim), software = "ALBI", VC_type = "kinship", Est = NA, low = NA, high = NA)


for (i in 1:n.sim){
	her <- #### load the estimated heritability for this simulation

	res[i,"Est"] <- her
	
	CI.line <- which.min(abs(CIs$Estimate - res[i,"Est"]))
	
	res[i,"low"]  <- CIs[CI.line,"CI_lower_bound"]
	res[i,"high"] <- CIs[CI.line,"CI_upper_bound"] 

}
