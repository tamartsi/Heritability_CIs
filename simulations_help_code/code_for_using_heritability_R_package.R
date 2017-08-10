# This is the command I used to submit the job to the cluster. Should be different for other people who use different clusters. 
#  qsub -q olga-long.q -pe local 4 -p -2 -j y -o trash -N herpL_k -t 1:1000 /projects/users/tsofer/runRscript_array.sh 20170228_large_samp.R
args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1]) 

sim.dir <- ""
require(GWASTools)
require(CompQuadForm)
require(heritability)

output.folder <- paste0(sim.dir, "/output/")
covMatList <- getobj(paste0(sim.dir, "/simulated_data/covMatList.RData"))
sqrt.covMat <- getobj((paste0(sim.dir, "simulated_data/sqrt_covMat.RData"))
scanAnnot <- getobj((paste0(sim.dir, "simulated_data/scanAnnot_for_sim.RData"))


scans.scanAnnot <- pData(scanAnnot)$scanID

beta <- c(2, 3)


	set.seed(i)
	ind.err <- rnorm(nrow(covMatList[[1]]))
	err <- sqrt.covMat %*% ind.err
	rownames(err) <- rownames(covMatList[[1]])
	
	X <- cbind(rep(1, length(err)), scanAnnot$EV1)
	y <- err + X %*% beta
	rownames(y) <- rownames(covMatList[[1]])	

	res <- marker_h2(data.vector = y, geno.vecto = rownames(y), covariates = X[,2], K = covMatList[[1]])

	

save(res , file = paste0(output.folder,  "her_package_res_num_", i, ".RData"))

### for this package I submitted one job at a time and saved results in RData files. 

## to collect results I used this script:

require(GWASTools)

outfiles <- list.files(output.folder, full.names = TRUE)

res <- data.frame(type = "REML", software = "heritability", VC_type = "kinship", transform = rep(c("identity", "log"), each = length(outfiles)), Est = NA, low = NA, high = NA)

first.ind.identity <- 1
first.ind.log <- 1+ length(outfiles)


for (i in 1:length(outfiles)){
	temp <- getobj(outfiles[i])
	res[first.ind.identity + i - 1,"Est"] <- temp$h2
	res[first.ind.log + i - 1,"Est"] <- temp$h2


	res[first.ind.identity + i - 1,"low"] <- temp$conf.int1[1]
	res[first.ind.log + i - 1,"low"] <- temp$conf.int2[1]


	res[first.ind.identity + i - 1,"high"] <- temp$conf.int1[2] 
	res[first.ind.log + i - 1,"high"] <- temp$conf.int2[2] 



}

write.table(res, file = "")
