# qsub -q olga-long.q -p -2 -j y -o trash -N small_her -t 1:10 /projects/users/tsofer/runRscript_array.sh 20170710_small_samp.R
args <- commandArgs(trailingOnly=TRUE)
start.seed <- (as.numeric(args[1])-1)*100 + 1
n.sim <- 100

sim.dir <- ""
require(GWASTools)
require(CompQuadForm)

source('20161010_varComp_CI_functions_with_eigen.R', chdir = TRUE)
output.folder <- paste0(sim.dir, "/output/")
covMatList <- getobj(paste0(sim.dir, "/simulated_data/covMatList.RData"))
sqrt.covMat <- getobj((paste0(sim.dir, "simulated_data/sqrt_covMat.RData"))
scanAnnot <- getobj((paste0(sim.dir, "simulated_data/scanAnnot_for_sim.RData"))


scans.scanAnnot <- pData(scanAnnot)$scanID

beta <- c(2, 3)

for (i in start.seed:(start.seed+99)){
	set.seed(i)
	ind.err <- rnorm(nrow(covMatList[[1]]))
	err <- sqrt.covMat %*% ind.err
	rownames(err) <- rownames(covMatList[[1]])
	
	X <- cbind(rep(1, length(err)), scanAnnot$EV1)
	y <- err + X %*% beta
	rownames(y) <- rownames(covMatList[[1]])
	n <- length(y)

	pData(scanAnnot)$y <- drop(y)

  res <- estVarComp(scanAnnot=scanAnnot,
             covMatList=covMatList,
             outcome="y",
            covar.vec="EV1")

	CIs <- estVarCompCI(res)




	
	row <- as.numeric(c(res$varComp, CIs["V_grm",], CIs["V_hh",]))
	write.table(t(row), file = paste0(output.folder, "estimated_VCs_reml.txt"), append = T, col.names = F)
	print(i)

}




