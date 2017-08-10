# The command that I used to submit to my cluster. Use your own!
# qsub -q olga-long.q -p -2 -j y -o trash -N small_her -t 1:10 /projects/users/tsofer/runRscript_array.sh 20170710_small_samp.R
args <- commandArgs(trailingOnly=TRUE)
start.seed <- (as.numeric(args[1])-1)*100 + 1
n.sim <- 100

sim.dir <- ""  ## put your own!
require(GWASTools)
require(CompQuadForm)

source('Heritability_CIs/Tsovar_varComp_CI_functions.R')
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

	VC.est <- calc.vc(Y = y, W = X, covMatList = covMatList)
	
	vc.inference.k <- testVCcalcCIs(covMatList, "grm", VC.est$VC, VC.est$XtXinv, VC.est$var.resids)
	vc.inference.h <- testVCcalcCIs(covMatList, "hh", VC.est$VC, VC.est$XtXinv, VC.est$var.resids)


	
	row <- c(VC.est$VC, vc.inference.k$propVar, vc.inference.k$VCpval, vc.inference.k$VC_CI, vc.inference.k$VC_ratio_CI, vc.inference.h$propVar, vc.inference.h$VCpval, vc.inference.h$VC_CI, vc.inference.h$VC_ratio_CI)
	
	print(i)

	if (i == start.seed) rows <- row else rows <- rbind(rows, row)
}


write.table(rows, file = paste0(output.folder, "estimated_VCs.txt"), append = T, col.names = F)


#### If we used only a single correlation matrix (e.g. if length(covMatList) = 1), 
####  we would estimate confidence intervals using these functions: 
prep.dat <- prepareDataForHeritabilityEst(covMatList[[1]], VC.est$VC, names(covMatList)[1])
her.CI <- heritability.CI(prep.dat)
## this is instead of using testVCcalcCIs()!

