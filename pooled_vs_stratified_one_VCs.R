####### Created on September 15 2016 by Tamar Sofer  ###################################
##############      tsofer@uw.edu  or tsofer.uw@gmail.com  #############################
source("TSofer_varComp_CI_functions.R")

##### Construct a positive semidefinite matrice, as a sources of correlation, and normalized it.
##### althought this is not  a realistic matrices to use in genetic association studies, 
##### mathematically we achieve the same goal when using them, so they are useful for demonstrating the algorithms. 


#### Demonstrate both pooled- and meta-analysis
n <- 2000
A <- matrix(rnorm(n*n), nrow = n, ncol = n)
A <- crossprod(A)
A <- diag(1/sqrt(diag(A))) %*% A %*% diag(1/sqrt(diag(A))) 




### 
varComp <- c(100, 40)
names(varComp) <- c("error", "A")

Sigma <- varComp[1]*diag(n)
for (i in 2:length(varComp)){
	Sigma <- Sigma + varComp[i]*get(names(varComp)[i])
}

svd.Sigma <- svd(Sigma)
sqrt.Sigma <- svd.Sigma$v %*% diag(sqrt(svd.Sigma$d)) %*% t(svd.Sigma$v)

covEff <- c(2,3)

ind.err <- rnorm(n)
err <- sqrt.Sigma %*% ind.err
X <- cbind(rep(1, length(err)), rnorm(n))
y <- err + X %*% covEff

#### First, consider the pooled analysis (assume we have all data)"
covMatList <- list(A)
names(covMatList) <- c("A")


########  Estimate variance components, return a few quantities that are helpful for calculating confidence
####### intervals
estVC <- calc.vc(Y = y, W = X, covMatList = covMatList)

estVC$VC
estVC$propVar

#######  Try using the algorithm for multiple correlation matrices:


system.time(A.VC <- testVCcalcCIs(covMatList = covMatList, test.mat.name = "A", VC.est = estVC$VC, XtXinv = estVC$XtXinv, var.resids = estVC$var.resids ))
A.VC

#######  Because we have only one correlation matrix, we can use the fast algorithm. First prepare data:

system.time(prep.dat <- prepareDataForHeritabilityEst(kinshipMat = A, vc = estVC$VC, kinship.vc.name = "A"))
system.time(A.VC.quick <- heritability.CI(prep.dat))
### much faster!
A.VC.quick


########  Demonstration of meta-analysis: divide the data to 4 equal part, prepare each study information and write to file, 
####### to mimic 4 studies. 
####### Meta-analysis function reads from files and prepares data for heritability.CI
output.dir <- file.path(getwd(), "Studies_for_meta")
dir.create(output.dir)
n.studies <- 4
study.n <- n/n.studies
for (i in 1:n.studies){
	study.file <- file.path(output.dir, paste0("study_", i, ".txt"))
	study.inds <- ((i-1)*study.n + 1):(i*study.n)
	
	A.study <- A[study.inds,study.inds]
	y.study <- y[study.inds]
	X.study <- X[study.inds,]
	n.study <- length(study.inds)
	estVC.study <- calc.vc(Y = y.study, W = X.study, covMatList = list(A = A.study))
	prep.dat.study <- prepareDataForHeritabilityEst(kinshipMat = A.study, vc = estVC.study$VC, kinship.vc.name = "A")
	write.table(t(c("lambda.0", prep.dat.study$lambda.0)), file = study.file, row.names = F, col.names = F)
	write.table(t(c("kappa", prep.dat.study$kappa)), file = study.file, row.names = F, col.names = F, append = T)
	write.table(t(c("kinship.var", prep.dat.study$kinship.var)), file = study.file, row.names = F, col.names = F, append = T)
	write.table(t(c("error.var", prep.dat.study$error.var)), file = study.file, row.names = F, col.names = F, append = T)
	write.table(t(c("n", n.study)), file = study.file, row.names = F, col.names = F, append = T)
}

prep.dat <- prepareDataForHeritabilityEst.META(output.dir)
A.VC.meta <- heritability.CI(prep.dat)
A.VC.meta
### wider confidence intervals


