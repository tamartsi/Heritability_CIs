####### Created on September 15 2016 by Tamar Sofer  ###################################
##############      tsofer@uw.edu    ###################################################
source("TSofer_varComp_CI_functions.R")

##### Construct positive semidefinite matrices, as sources of correlation, and normalized them.
##### althought these are not realistic matrices to use in genetic association studies, 
##### mathematically we achieve the same goal when using them, so they are useful for demonstrating the algorithms. 


#### Demonstrate the use the algorithm when multiple sources of correlation are present:
n <- 2000
A <- matrix(rnorm(n*n), nrow = n, ncol = n)
A <- crossprod(A)
A <- diag(1/sqrt(diag(A))) %*% A %*% diag(1/sqrt(diag(A))) 

B <- matrix(rnorm(n*n), nrow = n, ncol = n)
B <- crossprod(B)
B <- diag(1/sqrt(diag(B))) %*% B %*% diag(1/sqrt(diag(B))) 

C <- matrix(rnorm(n*n), nrow = n, ncol = n)
C <- crossprod(C)
C <- diag(1/sqrt(diag(C))) %*% C %*% diag(1/sqrt(diag(C))) 

###
varComp <- c(100, 40, 15, 10)
names(varComp) <- c("error", "A", "B", "C")

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

covMatList <- list(A,B,C)
names(covMatList) <- c("A", "B", "C")


########  Estimate variance components, return a few quantities that are helpful for calculating confidence
####### intervals
estVC <- calc.vc(Y = y, W = X, covMatList = covMatList)

estVC$VC
estVC$propVar

#######  Calculate confidence intervals for the variance component corresponding to A, 
####### and for the proportion of variance of A.


A.VC <- testVCcalcCIs(covMatList = covMatList, test.mat.name = "A", VC.est = estVC$VC, XtXinv = estVC$XtXinv, var.resids = estVC$var.resids )
A.VC

#######  Calculate confidence intervals for the variance components corresponding to B and C together (sigma^2_b + sigma^2_c), 
####### and for the proportion of variance of B + C.

BplusC.VC <- testVCcalcCIs(covMatList = covMatList, test.mat.name = c("B", "C"), VC.est = estVC$VC, XtXinv = estVC$XtXinv, var.resids = estVC$var.resids )
BplusC.VC







