sim.dir <- ""
require(GWASTools)
n <- 3000 ## or whatever n you want!

## from setting 2 in the supplementary material
grm.family <- matrix(c(1, 0.05, 0.05, 0.05, 1, 0.1, 0.05, 0.1, 1), nrow = 3)
## from setting 1 in the supplementary material
grm.family <- matrix(c(1, 0.4, 0.5, 0.4, 1, 0.6, 0.5, 0.6, 1), nrow = 3)

house.family <- matrix(c(1,1,1,1,1,1, 1,1,1), nrow = 3)
grm <- matrix(0, nrow = n, ncol = n)
hh.mat <- matrix(0, nrow = n, ncol = n)
n.family <- n/3
n.in.family <- 3
for (i in 1:n.family){
	inds.family <- ((i-1)*n.in.family + 1):(i*n.in.family)
	grm[inds.family,inds.family] <- grm.family
	hh.mat[inds.family,inds.family] <- house.family
}

colnames(grm) <- rownames(grm) <- colnames(hh.mat) <- rownames(hh.mat) <- paste0("p", 1:n)

### another way to simulate a positive definite matrix that resembles a grm:
# grm <- matrix(rnorm(n*n, sd = 0.1), n,n )
# grm <- crossprod(grm)
# diag(grm) <- 1
# in the paper I used a kinship matrix estimated using the HCHS/SOL data in some simulations. 
# this matrix is not provided on github. 



covMatList <- list(grm = grm, hh = hh.mat)


varComp <- c(40, 15, 100)
covMat <- covMatList[[1]]*varComp[1] + covMatList[[2]]*varComp[2] + diag(n)*varComp[3]

### if want to simulate the covariance with only a kinship matrix, use this code: 
### The packages heritability and ALBI can only work with a single matrix. 
# covMat <- covMatList[[1]]*varComp[1] + diag(nrow(covMatList[[1]]))*varComp[3]



svd.covMat <- svd(covMat)
sqrt.covMat <- svd.covMat$v %*% diag(sqrt(svd.covMat$d)) %*% t(svd.covMat$v)



dat <- data.frame(scanID = rownames(grm), EV1 = rnorm(n), stringsAsFactors = FALSE)
scanAnnot <- ScanAnnotationDataFrame(dat)




save(scanAnnot, file = paste0(sim.dir, "/simulated_data/scanAnnot_for_sim.RData"))
save(covMatList, file = paste0(sim.dir, "/simulated_data/covMatList.RData" ))
save(sqrt.covMat, file = paste0(sim.dir, "/simulated_data/sqrt_covMat.RData" ))

