


n <- 3000 ## or whatever n you want!

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
covMatList <- list(grm = grm, hh = hh.mat)
varComp <- c(40, 15, 100)
covMat <- covMatList[[1]]*varComp[1] + covMatList[[2]]*varComp[2] + diag(n)*varComp[3]
svd.covMat <- svd(covMat)
sqrt.covMat <- svd.covMat$v %*% diag(sqrt(svd.covMat$d)) %*% t(svd.covMat$v)
dat <- data.frame(scanID = rownames(grm), EV1 = rnorm(n), stringsAsFactors = FALSE)

beta <- c(2, 3)
ind.err <- rnorm(nrow(covMatList[[1]]))
err <- sqrt.covMat %*% ind.err
rownames(err) <- rownames(covMatList[[1]])

X <- cbind(rep(1, length(err)), dat$EV1)
y <- err + X %*% beta
rownames(y) <- rownames(covMatList[[1]])


est_vcs <- calc.vc(Y=y, W=X, covMatList = covMatList)

est_CIs <- calcCIsForVCs(covMatList, 
                         eval.vc.name = "grm", 
                         VC.est = est_vcs$VC, 
                         XtXinv = est_vcs$XtXinv, 
                         var.resids = est_vcs$var.resids)
est_CIs



est_CIs <- calcCIsForVCs(covMatList, 
                         eval.vc.name = "hh", 
                         VC.est = est_vcs$VC, 
                         XtXinv = est_vcs$XtXinv, 
                         var.resids = est_vcs$var.resids)
est_CIs
