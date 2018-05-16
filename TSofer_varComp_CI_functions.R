require(CompQuadForm)



### calculate variance commponents from Y, W, and covMatList

calc.vc <- function(Y, W, covMatList, verbose = TRUE, max.iter = 20, eps = 1e-4){
	
	### assume no missing data for now
	
	
	
	n <- length(Y)
	
	if (verbose) message("Transforming covariance matrices in preparation for analysis...")
	
	n.cov.mat <- length(covMatList)
	if (is.null(names(covMatList))){
		message(paste("covMatList does not have names, naming matrices by ", paste0("A", 1:n.cov.mat)))
		names(covMatList) <-  paste0("A", 1:n.cov.mat)
	}
	names.cov.mat <- names(covMatList)
	
	covMatVecList <- vector(mode = "list", length = n.cov.mat)
	names(covMatVecList) <- names.cov.mat
	for (i in 1:n.cov.mat){
		covMatVecList[[i]] <- covMatList[[i]][upper.tri(covMatList[[i]], diag = FALSE)]
	}

			

	XtX <- matrix(n, nrow = 1 + n.cov.mat, ncol = 1 + n.cov.mat)
	for (i in 1:n.cov.mat){
	        for (j in 1:n.cov.mat){
	                XtX[i + 1, j + 1] <- sum(covMatVecList[[i]]*covMatVecList[[j]])  + n
	
	        }
	}
	XtXinv <- solve(XtX)


	if (verbose) message("Fitting the model...")
	
	########  fitting the model:
	converged <- FALSE
	n.iter <- 0

	mod.lm <- lm(Y ~ -1 + W)
	beta <- coef(mod.lm)
	fits <- tcrossprod(W, t(beta))
	residM <- as.vector(Y - fits)
	sum.residM.sq <- sum(residM^2)
	resids.mat <- residM %*% t(residM)
	epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = FALSE)]

	Xerr <- matrix(sum.residM.sq, nrow = 1 + n.cov.mat, ncol = 1)
	for (i in 1:n.cov.mat){
		Xerr[i + 1,1] <- Xerr[i + 1,1] + sum(covMatVecList[[i]]*epsilon.vec )
	}
	vc.est.cur <- XtXinv  %*% Xerr
	vc.est.cur <- pmax(vc.est.cur, 0)
	rownames(vc.est.cur) <- c("error", names.cov.mat )
	
	if (verbose) {
		message("Initial variance component estimates ")
		print(vc.est.cur)	
	}

    while(!converged){
		V <-   diag(rep(vc.est.cur[1], n)) 
		for (i in 1:n.cov.mat){
			V <- V + vc.est.cur[i + 1]*covMatList[[i]]
		}
		# cholesky decomposition
		cholV <- chol(V)
		# inverse
		Vinv <- chol2inv(cholV)
		VinvW <- crossprod(Vinv,W)
		cholWtVinvW <- chol(crossprod(W, VinvW))
		WtVinvWInv <- chol2inv(cholWtVinvW)
		beta <- crossprod(WtVinvWInv, crossprod(VinvW,Y))
		fits <- tcrossprod(W, t(beta))
		residM <- as.vector(Y - fits)
		sum.residM.sq <- sum(residM^2)
		resids.mat <- residM %*% t(residM)
		epsilon.vec <- resids.mat[upper.tri(resids.mat, diag = FALSE)]

		Xerr <- matrix(sum.residM.sq, nrow = 1 + n.cov.mat, ncol = 1)
		for (i in 1:n.cov.mat){
			Xerr[i + 1,1] <- Xerr[i + 1,1] + sum(covMatVecList[[i]]*epsilon.vec )
		}

		vc.est.new <- XtXinv  %*% Xerr
		vc.est.new <- pmax(vc.est.new, 0)
		rownames(vc.est.new) <- c("error", names.cov.mat )
		n.iter <- n.iter + 1
        if (verbose) cat("iteration count ", n.iter, "current estimates ", vc.est.new, "\n")

		if (max(abs(vc.est.new - vc.est.cur) < eps) | n.iter == max.iter) converged <- TRUE else{
			vc.est.cur <- vc.est.new
			}
        }
        
        if (verbose) cat("Finished estimation, preparing return values...", "\n")
        ### convenient return values:
        
        #### now to test variance component: Q = vc.est.cur[i],  we need eigenvalues of matrix K, that has the matrix from 
        ### the quadratic form, times the variance of the residuals. 
        ## let My = (y -fits). Then var(My) = M %*% V %*% t(M), which comes to:
        V.My <- V  -  W %*% tcrossprod(WtVinvWInv, W)
	
		vc.est.cur <- drop(vc.est.cur)
		prop.var <- vc.est.cur/sum(vc.est.cur)
		names(prop.var) <- names(vc.est.cur)
        
        return(list(VC = vc.est.cur, propVar = prop.var, var.resids = V.My, XtXinv = XtXinv))
}


#### test variance components, and calculate confidence intervals for variance components and ratio of variance component with total variance. 
#### can calculate test and calculate confidence intervals for sums of variance components. In this case test.mat.name should have more than one string. 
#### (i.e. names of all variance components to test).
testVCcalcCIs <- function( covMatList, test.mat.name , VC.est, XtXinv, var.resids, verbose = TRUE, CI.prob = 0.95, prob.eps = 0.01,  end.point.eps = 1e-3){
	stopifnot(all(is.element(test.mat.name, names(covMatList))))
	
	if (verbose) message("Preparing matrices for analysis...")
	
	total.var <- sum(VC.est)
	n <- nrow(var.resids)
	
	n.cov.mat <- length(covMatList)
	names.cov.mat <- names(covMatList)
	covMatNoDiag <- vector(mode = "list", length = n.cov.mat)
	names(covMatNoDiag) <- names.cov.mat
	for (i in 1:length(covMatNoDiag)){
		covMatNoDiag[[i]] <- covMatList[[i]]
		diag(covMatNoDiag[[i]]) <- 0
	}

	
	ind.mats <- match(test.mat.name, names.cov.mat)
	ind.in.XtX <- ind.mats + 1
	
	if (verbose) message("Computing the appropriate quadratic form of the variance component (or sums of) estimator...")
	
	quad.mat <- matrix(0, ncol = ncol(covMatList[[1]]), nrow = nrow(covMatList[[1]]))
	for (i in 1:n.cov.mat){
		quad.mat <- quad.mat + sum(XtXinv[ind.in.XtX,i + 1])/2*covMatNoDiag[[i]]	
	}
	
	mat <- quad.mat %*% var.resids
	lambda <- eigen(mat, symmetric=TRUE, only.values = TRUE)$values
	VC.pval <- 1-davies(0, lambda, acc=1e-6)$Qq

	
	if (verbose) message("Performing binary search for endpoints of the confidence interval of the variance component (or sum of VC)...")
	
	low.CI.prob <- 1-(1-CI.prob)/2
	high.CI.prob <- (1-CI.prob)/2
	
	########  binary search for low.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- total.var; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
	while (abs(surv.mid - low.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > low.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
		}
		if (val.high == val.low) break
	}

	low.CI.val <- val.mid
	
	
	########  binary search for high.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- total.var; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
	while (abs(surv.mid - low.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > high.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- davies(val.mid, lambda, acc=1e-6)$Qq
		}
		if (val.high == val.low) break
	}

	high.CI.val <- val.mid

	#######  calculate confidence interval for the ratio:
	if (verbose) message("Performing binary search for endpoints of the confidence interval of the ratio between the variance component (or sum of VC) and the total variance...")	
	var.mat <- var.resids/n

	########  binary search for low.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- 1; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
	while (abs(surv.mid - low.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > low.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
		}
		if (val.high == val.low) break
	}

	low.ratio.CI.val <- val.mid
	
	
	########  binary search for high.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- 1; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
	while (abs(surv.mid - low.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > high.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- surv.ratio(const = val.mid, mat, var.mat)$surv
		}
		if (val.high == val.low) break
	}

	high.ratio.CI.val <- val.mid
	
	if (verbose) message("Preparing return values...")	

	prop.var <- VC.est/sum(VC.est)
	names(prop.var) <- names(VC.est)
	VC_CI <- c(low.CI.val, high.CI.val)
	VC_ratio_CI <- c(low.ratio.CI.val, high.ratio.CI.val)
	
	
	return(list(VC = sum(VC.est[ind.in.XtX ]), propVar = sum(VC.est[ind.in.XtX ])/sum(VC.est), VCpval = VC.pval, VC_CI = VC_CI, VC_ratio_CI = VC_ratio_CI, confidence = CI.prob))
	
	
	
}



returnEigen <- function(kin.est.mat, var.mat, const){
	D <- kin.est.mat - const*var.mat
	
	lambda <- eigen(D, symmetric = TRUE, only.values = TRUE)$values	
}



findSaddlePoint <- function(lambda, c1 = NULL, c2 = NULL, eps = 1e-6){
	all.positive <- all(lambda > 0)
	all.negative <- all(lambda < 0)
	if (all.negative) return(list(omega = NA, message = "all eigenvalues are negative, no finite solution for saddelpoint defining equation"))
	if (all.positive) return(list(omega = NA, message = "all eigenvalues are positive, no finite solution for saddelpoint defining equation"))
	
	
	c1 <- 1/(2*min(lambda))
	c2 <- 1/(2*max(lambda))	
	
	f.val <- function(lambda, cur.omega){
		return(sum(lambda/(1-2*lambda*cur.omega)))
	}
	
	
	cur.low <- c1
	cur.high <- c2
	
	f.cur.low <- f.val(lambda, cur.low)
	f.cur.high <- f.val(lambda, cur.high)
	
	if (f.cur.low > 0 | f.cur.high < 0) stop("Endpoints have the same function sign")
	
	converged <- FALSE
	
	while(!converged){
		mid <- (cur.low +  cur.high)/2
		f.mid <- f.val(lambda, mid)
		if (f.mid < 0){
			cur.low <- mid
			f.cur.low <- f.mid
		} else{
			cur.high <- mid
			f.cur.high <- f.mid
		}
		
		if (abs(f.mid) < eps) converged <- TRUE

	}
	
	return(mid)
	
}




calc.z.hat <- function(lambda, saddlepoint){
	vec1 <- lambda^2
	vec2 <- 2*saddlepoint*lambda
	vec3 <- (1-vec2)^2
	vec4 <- vec1/vec3
	z.hat <- saddlepoint*sqrt(2*sum(vec4)) 
	return(z.hat)
	
}

calc.zeta.hat <- function(lambda, saddlepoint){
	vec1 <- 2*saddlepoint*lambda
	vec2 <- log(1-vec1)
	zeta.hat <- sign(saddlepoint)*sqrt(sum(vec2))
	return(zeta.hat)
}



##### calculates the survival probability of const under the distribution of the ratio between the quadratic 
##### forms used to estimate the kinship variance and the total variance. 
surv.ratio <- function(const, kin.est.mat, var.mat){
	lambda <- returnEigen(kin.est.mat, var.mat, const)
	omega <- findSaddlePoint(lambda)
	if (is.na(omega)){
		if (lambda[1] < 0) reutrn(1) else return(0)
	}
	z.hat <- calc.z.hat(lambda, omega)
	zeta.hat <- calc.zeta.hat(lambda, omega)
	surv <- 1 - pnorm(zeta.hat) + dnorm(zeta.hat)*(1/z.hat - 1/zeta.hat) 
	return(list(surv= surv, saddlepoint = omega, z.hat = z.hat, zeta.hat = zeta.hat))
}




##### here we prepare data for estimating confidence intervals for heritability quickly
#### kinshipMat is the matrix of correlations that we want to test its variance component
#### vc is a vector with estimated error variance, and estimated kinship variance
### kinship.vc.name is the name of the estimated variance component in the vector 
prepareDataForHeritabilityEst <- function(kinshipMat, vc, kinship.vc.name){
	if (length(vc) != 2) stop("Two variance components required, error variance and variance corresponding to the correlation matrix argument")
	if (!is.element(kinship.vc.name, names(vc))) stop("Provided the name kinship variance component does not match any name in the vc vector")
	
	kinshipMat.no.diag <- kinshipMat
	diag(kinshipMat.no.diag) <- 0
	lambda.0 <- eigen(kinshipMat.no.diag, symmetric = TRUE, only.values = TRUE)$values
	kappa <- sum(kinshipMat.no.diag^2)
	ind.kinshipVar.vc <- which(names(vc) == kinship.vc.name)
	kinship.var <- vc[ind.kinshipVar.vc]
	error.var <- vc[-ind.kinshipVar.vc]
	
	return(list(lambda.0 = lambda.0, kappa = kappa, kinship.var = kinship.var, error.var = error.var))
}


#### A quick implementation of the heritability confidence intervals algorithms
#### to be used when kinship matrix is the only correlation matrix between individuals
#### (other matrices can be used as well)
#### This function recived a list of "prep.dat" that is either the output of the function prepareDataForHeritabilityEst
#### or the function prepareDataForHeritabilityEst.META. The later function pepare data from multiple studies. 
heritability.CI <- function(prep.dat, CI.prob= 0.95, prob.eps = 0.01, end.point.eps = 1e-3){
	if (CI.prob > 1) stop("Confidence intervals probability should be smaller than 1")
	if (CI.prob <= 0) stop("Confidence intervals probability cannot be zero or less")
	
	
	lambda.0 <- prep.dat$lambda.0
	kappa <- prep.dat$kappa
	kinship.var <- prep.dat$kinship.var 
	error.var <- prep.dat$error.var
	
	n <- length(lambda.0)
	total.var <- kinship.var + error.var
	
	calc.surv.for.const <- function(const){
		cur.lambda <- (lambda.0*(lambda.0*kinship.var + error.var + kinship.var))/kappa - const*(lambda.0*kinship.var + kinship.var + error.var)/n
		if (all(cur.lambda < 0)) surv <- 1
		if (all(cur.lambda > 0)) surv <- 0
		saddlepoint <- findSaddlePoint(cur.lambda)
		z.hat <- calc.z.hat(cur.lambda, saddlepoint )
		zeta.hat <- calc.zeta.hat(cur.lambda, saddlepoint )
		surv <- 1 - pnorm(zeta.hat) + dnorm(zeta.hat)*(1/z.hat - 1/zeta.hat) 
		return(surv)
	}
		
	
	
	
	low.CI.prob <- 1-(1-CI.prob)/2
	high.CI.prob <- (1-CI.prob)/2
	
	########  binary search for low.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- 1; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid <- calc.surv.for.const(val.mid)
	while (abs(surv.mid - low.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > low.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- calc.surv.for.const(val.mid)
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- calc.surv.for.const(val.mid)
		}
		if (val.high == val.low) break
	}

	low.CI.val <- val.mid
	
	
	########  binary search for high.CI point
	val.low <- 0 ; surv.low <- 1
	val.high <- 1; surv.high <- 0
	val.mid <- (val.high + val.low)/2 ; surv.mid  <- calc.surv.for.const(val.mid)
	while (abs(surv.mid - high.CI.prob) > prob.eps & abs(val.mid - val.low) > end.point.eps & abs(val.mid - val.high) > end.point.eps) {
		if (surv.mid > high.CI.prob){
			val.low <- val.mid; surv.low <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- calc.surv.for.const(val.mid)
		} else{ # surv.mid < low.CI.prob
			val.high <- val.mid; surv.high <- surv.mid
			val.mid <- (val.high + val.low)/2 ; surv.mid <- pmax(0, calc.surv.for.const(val.mid))
		}
		if (val.high == val.low) break
	}

	high.CI.val <- val.mid

	CI <- c(low.CI.val, high.CI.val)
	
	return(list(CI = CI, CI.prob = CI.prob))

}


##### reads in results from study-specific files and prepare data for calculating heritability and its confidence interval.
prepareDataForHeritabilityEst.META <- function(file.dir){
 	files <- list.files(file.dir, full.names = TRUE)
	lambda.0.list <- vector(mode = "list", length = length(files))
	kappa.vec <- kinship.var.vec <- error.var.vec <- n.vec <- rep(NA, length(files))
	
	for (i in 1:length(files)){
		for (j in 1:5){ ### assume the correct number of rows in the file!
			temp <- read.table(files[i], nrow  = 1, skip = j-1)
			if (temp[1] == "lambda.0") lambda.0.list[[i]] <- as.numeric(temp[-1])
			if (temp[1] == "kappa") 	kappa.vec[i] <- as.numeric(temp[-1])
			if (temp[1] == "kinship.var") kinship.var.vec[i] <- as.numeric(temp[-1])
			if (temp[1] == "error.var") error.var.vec[i] <- as.numeric(temp[-1])
			if (temp[1] == "n") n.vec[i] <- as.numeric(temp[-1])
		}
	}
	
	lambda.0 <- do.call( c, lambda.0.list)
	kappa <- sum(kappa.vec)
	kinship.var <- sum(kappa.vec*kinship.var.vec)/sum(kappa.vec)
	error.var <- sum(n.vec*error.var.vec)/sum(n.vec)
	
	return(list(lambda.0 = lambda.0, kappa = kappa, kinship.var = kinship.var, error.var = error.var))
	
}






