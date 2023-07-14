## attempting to revise the heritability.CI function to address the setting when 
## heritability is estiamted as 0.
## lower confidence interval is assumed 0, and the funciton only uses binary 
## search for the high point of the confidence interval.


#### A quick implementation of the heritability confidence intervals algorithms
#### to be used when kinship matrix is the only correlation matrix between individuals
#### (other matrices can be used as well)
#### This function receives a list of "prep.dat" that is either the output of the function prepareDataForHeritabilityEst
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
  
  
  if (kinship.var == 0){
    low.CI.prob <- 0
    high.CI.prob <- (1-CI.prob)
    
    low.CI.val <- 0
    
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
    
    
  } else{
    
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
  }
  CI <- c(low.CI.val, high.CI.val)

  
  return(list(CI = CI, CI.prob = CI.prob))
  
}



# Add the same option for the other function
####  calculate confidence intervals for variance components and ratio of variance component with total variance. 
#### can calculate confidence intervals for sums of variance components. In this case eval.vc.name should have more than one string. 
#### (i.e. names of all variance components to estimate confidence intervals for).
calcCIsForVCs <- function( covMatList, eval.vc.name , VC.est, XtXinv, var.resids, verbose = TRUE, CI.prob = 0.95, prob.eps = 0.01,  end.point.eps = 1e-3){
  stopifnot(all(is.element(eval.vc.name, names(covMatList))))
  
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
  
  
  ind.mats <- match(eval.vc.name, names.cov.mat)
  
  est.var.by.eval.vc.name <- sum(VC.est[ind.mats])
  
  ind.in.XtX <- ind.mats + 1
  
  if (verbose) message("Computing the appropriate quadratic form of the variance component (or sums of) estimator...")
  
  quad.mat <- matrix(0, ncol = ncol(covMatList[[1]]), nrow = nrow(covMatList[[1]]))
  for (i in 1:n.cov.mat){
    quad.mat <- quad.mat + sum(XtXinv[ind.in.XtX,i + 1])/2*covMatNoDiag[[i]]	
  }
  
  mat <- quad.mat %*% var.resids
  lambda <- eigen(mat, symmetric=TRUE, only.values = TRUE)$values
  
  # probability of having value lower than the one observed for the variance
  prob.low.var <- 1 - davies(vest.var.by.eval.vc.name, lambda, acc=1e-6)$Qq
  
  
  
  if (verbose) message("Performing binary search for endpoints of the confidence interval of the variance component (or sum of VC)...")
  
  if (est.var.by.eval.vc.name == 0){
    low.CI.prob <- 0
    high.CI.prob <- (1-CI.prob)
    
    # low point of confidence interval is assumed to be zero
    low.CI.val <- 0
  } else{
    low.CI.prob <- 1-min((1-CI.prob)/2, prob.low.var)
    high.CI.prob <- (1-CI.prob) - min((1-CI.prob)/2, prob.low.var)
    
    if (low.CI.prob > 1-(1-CI.prob)/2){
      low.CI.val <- 0
    }  else{
      
      
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
    }
      
  }
  
 
  
  
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
  prob.low.ratio <- 1 - surv.ratio(const = est.var.by.eval.vc.name/total.var, mat, var.mat)$surv
  
  if (est.var.by.eval.vc.name == 0){
    low.CI.prob <- 0
    high.CI.prob <- (1-CI.prob)
    # low point of confidence interval is assumed to be zero
    low.ratio.CI.val <- 0
  } else{
    ########  binary search for low.CI point
    low.CI.prob <- 1-min((1-CI.prob)/2, prob.low.ratio)
    high.CI.prob <- (1-CI.prob) - min((1-CI.prob)/2, prob.low.ratio)
    
    if (low.CI.prob > 1-(1-CI.prob)/2){
      low.ratio.CI.val <- 0
    }  else{
    
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
      }
    }
  
  
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
  
  
  return(list(VC = sum(VC.est[ind.in.XtX ]), propVar = sum(VC.est[ind.in.XtX ])/sum(VC.est),  VC_CI = VC_CI, VC_ratio_CI = VC_ratio_CI, confidence = CI.prob))
  
  
  
}



