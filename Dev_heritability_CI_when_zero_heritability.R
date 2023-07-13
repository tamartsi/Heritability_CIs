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