size <- function(low, high){
	inds <- which(low == high)
	if (length(inds) > 0){
		high <- high[-inds]
		low <- low[-inds]
	}
	x <- mean(high - low )		

	formatC(x, digits = 2, format = "f")
}

RMSE <- function(est, true.val){
	
	formatC(sqrt(mean((est - true.val)^2)), digits =2, format = "f")
	}
	

coverage <- function(low, high, true.val){
	formatC(mean(low <= true.val & true.val <= high ), digits =2, format = "f")
}


n <- 

true.hh.vc <- 15
true.kin.vc <- 40
true.resid.vc <- 100
all.var <- true.hh.vc + true.kin.vc + true.resid.vc

sim.dir <- ""
dir <-  paste0(sim.dir, "/output/")

res.he <- read.table(file.path(dir, "estimated_VCs.txt"))
colnames(res.he) <- c("V1", "error", "kinship", "hh", "est.her", "her.pval", "kin.low", "kin.high", "her.low", "her.high", "est.hh.prop", "hh.pval", "hh.low", "hh.high", "hh.prop.low" ,"hh.prop.high")

res.genesis <- read.table(file.path(dir, "estimated_VCs_reml.txt"))
colnames(res.genesis) <- c("V1",  "kinship", "hh", "error", "est.her", "her.low", "her.high", "est.hh.prop", "hh.prop.low" ,"hh.prop.high")
inds <- which(res.genesis$kinship == 0)
if (length(inds) > 0){
	res.genesis$her.low[inds] <- 0
	res.genesis$her.high[inds] <- 0
}

inds <- which(res.genesis$hh == 0)
if (length(inds) > 0){
	res.genesis$hh.prop.low[inds] <- 0
	res.genesis$hh.prop.high[inds] <- 0
}

res.genesis$her.low <- pmax(res.genesis$her.low, 0)
res.genesis$her.high <- pmin(res.genesis$her.high, 1)
res.genesis$hh.prop.low <- pmax(res.genesis$hh.prop.low, 0)
res.genesis$hh.prop.high <- pmin(res.genesis$hh.prop.high, 1)


tab <- data.frame(n = NA, coverage = NA, width = NA, RMSE = NA, stringsAsFactors = FALSE)

tab["HE-kin",] <- c(n, coverage(res.he$her.low, res.he$her.high, true.kin.vc/all.var), size(res.he$her.low, res.he$her.high), RMSE(res.he$est.her , true.kin.vc/all.var))

tab["HE-hh",] <- c(n, coverage(res.he$hh.prop.low, res.he$hh.prop.high, true.hh.vc/all.var), size(res.he$hh.prop.low, res.he$hh.prop.high), RMSE(res.he$est.hh.prop , true.hh.vc/all.var))


tab["GENESIS-kin",] <- c(n, coverage(res.genesis$her.low, res.genesis$her.high, true.kin.vc/all.var), size(res.genesis$her.low, res.he$her.high), RMSE(res.genesis$est.her , true.kin.vc/all.var))

tab["GENESIS-hh",] <- c(n, coverage(res.genesis$hh.prop.low, res.genesis$hh.prop.high, true.hh.vc/all.var), size(res.genesis$hh.prop.low, res.genesis$hh.prop.high), RMSE(res.genesis$est.hh.prop , true.hh.vc/all.var))


tab <- tab[-1,]

tab.all <- tab



### code to restructure and put in latex table: 
tab.all.kin <- tab.all[grep("kin", rownames(tab.all)),]
tab.all.hh <- tab.all[grep("hh", rownames(tab.all)),]
colnames(tab.all.kin) <- c("n", paste("kin",colnames(tab.all.kin)[2:4], sep = "-"))
colnames(tab.all.hh) <- c("n", paste("hh",colnames(tab.all.hh)[2:4], sep = "-"))

tab.all.re <- cbind(tab.all.kin, tab.all.hh[,2:4])
tab.all.re$method <- rep(c("HE", "GENESIS"), 4) ## depends on the actual number of simulations used
tab.all.re <- tab.all.re[,c("n", "method", colnames(tab.all.re[,c(2:7)]))]

write.csv(tab.all.re, file = paste0(dir, "/all_results_small_vals_wide.csv"))
require(xtable)
xres <- xtable(tab.all.re)
print.xtable(xres, type = "latex", file = paste0(dir, "/all_results_small_vals_wide.tex"))

