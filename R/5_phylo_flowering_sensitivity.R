# load packages
library(phytools)
library(geiger)
# devtools::install_github("uyedaj/bayou")
require(bayou)

# read in phylo and climate vars ####
phy <- read.tree("Poaceae_hyp1_allC4_maxcred.phy")
clim1 <- read.csv("data/slopes/avetemp_month_prior_slopes.csv")
clim2 <- read.csv("data/slopes/avetemp_2month_prior_slopes.csv")
clim3 <- read.csv("data/slopes/ppt_2month_prior_slopes.csv")

# estimate phylogenetic signal ####
clim <- list(clim1, clim2, clim3)
climvar <- c("TempMonthPrior", "Temp2MonthPrior", "PPT2MonthPrior")
psig_results <- vector()
pssig_results <- vector()

# skip for loop to avoid running MCMC (unhash line 41 and 46 to run)
for (i in 1:3) {
  sens <- as.data.frame(clim[i])
  rownames(sens) <- sens$binomial
  td <- treedata(phy, sens)
  x <- sens$estimate
  names(x) <- sens$binomial
  se <- sens$std.error
  names(se) <- sens$binomial
  ps <- phylosig(td$phy, x = x, method="K", test=TRUE, nsim=1000, se = se, start=NULL,control=list())
  ps$climvar <- paste(climvar[i])
  pss <- phylosig(td$phy, x = x, method="lambda", test=TRUE, se = se, start=NULL,control=list())
  pss$climvar <- paste(climvar[i])
  psig_results <- cbind(ps, psig_results)
  pssig_results <- cbind(pss, pssig_results)
## results for phylogenetic signal of 3 clim variables in 'psig_results' (method = K) and 'pssig_results' (method = lambda)

# estimate phylogenetic halflife ####
priorOU <- make.prior(td$phy,dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dk="cdpois", dtheta="dnorm"),param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1),dtheta=list(mean=mean(x), sd=1.5*sd(x))))
startpars <- priorSim(priorOU, td$phy, plot=TRUE)$pars[[1]]
priorOU(startpars)
set.seed(1)
mcmcOU <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodeli", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
mcmcOU$run(1000000)
chainOU <- mcmcOU$load()
chainOU <- set.burnin(chainOU, 0.3)
save(chainOU, file=paste("data/bayouchaini", climvar[i], ".Rdata", sep=""))
mcmcOU2 <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodelii", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
mcmcOU2$run(1000000)
chainOU2 <- mcmcOU2$load()
chainOU2 <- set.burnin(chainOU2, 0.3)
save(chainOU2, file=paste("data/bayouchainii", climvar[i], ".Rdata", sep="")) }

## read in different chains
## function to load and rename RData file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])   }
## chain i and ii for temperature month prior
chaini_tmonth <- loadRData("data/bayouchainiTempMonthPrior.Rdata")
chainii_tmonth <- loadRData("data/bayouchainiiTempMonthPrior.Rdata")
## chain i and ii for temperature 2 months prior
chaini_t2month <- loadRData("data/bayouchainiTemp2MonthPrior.Rdata")
chainii_t2month <- loadRData("data/bayouchainiiTemp2MonthPrior.Rdata")
## chain i and ii for ppt 2 months prior
chaini_ppt <- loadRData("data/bayouchainiPPT2MonthPrior.Rdata")
chainii_ppt <- loadRData("data/bayouchainiiPPT2MonthPrior.Rdata")

# check chains converged & summarise
summary(chaini_t2month)
summary(chainii_t2month)
RlnL <- gelman.R("lnL", chain1=chaini_t2month, chain2=chainii_t2month, plot=TRUE, type="n", ylim=c(0.9, 10))
Ralpha <- gelman.R("alpha", chain1=chaini_t2month, chain2=chainii_t2month, plot=TRUE, type="n", ylim=c(0.9, 10))
Rsig2 <- gelman.R("sig2", chain1=chaini_t2month, chain2=chainii_t2month, plot=TRUE, type="n", ylim=c(0.9, 10))

# halflife
log(2)/summary(chaini_t2month$alpha)
log(2)/summary(chainii_t2month$alpha)

summary(chaini_tmonth)
summary(chainii_tmonth)
RlnL <- gelman.R("lnL", chain1=chaini_tmonth, chain2=chainii_tmonth, plot=TRUE, type="n", ylim=c(0.9, 10))
Ralpha <- gelman.R("alpha", chain1=chaini_tmonth, chain2=chainii_tmonth, plot=TRUE, type="n", ylim=c(0.9, 10))
Rsig2 <- gelman.R("sig2", chain1=chaini_tmonth, chain2=chainii_tmonth, plot=TRUE, type="n", ylim=c(0.9, 10))

# halflife
log(2)/summary(chaini_tmonth$alpha)
log(2)/summary(chainii_tmonth$alpha)

summary(chaini_ppt)
summary(chainii_ppt)
RlnL <- gelman.R("lnL", chain1=chaini_ppt, chain2=chainii_ppt, plot=TRUE, type="n", ylim=c(0.9, 10))
Ralpha <- gelman.R("alpha", chain1=chaini_ppt, chain2=chainii_ppt, plot=TRUE, type="n", ylim=c(0.9, 10))
Rsig2 <- gelman.R("sig2", chain1=chaini_ppt, chain2=chainii_ppt, plot=TRUE, type="n", ylim=c(0.9, 10))

# halflife
log(2)/summary(chaini_ppt$alpha)
log(2)/summary(chainii_ppt$alpha)

plotSimmap.mcmc(chaini_tmonth, burnin = 0.3, pp.cutoff = 0.3)
plotSimmap.mcmc(chaini_t2month, burnin = 0.3, pp.cutoff = 0.3)
plotSimmap.mcmc(chaini_ppt, burnin = 0.3, pp.cutoff = 0.3)