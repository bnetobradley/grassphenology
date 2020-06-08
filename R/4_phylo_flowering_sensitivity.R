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

## plot phylo
#clim1$binomial <- factor(clim1$binomial,levels = td$phy$tip.label)

# estimate phylogenetic signal ####
clim <- list(clim1, clim2, clim3)
climvar <- c("TempMonthPrior", "Temp2MonthPrior", "PPT2MonthPrior")
psig_results <- vector()
pssig_results <- vector()

# skip for loop to avoid running MCMC (unhash line 41 and 46 to run)
for (i in 3) {
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
# priors set as in: https://github.com/uyedaj/EQG2017/blob/master/R/bayouWorkshopNotebook.Rmd
priorOU <- make.prior(td$phy,dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dk="cdpois", dtheta="dnorm"),param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1),dtheta=list(mean=mean(x), sd=1.5*sd(x))))
startpars <- priorSim(priorOU, td$phy, plot=TRUE)$pars[[1]]
priorOU(startpars)
set.seed(1)
mcmcOU <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodeli", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
#mcmcOU$run(1000000)
chainOU <- mcmcOU$load()
chainOU <- set.burnin(chainOU, 0.3)
save(chainOU, file=paste("data/bayouchaini", climvar[i], ".Rdata", sep=""))
mcmcOU2 <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodelii", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
#mcmcOU2$run(1000000)
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

### only significant estimates

clim1s <- clim1 %>% filter(p.value < 0.05)
clim2s <- clim2 %>% filter(p.value < 0.05)
clim3s <- clim3 %>% filter(p.value < 0.05)

# estimate phylogenetic signal ####
clims <- list(clim1s, clim2s, clim3s)
climvar <- c("SigTempMonthPrior", "SigTemp2MonthPrior", "SigPPT2MonthPrior")
spsig_results <- vector()
spssig_results <- vector()

# skip for loop to avoid running MCMC (unhash line 41 and 46 to run)
for (i in 3) {
  sens <- as.data.frame(clims[i])
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
  spsig_results <- cbind(ps, spsig_results)
  spssig_results <- cbind(pss, spssig_results)
  ## results for phylogenetic signal of 3 clim variables in 'psig_results' (method = K) and 'pssig_results' (method = lambda)
  
  # estimate phylogenetic halflife ####
  priorOU <- make.prior(td$phy,dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dk="cdpois", dtheta="dnorm"),param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),dk=list(lambda=10, kmax=50), dsb=list(bmax=1, prob=1),dtheta=list(mean=mean(x), sd=1.5*sd(x))))
  startpars <- priorSim(priorOU, td$phy, plot=TRUE)$pars[[1]]
  priorOU(startpars)
  set.seed(1)
  mcmcOU <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodeliSig", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
  # mcmcOU$run(1000000)
  chainOU <- mcmcOU$load()
  chainOU <- set.burnin(chainOU, 0.3)
  save(chainOU, file=paste("data/bayouchainiSig", climvar[i], ".Rdata", sep=""))
  mcmcOU2 <- bayou.makeMCMC(td$phy, x, SE=se, prior=priorOU, new.dir=F, outname=paste("OUmodeliiSig", climvar[i], sep=""), plot.freq=NULL) # Set up the MCMC
  # mcmcOU2$run(1000000)
  chainOU2 <- mcmcOU2$load()
  chainOU2 <- set.burnin(chainOU2, 0.3)
  save(chainOU2, file=paste("data/bayouchainiiSig", climvar[i], ".Rdata", sep="")) }

# read chains back in
## chain i and ii for temperature month prior
chainisig_tmonth <- loadRData("data/bayouchainiSigSigTempMonthPrior.Rdata")
chainiisig_tmonth <- loadRData("data/bayouchainiiSigSigTempMonthPrior.Rdata")
## chain i and ii for temperature 2 months prior
chainisig_t2month <- loadRData("data/bayouchainiSigSigTemp2MonthPrior.Rdata")
chainiisig_t2month <- loadRData("data/bayouchainiiSigSigTemp2MonthPrior.Rdata")
## chain i and ii for ppt 2 months prior
chainisig_ppt <- loadRData("data/bayouchainiSigSigPPT2MonthPrior.Rdata")
chainiisig_ppt <- loadRData("data/bayouchainiiSigSigPPT2MonthPrior.Rdata")

####
library(grid) 
colfunc <- colorRampPalette(c("#856D89","#FFFFFF", "#C0E9E1"))
g <- make_gradient(deg = 180, n = 500, cols = colfunc(10))
gp <- ggplot(clim1, aes(x = estimate, y = binomial)) + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + geom_errorbar(data=clim1, mapping=aes(x=estimate, xmax=(estimate + std.error), xmin=(estimate - std.error)), alpha = 0.25) + geom_point() + geom_vline(xintercept = 0, linetype = 3) + xlim(-10,10) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank())

clim1$binomial <- factor(clim1$binomial,levels = td$phy$tip.label)
ggplot(clim1, aes(y=binomial, x=estimate)) + geom_point(na.rm = F) + coord_polar(theta = "y") + theme(legend.text = element_blank()) + xlim(20,-16)
subvp <- viewport(width = 0.2, height = 0.78, x = 0.49, y = 0.50)
plot.phylo(td$phy,  cex = 0.76, label.offset = 40)
print(gp, vp = subvp)