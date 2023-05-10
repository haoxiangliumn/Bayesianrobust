## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)
library(Bayesianrobust) # load our package

## -----------------------------------------------------------------------------
# 1 not monotone
a <- c(NA, 1, 4, NA, 5, 9)
b <- c(2, NA, 3, 5, 6, NA)
c <- c(2, 4, 6, NA, 5, 6)
x <- c(1, 2, 3, 4, 5, 6)
yobs <- data.frame(a, b, c)
ismonotone(yobs, x)

## -----------------------------------------------------------------------------
# get data
x <- cbind(1, creatinine_data["In140_AGE"]) # include intercept features
y <- creatinine_data[c("InSC", "InWT", "InCR")] # responses
mon <- ismonotone(y, x)

## -----------------------------------------------------------------------------
mon$monotone

## -----------------------------------------------------------------------------
mon$permute

## -----------------------------------------------------------------------------
yobs <- mon$ynew
xobs <- mon$xnew

## -----------------------------------------------------------------------------
mon$ynew_column

## -----------------------------------------------------------------------------
mon$ynew_row

## -----------------------------------------------------------------------------
iter <- 5000
mcmc <- dageneral(xobs, yobs, m=ncol(yobs), A=diag(0,ncol(yobs)), cw=cw_gamma, iter=iter) 

## -----------------------------------------------------------------------------
library(ggplot2)
lower <- 10 # burn-in
upper <- iter + 1
dim <- 2 # beta(2,1) 
sim <- data.frame(Model=rep(c("gamma"), each=(upper - lower + 1)), beta= mcmc$samples[lower:upper, dim])
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## -----------------------------------------------------------------------------
# constant
mcmc_c <- dageneral(xobs, yobs, m=ncol(yobs),A=diag(0,ncol(yobs)), cw=cw_constant, iter=iter)

## -----------------------------------------------------------------------------
# gamma(a=4, b=4)
cw_gamma4 <- function(di, ri) {
  return(cw_gamma0(d=di, ri=ri, a=4, b=4))
}
# where
#cw_gamma0 <- function(di, ri, a=2, b=2) {
#  return(rgamma(n=1, shape=di / 2 + a, rate=ri / 2 + b))
#}
mcmc_g4 <- dageneral(xobs, yobs, cw=cw_gamma4, iter=iter) 

## -----------------------------------------------------------------------------
xobs_r <- xobs[-27, ]
yobs_r <- yobs[-27, ]

## -----------------------------------------------------------------------------
# constant
mcmc_rc <- dageneral(xobs_r, yobs_r, cw=cw_constant, iter=iter)

## -----------------------------------------------------------------------------
sim <- data.frame(Model=rep(c("gamma2", "p", "gamma4", "rp"), each=(upper - lower + 1)), beta= c(mcmc$samples[lower:upper, dim], mcmc_c$samples[lower:upper, dim], mcmc_g4$samples[lower:upper, dim], mcmc_rc$samples[lower:upper, dim]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## -----------------------------------------------------------------------------
ym <- creatinine_data[c("InSC", "InWT", "InCR")]
ym[29, 1] = NA
ym[30, 1] = NA

monm <- ismonotone(ym, x)

## -----------------------------------------------------------------------------
monm$monotone

## -----------------------------------------------------------------------------
monm$permute

## -----------------------------------------------------------------------------
yobs1 <- monm$ynew
xobs1 <- monm$xnew

## -----------------------------------------------------------------------------
monm$ynew_column

## -----------------------------------------------------------------------------
monm$ynew_row

## -----------------------------------------------------------------------------
mcmcm <- dageneral(xobs1, yobs1, m=ncol(yobs1),A=diag(0,ncol(yobs1)), cw=cw_gamma, iter=iter) 

## ---- eval = FALSE------------------------------------------------------------
#  set.seed(1)
#  library(GeneralizedHyperbolic)
#  library(LaplacesDemon)
#  library(matlib)
#  library(MASS)
#  library(matrixNormal)
#  
#  # generate data
#  library(mvtnorm)
#  n <- 50
#  X <- cbind(1, as.matrix(rnorm(n, 0, 1)))
#  beta <- matrix(c(0.8, 0.4, 0.5, 0.3), nrow=2, ncol=2)
#  var <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
#  error <- matrix(0, nrow=n, ncol=2)
#  w <- rgamma(n, shape=2, rate=2)
#  for (i in 1:n) {
#    error[i, ] <- rmvnorm(1, rep(0, nrow(var)), var / w[i])
#  }
#  yobs <- X %*% beta + error
#  data <- cbind(yobs, X[, 2])
#  data_d <- cbind(X[, 2], yobs)
#  colnames(data_d) <- c("x", "yr1", "yr2")
#  
#  # display data
#  library(GGally)
#  ggpairs(data.frame(data_d))
#  # ggsave("scatter.png", width=150, height=100, dpi=700, units="mm")
#  rm(data_d)

## ---- eval = FALSE------------------------------------------------------------
#  # Parameters
#  d <- ncol(yobs)
#  m <- d # prior parameter
#  A <- matrix(0, nrow=d, ncol=d) # matrix in the prior
#  
#  iter <- 1000
#  lower <- 10 # burn-in
#  upper <- iter + 1
#  dim <- 2 # beta_21

## ---- eval = FALSE------------------------------------------------------------
#  # Bayesian robust multivariate linear regression DA algorithm
#  ## There is no missing value. In this situation, the mda_gibbs and dami_gibbs coincide.
#  ## cw_gamma default is gamma(2, 2)
#  resdag <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # da monotone
#  resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma, Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter) # daimf fully observed
#  
#  # plot histogram
#  sim <- data.frame(Model=rep(c("dag", "daimfg"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
#  ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## ---- eval = FALSE------------------------------------------------------------
#  # create monotone missing data
#  for (i in c(47:50)) {
#    yobs[i, 1] <- NA
#  }
#  
#  # implement
#  resdag <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # da monotone
#  resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma,
#  Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter) # daimf monotone and impute to fully observed
#  
#  # plot histogram
#  sim <- data.frame(Model=rep(c("dag", "daig"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
#  ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## ---- eval = FALSE------------------------------------------------------------
#  # impute missing response values
#  resdagy <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter, yfull=TRUE)
#  resdaigy <- dageneral(X, yobs, m, A, cw=cw_gamma,
#  Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter, yfull=TRUE) # daimf# da
#  id <- ncol(resdagy$samples)
#  # plot histogram
#  sim <- data.frame(Model=rep(c("dag", "daimfg"), each=(upper - lower + 1)), beta=c(resdagy$samples[lower:upper, id], resdaigy$samples[lower:upper, id]))
#  ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("yimputed") + ylab("Density")

## ---- eval = FALSE------------------------------------------------------------
#  # create non-monotone missing data
#  for (i in c(47:50)) {
#    yobs[i, 1] <- NA
#  }
#  
#  yobs[45, 2] <- NA
#  Ik <- !is.na(yobs)
#  Ik[45, 2] <- TRUE
#  
#  # implement
#  resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma, Ik=Ik, iter=iter) # dai
#  resdaifg <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # daif
#  
#  # plot histogram
#  sim <- data.frame(Model=rep(c("daig", "daifg"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
#  ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

