## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)
library(Bayesianrobust) # load our package

library(GeneralizedHyperbolic)
library(LaplacesDemon)
library(matlib)
library(MASS)
library(matrixNormal)
library(ggplot2)
# generate data
library(mvtnorm)
n <- 50
X <- cbind(1, as.matrix(rnorm(n, 0, 1)))
beta <- matrix(c(0.8, 0.4, 0.5, 0.3), nrow=2, ncol=2)
var <- matrix(c(1, 0.6, 0.6, 2), 2, 2)
error <- matrix(0, nrow=n, ncol=2)
w <- rgamma(n, shape=2, rate=2)
for (i in 1:n) {
  error[i, ] <- rmvnorm(1, rep(0, nrow(var)), var / w[i])
}
yobs <- X %*% beta + error
data <- cbind(yobs, X[, 2])
data_d <- cbind(X[, 2], yobs)
colnames(data_d) <- c("x", "yr1", "yr2")

# display data
library(GGally)
ggpairs(data.frame(data_d))
# ggsave("scatter.png", width=150, height=100, dpi=700, units="mm")
rm(data_d)

## -----------------------------------------------------------------------------
# Parameters
d <- ncol(yobs)
m <- d # prior parameter
A <- matrix(0, nrow=d, ncol=d) # matrix in the prior

iter <- 1000
lower <- 1 # burn-in
upper <- iter + 1
dim <- 2 # beta_21 

## -----------------------------------------------------------------------------
# Bayesian robust multivariate linear regression DA algorithm
## There is no missing value. In this situation, the mda_gibbs and dami_gibbs coincide.
## cw_gamma default is gamma(2, 2)
resdag <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # da monotone
resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma, Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter) # daimf fully observed

# plot histogram
sim <- data.frame(Model=rep(c("dag", "daimfg"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## -----------------------------------------------------------------------------
# create monotone missing data
for (i in c(47:50)) {
  yobs[i, 1] <- NA
}

# implement
resdag <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # da monotone
resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma,
Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter) # daimf monotone and impute to fully observed

# plot histogram
sim <- data.frame(Model=rep(c("dag", "daig"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

## -----------------------------------------------------------------------------
# impute missing response values
resdagy <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter, yfull=TRUE)
resdaigy <- dageneral(X, yobs, m, A, cw=cw_gamma,
Ik=matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs)),  iter=iter, yfull=TRUE) # daimf# da
id <- ncol(resdagy$samples)
# plot histogram
sim <- data.frame(Model=rep(c("dag", "daimfg"), each=(upper - lower + 1)), beta=c(resdagy$samples[lower:upper, id], resdaigy$samples[lower:upper, id]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("yimputed") + ylab("Density")

## -----------------------------------------------------------------------------
# create non-monotone missing data
for (i in c(47:50)) {
  yobs[i, 1] <- NA
}

yobs[45, 2] <- NA
Ik <- !is.na(yobs)
Ik[45, 2] <- TRUE

# implement
resdaig <- dageneral(X, yobs, m, A, cw=cw_gamma, Ik=Ik, iter=iter) # dai
resdaifg <- dageneral(X, yobs, m, A, cw=cw_gamma, iter=iter) # daif

# plot histogram
sim <- data.frame(Model=rep(c("daig", "daifg"), each=(upper - lower + 1)), beta=c(resdag$samples[lower:upper, dim], resdaig$samples[lower:upper, dim]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density")

