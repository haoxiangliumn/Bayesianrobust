\name{DA_functions}
\alias{DA_functions}
\alias{daimf}
\alias{da}
\alias{daif}
\alias{dai}
\title{Data augmentation samplers introduction}
\description{
DA/DAI samplers for Bayesian robust multivariate linear regression with incomplete response.
}
\usage{
da(X,yobs,m=ncol(yobs),A=diag(0,ncol(yobs)),cw=cw_gamma,iter=100,yfull=FALSE)
daimf(X,yobs,m=ncol(yobs),A=diag(0,ncol(yobs)),cw=cw_gamma,iter=100,yfull=FALSE)
daif(X,yobs,m=ncol(yobs),A=diag(0,ncol(yobs)),cw=cw_gamma,iter=100,yfull=FALSE)
dai(X,yobs,m=ncol(yobs),A=diag(0,ncol(yobs)),cw=cw_gamma,Ik,iter=100,yfull=FALSE)
}
\arguments{
  \item{X}{ matrix. Feature matrix.}
  \item{yobs}{ matrix. Response matrix (with missing values NA).}
  \item{m}{ real number. Prior parameter. Default \code{ncol(yobs)}. } 
  \item{A}{ matrix. Positive semidefinite matrix in the prior. Default \code{diag(0, ncol(yobs))}.}   
  \item{cw}{function. Function to sample from the conditional weight distribution. Default \code{cw_gamma}. }
  \item{Ik}{matrix and all elements should be logical values. Missing structure matrix. Indicate the imputed data in the I step. Default fully observed.}
  \item{iter}{ positive integer. The length of the  Markov chain. Default \code{100}. }
  \item{yfull}{ logical. \code{TRUE}: impute all missing components in the response and output them. \code{FALSE}: do not output or calculate (if non-necessary) the missing components in the response. Default \code{FALSE}.
}
}
\value{
Returns a list with two elements.
\item{samples}{ matrix. MCMC samples with \code{iter+1} rows and \code{p*d+d*d} columns. Each row is one iteration of the Markov chain (from time 0 to time \code{iter}). The sampled matrices are vectorized by stacking the columns. The \code{1}st to \code{p*d}th columns are the sampled \eqn{\beta} coefficients, and the \code{p*d+1}th to \code{p*d+d*d}th columns are the sampled \eqn{\Sigma} coefficients. If \code{yfull} is \code{TRUE}, the rest of the elements in each row are imputed missing components of \code{yobs}.}

\item{time}{ positive real number. Time to sample \code{iter} number of samples.}

}
\details{
Use function \code{\link{dageneral}} to implement these algorithms. It will automatically find the best fit to sample from the posterior.
 
Function \code{da} implements the DA algorithm and requires \code{yobs} to be monotone.
 
Function \code{daimf} implements the DAI algorithm, requires \code{yobs} to be monotone and imputes data to be fully observed in the DAI I step.
  
Function \code{daif} implements the DAI algorithm and imputes data to be fully observed in the DAI I step.

Function \code{daif} implements the DAI algorithm.
}

\examples{
library(Bayesianrobust)
library(ggplot2)
data(sim_data)
X <- cbind(1, sim_data[, ncol(sim_data)])
yobs <- sim_data[, 1:(ncol(sim_data) - 1)]
iter <- 1000
r1 <- daimf(X, yobs, iter=iter) 
r2 <- da(X, yobs, iter=iter) 
r3 <- daif(X, yobs, iter=iter) 
r4 <- dai(X, yobs, iter=iter) 
dim <- 2  # Beta_{21}

# draw histogram for the (1-2) coefficient
burn <- 6 # burnin burn-1
sim <- data.frame(
Model=rep(c("daimf", "da", "daif", "dai"), each=(iter-burn+2)), 
beta=c(r1$samples[burn:(iter+1), dim], r2$samples[burn:(iter+1), dim], 
r3$samples[burn:(iter+1), dim], r4$samples[burn:(iter+1), dim]))
ggplot(sim, aes(x=beta, color=Model)) + geom_density() + xlab("Beta") + ylab("Density") 
}
