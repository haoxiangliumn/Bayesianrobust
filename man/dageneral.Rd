\name{DA_general}
\alias{dageneral}
\title{Data augmentation samplers}
\description{
Data augmentation algorithms for Bayesian robust multivariate linear regression with incomplete response.
}
\usage{
dageneral(X,yobs,m=ncol(yobs),A=diag(0,ncol(yobs)),cw=cw_gamma,Ik,iter=100,yfull=FALSE)
}
\arguments{
  \item{X}{ matrix. Feature matrix.}
  \item{yobs}{ matrix. Response matrix (with missing values NA).}
  \item{m}{ real number. Default \code{ncol(yobs)}. Prior parameter. See details.} 
  \item{A}{ matrix. Default \code{diag(0, ncol(yobs))}. Positive semidefinite matrix in the prior. See details.}   
  \item{cw}{function. Default \code{cw_gamma}. Function to generate samples from conditional weight distribution. See details.}
  \item{Ik}{optional. Matrix and all elements should be logical values.  See details.}
  \item{iter}{ positive integer. Default \code{100}. Number of MCMC samples}
  \item{yfull}{ logical, default \code{FALSE}. \code{TRUE}: impute all missing componnents in the response and output them. \code{FALSE}: do not output or calculate (if non-necessary) the missing components in the response.
}
}
\value{
Returns a list with two elements.
\item{samples}{ matrix. MCMC samples with \code{iter+1} rows and \code{p*d+d*d} columns. Each row is one iteration of the sampler (from \code{t=0} to \code{t=iter}). The sampled matrices are vectorized by stacking the columns. The \code{1}st to \code{p}th columns are the sampled \eqn{\beta} coefficients, and the \code{p+1}th to \code{p*d+d*d}th coulmns are the sampled \eqn{\Sigma} coeffficients. If \code{yfull} is \code{TRUE}, the rest of elements in each row are imputed missing componnets of \code{yobs}. See details.}

\item{time}{ positive real number. Time to sample \code{iter} number of samples.}

It will print the corresponding DA/DAI algorithms used. See \code{\link{DA_functions}} for the functions implementing these algorithms.
}

\details{
Data augmentation algorithms to sample from the posterior in Bayesian robust multivariate linear regression with incomplete response. For the full model and algorithm details. See \eqn{paper}.

\code{cw} function generates samples from the conditional weight distribution, see \code{\link{cw}} page for how to define such functions manually and the choices of distributions provided. 

{monotone:} 
Assume that \code{yobs} have \code{n} observations and each has \code{d} components. \code{yobs} is monotone if for the \code{i}th observation, \code{i = 1,...,n}, with \code{j_i} the first observed component, the \code{(j_i+1)}th to \code{d}th elements are all observed, and \code{j_i} is no smaller than \code{j_{i+1}}. Note that in practice, there are cases where the observed missing structure can be rearranged to become monotone by permuting the rows and columns of the response matrix.
If there are no missing data, the corresponding missing structure is monotone.

{prior:} We assume the following commonly used groups of priors: \deqn{prior(\beta,\Sigma)\propto |\Sigma|^{-(m+1)/2} exp(-0.5 tr(\Sigma^{-1} A))}
The default choices of \eqn{m} and \eqn{\Sigma} are the independence Jeffrey's prior.

\code{Ik} represents the missing structure after mandatory imputation in the data augmentation I step. If \code{yobs} is monotone, then by default, there is no imputation in the DA I step. If \code{yobs} is not monotone, the default \code{Ik} is fully observed, i.e., \code{matrix(TRUE, nrow=nrow(yobs), ncol=ncol(yobs))}. 

Let \code{I} be the missing structure of \code{yobs}, i.e., \code{I_{ij}=1} if \code{yobs_{ij}} is not \code{NA}, and  
\code{I_{ij}=0} if \code{yobs_{ij}} is \code{NA}. Let \code{Ik} be another missing structure. \code{Ik} is larger than \code{I} if for all nonzero elements in \code{I}, the corresponding elements in \code{Ik} are also nonzero.

\code{Ik} should be larger than \code{I}, where \code{I} is the missing pattern of \code{yobs}.

\code{samples}: Let \code{p} be the number of predictors (\code{ncol(X)}), and \code{d} be the number of response componnets (\code{ncol(yobs)}).  \code{samples}
 stores \code{iter+1} MCMC samples. The sampled matrices are vectorized by stacking the columns. For the ith sample, \eqn{\beta} is \code{as.matrix(samples[i+1, 1:p*d], nrow=p, ncol=d)}, and \eqn{\Sigma} is \code{as.matrix(samples[i+1, (p*d+1):p*d+d^2], nrow=d, ncol=d)}. If \code{yfull} is \code{TRUE}, the rest of the elements in each row are imputed missing components of \code{yobs}, and they are \code{!I} elemets of \code{yobs}. 


To effectively perform the algorithm, there should exist \code{(X', yobs')}, a subset of \code{(X, yobs)}, satisfying assumption H1 and H2, where \code{yobs}' is monotone and \code{X}' is the corresponding observed features. See \code{paper} for the details.

H2: The mixing distribution \eqn{P_{mix}(dw)} should satisfy
\deqn{\int_{0}^\infty w^{d/2}P_{mix}(dw) <\infty,} where d is the number of response components.

A special case of H1 is to pick fully observed part \code{yobs}'. Then  H1 is \deqn{r(yobs' : X') = p + d, n > p + 2d - m - 1,} where \eqn{yobs':X'} is a matrix combining \eqn{X'} on the right side of \eqn{yobs'}, and \code{r} is the rank of the matrix.
 }

\examples{
set.seed(2)
library(Bayesianrobust)
library(ggplot2)
data(sim_data)
X <- cbind(1, sim_data[, ncol(sim_data)])
yobs <- sim_data[, 1:(ncol(sim_data) - 1)]
iter <- 1000

r1 <- dageneral(X, yobs, iter=iter) 
dim <- 2  # beta_21
r1$time # calculation time 
plot(density(r1$samples[, dim]), main="density")
}
\keyword{DA}