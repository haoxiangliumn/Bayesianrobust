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
\item{samples}{ matrix. MCMC samples with \code{iter+1} rows and \code{p*d+d*d} columns. Each row is one iteration of the Markov chain (from time 0 to time \code{iter}). The sampled matrices are vectorized by stacking the columns. The \code{1}st to \code{p*d}th columns are the sampled \eqn{\beta} coefficients, and the \code{p*d+1}th to \code{p*d+d*d}th coulmns are the sampled \eqn{\Sigma} coeffficients. If \code{yfull} is \code{TRUE}, the rest of elements in each row are imputed missing components of \code{yobs}.}

\item{time}{ positive real number. Time to sample \code{iter} number of samples.}

It will print the corresponding DA/DAI algorithms used. See \code{\link{DA_functions}} for introductions of these algorithms.
}

\details{
Data augmentation algorithms to sample from the posterior in Bayesian robust multivariate linear regression with incomplete response. 

For the full model and algorithm details, see \url{http://arxiv.org/abs/2212.01712}. 

The \code{cw} function generates samples from the conditional weight distribution, see \code{\link{cw}} page for the choices of distributions provided and how to define such functions manually. 

%{monotone:} 
%Assume that \code{yobs} have \code{n} observations and each has \code{d} components. \code{yobs} is monotone if for the \code{i}th observation, \code{i = 1,...,n}, with \code{j_i} the first observed component, the \code{(j_i+1)}th to \code{d}th elements are all observed, and \code{j_i} is no smaller than \code{j_{i+1}}. Note that in practice, there are cases where the observed missing structure can be rearranged to become monotone by permuting the rows and columns of the response matrix.
%If there are no missing data, the corresponding missing structure is monotone.

{Prior:} We assume the following commonly used groups of priors: \deqn{\text{prior}(\beta,\Sigma)\propto |\Sigma|^{-(m+1)/2} \exp(-0.5 \text{tr}(\Sigma^{-1} A))}
The default values of \eqn{m} and \eqn{\Sigma} are the independence Jeffrey's prior.

{Monotone:} See the paper for monotone missing structure definition and standard format. See \code{\link{monotone}} function for monotone checking and permutation into standard format. 

\code{Ik}: If \code{yobs} is monotone, then by default, there is no imputation in the DA I step. 
The missing structure \code{Ik} represents the missing structure after mandatory imputation in the DAI I step. 
If \code{yobs} is not monotone, the default \code{Ik} is fully observed, i.e.,  all elements in \code{Ik} are \code{TRUE}. 
Clearly, elements observed in \code{yobs} should also be non-missing under \code{Ik}.

%Let \code{I} be the missing structure of \code{yobs}, i.e., \code{I_{ij}=1} if \code{yobs_{ij}} is not \code{NA}, and  
%\code{I_{ij}=0} if \code{yobs_{ij}} is \code{NA}. Let \code{Ik} be another missing structure. \code{Ik} is larger than \code{I} if for all nonzero elements in \code{I}, the corresponding elements in \code{Ik} are also nonzero.

%\code{Ik} should be larger than \code{I}, where \code{I} is the missing pattern of \code{yobs}.

\code{samples}: Let \code{p} be the number of predictors (\code{ncol(X)}), and \code{d} be the number of response componnets (\code{ncol(yobs)}).  \code{samples}
 stores \code{iter+1} MCMC samples. The sampled matrices are vectorized by stacking the columns. 
 For the ith sample, \eqn{\beta(i)} is \code{as.matrix(samples[i+1, 1:p*d], nrow=p, ncol=d)}, and \eqn{\Sigma(i)} is \code{as.matrix(samples[i+1, (p*d+1):p*d+d^2], nrow=d, ncol=d)}. If \code{yfull} is \code{TRUE}, the rest of the elements in each row are imputed missing components of \code{yobs}. 


%To effectively perform the algorithm, there should exist \code{(X', yobs')}, a subset of \code{(X, yobs)}, satisfying assumption H1 and H2, where \code{yobs}' is monotone and \code{X}' is the corresponding observed features. See \code{paper} for the details.

%H2: The mixing distribution \eqn{P_{mix}(dw)} should satisfy
%\deqn{\int_{0}^\infty w^{d/2}P_{mix}(dw) <\infty,} where d is the number of response components.

%A special case of H1 is to pick fully observed part \code{yobs}'. Then  H1 is \deqn{r(yobs' : X') = p + d, n > p + 2d - m - 1,} where \eqn{yobs':X'} is a matrix combining \eqn{X'} on the right side of \eqn{yobs'}, and \code{r} is the rank of the matrix.
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
\references{
Li, Haoxiang, Qian Qin, and Galin L. Jones. "Convergence Analysis of Data Augmentation Algorithms for Bayesian Robust Multivariate Linear Regression with Incomplete Data." arXiv preprint arXiv:2212.01712 (2022).
}