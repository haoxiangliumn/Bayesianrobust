\name{cw}
\alias{cw}
\alias{cw_gig}
\alias{cw_gamma}
\alias{cw_constant}
\alias{cw_mass}
\alias{cw_gig0}
\alias{cw_gamma0}
\alias{cw_mass0}
\title{Conditional weight distributions}
\description{
Function to draw sample weight from the Conditional weight distributions in the DA/DAI algorithms I step.
}
\usage{
cw_gig(di, ri,...)
cw_gamma(di, ri, ...)
cw_constant(di, ri)
cw_mass(di, ri, ...) 
cw_gig0(di, ri, a=1, b=1, p=-0.5)
cw_gamma0(di, ri, a=2, b=2)
cw_mass0(di, ri, gamma=0.6, lambda=0.8) 
}
\arguments{
  \item{di}{integer. The number of components observed in the response Y_i}
  \item{ri}{non-negative real number.}
  \item{...}{parameters in the mixing distribution.}
  \item{a, b, p, gamma, lambda}{parameters in the mixing distribution.}
}
\details{
Function 
\code{cw(di,ri,...)} draws samples from the conditional mixing distribution, where
\deqn{P(w) \propto w^{\code{di}/2}\exp(-\code{ri}*w/2)P_{mix}(dw).}

\eqn{P_{mix}(dw)} is the mixing distribution for robust regression. A well-known example is that when \eqn{P_{mix}} follows gamma\eqn{(v/2, v/2)}, the errors follow multivariate \eqn{t(v)} distribution.
See \url{http://arxiv.org/abs/2212.01712} for more details.
%\code{...} Parameters in the mixing distribution. 
 
\code{cw_gig(di,ri,...):}{\eqn{P_{mix}} follows GIG(a, b, p), and a=1, b=1, p=-0.5 by default. It calls \code{cw_gig0(di, ri, a=1, b=1, p=-0.5)} function.}

\code{cw_gamma(di,ri,...):}{\eqn{P_{mix}} follows gamma(a, b), and a=2, b=2 by default.  It calls \code{cw_gamma0(di, ri, a=2, b=2)} function.}

\code{cw_constant(di,ri,...):}{\eqn{P_{mix}} is a constant (w=1).} 

\code{cw_mass(di,ri,...):}{\eqn{P_{mix}} follows \code{w=lambda} with probability gamma; \code{w = 1} with probability \code{1-gamma}, and \code{gamma=0.6}, \code{lambda=0.8} by default. It calls \code{cw_mass0(di, ri, gamma=0.6, lambda=0.8)} function.}


\code{cw_gamma0(di, ri, a=2, b=2): } \code{a, b} are the parameters in the mixing distribution with \code{rgamma(n = 1, shape = a, rate = b)}.

\code{cw_gig0(di, ri, a=1, b=1, p=-0.5):} \code{a, b, p} are the parameters in the mixing distribution with \code{rgig(n=1, chi=b, psi=a, lambda=p)}.

\code{cw_mass0(di, ri, gamma=0.6, lambda=0.8): } \code{gamma, lambda} are the parameters in the mixing distribution, and the distribution \eqn{P_{mix}} follows \code{w=lambda} with probability gamma; \code{w = 1} with probability \code{1-gamma}.

%To effectively implement the algorithms, the mixing distribution \eqn{P_{mix}(dw)} should satisfy
%\deqn{\int_{0}^\infty w^{d/2}P_{mix}(dw) <\infty,} where d is the number of response components.
}
\value{
Positive real number. Sample weight from the conditional weight distribution. 
}

\examples{
# P_{mix}(dw) follows gamma(a, b), then
cw_gamma0 <- function(di, ri, a=2, b=2) {
  return(rgamma(n=1, shape=di / 2 + a, rate=ri / 2 + b))
}
# cw function:
cw_gamma <- function(di, ri, ...) {
  #   ...: Mixing distribution parameters. a and b in gamma(a=2, b=2)
  return(cw_gamma0(di=di, ri=ri, ...))
}
}
\keyword{DA}