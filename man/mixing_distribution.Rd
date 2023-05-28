\name{cw}
\alias{cw}
\alias{cw_gig}
\alias{cw_gamma}
\alias{cw_constant}
\alias{cw_mass}
\title{Conditional weight distributions}
\description{
Function to draw sample weight from the Conditional weight distributions in the DA/DAI algorithms I step.
}
\usage{
cw_constant(di, ri)
cw_gig(di, ri, a=1, b=1, p=-0.5)
cw_gamma(di, ri, a=2, b=2)
cw_mass(di, ri, gamma=0.6, lambda=0.8) 
}
\arguments{
  \item{di}{integer. The number of components observed in the response Y_i}
  \item{ri}{non-negative real number.}
  \item{a, b, p, gamma, lambda}{parameters in the mixing distribution.}
}
\details{
Function 
\code{cw(di,ri,...)} draws samples from the conditional mixing distribution, where
\deqn{P(w) \propto w^{\code{di}/2}\exp(-\code{ri}*w/2)P_{mix}(dw).}

\eqn{P_{mix}(dw)} is the mixing distribution for robust regression. A well-known example is that when \eqn{P_{mix}} follows gamma\eqn{(v/2, v/2)}, the errors follow multivariate \eqn{t(v)} distribution.
See \url{http://arxiv.org/abs/2212.01712} for more details.
%\code{...} Parameters in the mixing distribution. 
 
\code{cw_gig(di, ri, a=1, b=1, p=-0.5)):}\eqn{P_{mix}} follows GIG(a, b, p), and a=1, b=1, p=-0.5 by default. 

\code{cw_gamma(di, ri, a=2, b=2):}\eqn{P_{mix}} follows gamma(a, b), and a=2, b=2 by default.

\code{cw_constant(di,ri):}{\eqn{P_{mix}} is a constant (w=1).} 

\code{cw_mass(di, ri, gamma=0.6, lambda=0.8):}{\eqn{P_{mix}} follows \code{w=lambda} with probability gamma; \code{w = 1} with probability \code{1-gamma}, and \code{gamma=0.6}, \code{lambda=0.8} by default. }




%To effectively implement the algorithms, the mixing distribution \eqn{P_{mix}(dw)} should satisfy
%\deqn{\int_{0}^\infty w^{d/2}P_{mix}(dw) <\infty,} where d is the number of response components.
}
\value{
Positive real number. Sample weight from the conditional weight distribution. 
}

\examples{
# cw_gamma: P_{mix} follows gamma(a=2, b=2), then
cw_gamma <- function(di, ri, a=2, b=2) {
  return(rgamma(n=1, shape=di / 2 + a, rate=ri / 2 + b))
}
}
\keyword{DA}