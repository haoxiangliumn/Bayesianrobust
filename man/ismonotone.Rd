\name{ismonotone}
\alias{ismonotone}
\title{Monotone data}
\description{
Check whether the data is monotone or not, and permute the rows and columns of the data to the standard monotone format.
}
\usage{
ismonotone(yobs)
}
\arguments{
  \item{yobs}{ matrix. Response matrix (with missing values NA).}
}
\value{
Returns a list.
\item{monotone}{logic. \code{TRUE} if yobs is monotone, and \code{FALSE} otherwise.}

\item{ynew}{ matrix. New \code{yobs} in the standard monotone format or closest to the monotone format.}

\item{ynew_column}{vector. \code{ynew} (permuted) original column location in \code{yobs}.}

\item{ynew_row}{\code{ynew} (permuted) original row location in \code{yobs}.
}}



\examples{
# 1 not monotone
a = c(NA, 1,4,NA,5,9)
b = c(2, NA, 3,5,6,NA)
c = c(2,4,6,NA, 5,6)
data = data.frame(a,b,c)
ismonotone(data)

# 2 monotone
a = c(4, 1,4,NA,5,9)
b = c(2, NA, 3,NA,6,NA)
c = c(2,4,6,NA, 5,6)
data = data.frame(a,b,c)
ismonotone(data)
}
