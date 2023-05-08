\name{creatinine_data}
\title{Creatinine data}
\description{
The data for the creatinine clearance from Shih and Weisberg (1986).   
}
\usage{
creatinine_data
}
\details{
This data is used in an illustrative example by Liu (1996).

The data are from a clinical trial on 34 male patients with body weight (WT) in kg, serum creatinine (SC) concentration in mg/deciliter, age (AGE) in years, and endogenous creatinine (CR) clearance. 
InSC = In(SC), InWT = In(WT), InCR = In(CR), In140_AGE = In(140-AGE).
Of the 34 male patients, 2 had no record in WT and 4 were missing SC. The missingness is assumed to be ignorable.
There is one influential observation patient 27.
}



\examples{
data(creatinine_data)
}

\references{
Liu, Chuanhai. "Bayesian robust multivariate linear regression with incomplete data." Journal of the American Statistical Association 91.435 (1996): 1219-1227.

Shih, Weichung J., and Sanford Weisberg. "Assessing influence in multiple linear regression with incomplete data." Technometrics 28.3 (1986): 231-239.

}
