\name{rocboxcoxCI}
\alias{rocboxcoxCI}

\title{
Inference around the sensitivity at a given specificity (and vice versa) for a single Box-Cox transformed biomarker
}

\description{
This function applies the Box-Cox transformation and provides the ability to do inferences for the sensitivity at a given specificity (and vice versa).
}

\usage{
rocboxcoxCI(marker, D, givenSP, givenSE, alpha, plots)
}

\arguments{
  \item{marker}{A vector of length n that contains the biomarker scores of all individuals.}
  \item{D}{A vector of length n that contains the true disease status of an individual, where 0 denotes a healthy/control individual, and 1 denotes a diseased individual.}
  \item{givenSP}{A vector of specificity values that the user wants to fix/set, at which the sensitivity is to be estimated. In this case, the ‘givenSE’ argument needs to be set to NA.}
  \item{givenSE}{A vector of sensitivity values that the user want to fix/set, at which the specificity is to be estimated. In this case, the ‘givenSP’ argument needs to be set to NA.}
  \item{alpha}{Nominal level used to calculate the confidence intervals. A common choice is 0.05.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", it returns both (1) the Box-Cox based ROC plot along with pointwise 95\% confidence intervals for the full spectrum of FPRs and (2) a second plot that visualizes the confidence intervals at the given sensitivities or specificities.}
}

\value{
  \item{SPandCIs}{The specificity values and the CIs around them.}
  \item{SEandCIs}{The sensitivity values and the CIs around them.}
  \item{SEvalues}{The sensitivity values provided by the user at which the specificity was calculated. If the user did not provide any sensitivity values, this argument should be set to NA.}
  \item{SPvalues}{The specificity values provided by the user at which the sensitivity was calculated. If the user did not provide any specificity values, this argument should be set to NA.}
}

\references{
Bantis LE, Feng Z. (2016). Comparison of two correlated ROC curves at a given specificity or sensitivity level. Statistics in Medicine, 35(24):4352-4367. https://doi.org/10.1002/sim.7008

Box GEP, Cox DR. (1964). An Analysis of Transformations. Journal of the Royal Statistical Society. 26(2):211-252. https://www.jstor.org/stable/2984418
}

\author{
Leonidas Bantis
}

\examples{
set.seed(123)
x <- rgamma(100, shape=2, rate = 8)
y <- rgamma(100, shape=2, rate = 4)
scores <- c(x,y)
D=c(pracma::zeros(1,100), pracma::ones(1,100))
givenSP=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
givenSE=NA
out=rocboxcoxCI(marker=scores, D, givenSP=givenSP, givenSE=NA, alpha=0.05, plots="on")
}
