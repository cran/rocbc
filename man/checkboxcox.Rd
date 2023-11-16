\name{checkboxcox}
\alias{checkboxcox}

\title{
Tests whether Box-Cox is appropriate for the given dataset (one-marker version)
}

\description{
This function tests whether the Box-Cox transformation is able to achieve approximate normality for your data. That is, it will allow the user to know whether it is appropriate to use all the methods discussed later on in this package.
}

\usage{
checkboxcox(marker, D, plots, printShapiro = FALSE)
}

\arguments{
  \item{marker}{A vector of length n that contains the biomarker scores of all individuals.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", the user gets the histograms of the biomarker for both the healthy and the diseased group before and after the Box-Cox transformation. In addition, all four corresponding qq-plots are provided.}
  \item{printShapiro}{Boolean. When set to TRUE, the results of the Shapiro-Wilk test will be printed to the console. When set to FALSE, the results are suppressed. Default value is FALSE.}
}

\value{
  \item{res_shapiro}{A results table that contains the results of four Shapiro-Wilk tests for normality testing. Two of these refer to normality testing of the healthy and the diseased groups before the Box-Cox transformation, and the remaining two refer to the Box-Cox transformed biomarkers scores for the healthy and the diseased groups. Thus, this testing process produces four p-values. In addition, if the plots are set to 'on', then the output provides (1) the histograms of the biomarker for both the healthy and the diseased groups before and after the Box-Cox transformation, (2) all four corresponding qq-plots, and (3) a plot with the empirical ROC curve overplotted with the Box-Cox based ROC curve for visual comparison purposes.}
  \item{transformation.parameter}{The single transformation parameter, lambda, that is applied for both groups simultaneously.}
  \item{transx}{The Box-Cox transformed scores for the healthy.}
  \item{transy}{The Box-Cox transformed scores for the diseased.}
  \item{pval_x}{The p-value of the Shapiro Wilk test of normality for the healthy group (before the Box-Cox transformation).}
  \item{pval_y}{The p-value of the Shapiro Wilk test of normality for the diseased group (before the Box-Cox transformation).}
  \item{pval_transx}{The p-value of the Shapiro Wilk test of normality for the healthy group (after the Box-Cox transformation).}
  \item{pval_transy}{The p-value of the Shapiro Wilk test of normality for the diseased group (after the Box-Cox transformation).}
  \item{roc}{A function of the estimated Box-Cox ROC curve. You can use this to simply request TPR values for given FPR values.}
}

\references{
Bantis LE, Nakas CT, Reiser B. (2021). Statistical inference for the difference between two maximized Youden indices obtained from correlated biomarkers. Biometrical Journal, 63(6):1241-1253. https://doi.org/10.1002/bimj.202000128

Bantis LE, Nakas CT, Reiser B. (2018). Construction of confidence intervals for the maximum of the Youden index and the corresponding cutoff point of a continuous biomarker. Biometrical Journal, 61(1):138-156. https://doi.org/10.1002/bimj.201700107

Bantis LE, Feng Z. (2016). Comparison of two correlated ROC curves at a given specificity or sensitivity level. Statistics in Medicine, 35(24):4352-4367. https://doi.org/10.1002/sim.7008

Bantis LE, Nakas CT, Reiser B. (2014). Construction of confidence regions in the ROC space after the estimation of the optimal Youden index-based cut-off point. Biometrics, 70(1):212-223. https://doi.org/10.1111/biom.12107

Box GEP, Cox DR. (1964). An Analysis of Transformations. Journal of the Royal Statistical Society. 26(2):211-252. https://www.jstor.org/stable/2984418
}

\author{
Leonidas Bantis
}

\examples{
set.seed(123)
x <- rgamma(100, shape=2, rate = 8) # generates biomarker data from a gamma
                                 # distribution for the healthy group.
y <- rgamma(100, shape=2, rate = 4) # generates biomarker data from a gamma
                                 # distribution for the diseased group.
scores <- c(x,y)
D=c(pracma::zeros(1,100), pracma::ones(1,100))
out=checkboxcox(marker=scores, D, plots="on")
summary(out)
}
