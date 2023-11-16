\name{checkboxcox2}
\alias{checkboxcox2}

\title{
Tests whether Box-Cox is appropriate for the given dataset (two-marker version)
}

\description{
This function tests whether the Box-Cox transformation is able to achieve approximate normality for your data. That is, it will allow the user to know whether it is appropriate to use all the methods discussed later on in this package.
}

\usage{
checkboxcox2(marker1, marker2, D, plots, printShapiro = FALSE)
}

\arguments{
  \item{marker1}{A vector of length n that contains the biomarker scores of all individuals for the first marker.}
  \item{marker2}{A vector of length n that contains the biomarker scores of all individuals for the second marker.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", the user gets the histograms of the biomarker for both the healthy and the diseased group before and after the Box-Cox transformation for both marker1 and marker2. In addition, all eight corresponding qq-plots are provided.}
  \item{printShapiro}{Boolean. When set to TRUE, the results of the Shapiro-Wilk test will be printed to the console. When set to FALSE, the results are suppressed. Default value is FALSE.}
}

\value{
  \item{res_shapiro}{A results table that contains the results of eight Shapiro-Wilk tests for normality testing. Four of these refer to normality testing of the healthy and the diseased groups before the Box-Cox transformation, and the remaining four refer to the Box-Cox transformed biomarkers scores for the healthy and the diseased groups. Thus, this testing process produces eight p-values. In addition, if the plots are set to 'on', then the output provides (1) the histograms of the biomarker for both the healthy and the diseased groups before and after the Box-Cox transformation for both marker1 and marker2, (2) all eight corresponding qq-plots, and (3) two plots (one for marker1, one for marker2) with the empirical ROC curve overplotted with the Box-Cox based ROC curve for visual comparison purposes.}
  \item{transx1}{The Box-Cox transformed scores for the first marker and the healthy group.}
  \item{transy1}{The Box-Cox transformed scores for the first marker and the diseased group.}
  \item{transformation.parameter.1}{The estimated Box-Cox transformation parameter (lambda) for marker 1.}
  \item{transx2}{The Box-Cox transformed scores for the second marker and the healthy group.}
  \item{transy2}{The Box-Cox transformed scores for the second marker and the diseased group.}
  \item{transformation.parameter.2}{The estimated Box-Cox transformation parameter (lambda) for marker 2.}
  \item{pval_x_marker1}{The p-value of the Shapiro Wilk test of normality for the marker1 healthy group (before the Box-Cox transformation).}
  \item{pval_y_marker1}{The p-value of the Shapiro Wilk test of normality for the marker1 diseased group (before the Box-Cox transformation).}
  \item{pval_transx_marker1}{The p-value of the Shapiro Wilk test of normality for the marker1 healthy group (after the Box-Cox transformation).}
  \item{pval_transy_marker1}{The p-value of the Shapiro Wilk test of normality for the marker1 diseased group (after the Box-Cox transformation).}
  \item{pval_x_marker2}{The p-value of the Shapiro Wilk test of normality for the marker2 healthy group (before the Box-Cox transformation).}
  \item{pval_y_marker2}{The p-value of the Shapiro Wilk test of normality for the marker2 diseased group (before the Box-Cox transformation).}
  \item{pval_transx_marker2}{The p-value of the Shapiro Wilk test of normality for the marker2 healthy group (after the Box-Cox transformation).}
  \item{pval_transy_marker2}{The p-value of the Shapiro Wilk test of normality for the marker2 diseased group (after the Box-Cox transformation).}
  \item{roc1}{A function that refers to the ROC of the first marker. It allows the user to feed in FPR values and the corresponding TPR values.}
  \item{roc2}{A function that refers to the ROC of the first marker. It allows the user to feed in FPR values and the corresponding TPR values.}
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

nx <- 100
Sx <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)

mux <- c(X = 10, Y = 12)
X = mvtnorm::rmvnorm(nx, mean = mux, sigma = Sx)

ny <- 100
Sy <- matrix(c(1.1, 0.6, 0.6, 1.1), nrow = 2, ncol = 2)

muy <- c(X = 11, Y = 13.7)
Y = mvtnorm::rmvnorm(ny, mean = muy, sigma = Sy)

dx = pracma::zeros(nx,1)
dy = pracma::ones(ny,1)

markers = rbind(X,Y);
marker1 = markers[,1]
marker2 = markers[,2]
D = c(rbind(dx,dy))

out=checkboxcox2(marker1, marker2, D, plots = "on")
summary(out)
}
