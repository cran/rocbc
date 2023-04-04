\name{comparebcAUC}
\alias{comparebcAUC}

\title{
Performs inference to compare the AUCs of two correlated (or uncorrelated) Box-Cox transformed biomarkers at a given nominal level
}

\description{
This function provides a comparison of two correlated markers in terms of their AUCs (areas under the Box-Cox based ROC curves). Markers measurements are assumed to be taken on the same individuals for both markers.
}

\usage{
comparebcAUC(marker1, marker2, D, alpha, plots)
}

\arguments{
  \item{marker1}{A vector of length n that contains the biomarker scores of all individuals for the first marker.}
  \item{marker2}{A vector of length n that contains the biomarker scores of all individuals for the second marker.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{alpha}{Nominal level used to calculate the confidence intervals. A common choice is 0.05.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on" it returns the Box-Cox based ROC along with informative information about the two AUCs in the legend of the plot.}
}

\value{
  \item{resultstable}{A summary table of the comparison that contains the AUC of each marker, the p-value of the difference using the probit transformation, the p-value of the difference, and the confidence interval of the difference.}
  \item{AUCmarker1}{The AUC of the first marker.}
  \item{AUCmarker2}{The AUC of the second marker.}
  \item{pvalue_probit_difference}{The p-value for the comparison of the AUCs. It employs the probit transformation that has greater power than when not using it (as in pvalue_difference).}
  \item{pvalue_difference}{The p-value for the comparisons of the AUC (without the probit transformation). Simulations have shown that this is inferior to the ‘pvalue_probit_difference’.}
  \item{CI_difference}{The confidence interval for the difference of the AUCs.}
  \item{roc1}{A function that refers to the ROC of the first marker. It allows the user to feed in FPR values and the corresponding TPR values.}
  \item{roc2}{A function that refers to the ROC of the first marker. It allows the user to feed in FPR values and the corresponding TPR values.}
  \item{transx1}{The Box-Cox transformed scores for the first marker and the healthy group.}
  \item{transy1}{The Box-Cox transformed scores for the first marker and the diseased group.}
  \item{transformation.parameter.1}{The estimated Box-Cox transformation parameter (lambda) for marker 1.}
  \item{transx2}{The Box-Cox transformed scores for the second marker and the healthy group.}
  \item{transy2}{The Box-Cox transformed scores for the second marker and the diseased group.}
  \item{transformation.parameter.2}{The estimated Box-Cox transformation parameter (lambda) for marker 2.}
}

\references{
Bantis LE, Feng Z. (2016). Comparison of two correlated ROC curves at a given specificity or sensitivity level. Statistics in Medicine, 35(24):4352-4367. https://doi.org/10.1002/sim.7008
}

\author{
Leonidas Bantis
}

\examples{
#GENERATE SOME BIVARIATE DATA===

set.seed(123)

nx <- 100
Sx <- matrix(c(1,   0.5,
               0.5,  1),
             nrow = 2, ncol = 2)

mux <- c(X = 10, Y = 12)
X=mvtnorm::rmvnorm(nx, mean = mux, sigma = Sx)

ny <- 100
Sy <- matrix(c(1.1,   0.6,
               0.6,  1.1),
             nrow = 2, ncol = 2)

muy <- c(X = 11, Y = 13.7)
Y=mvtnorm::rmvnorm(ny, mean = muy, sigma = Sy)

dx=pracma::zeros(nx,1)
dy=pracma::ones(ny,1)

markers=rbind(X,Y);
marker1=markers[,1]
marker2=markers[,2]
D=c(rbind(dx,dy))

#===COMPARE THE AUCs of Marker 1 vs Marker 2

out=comparebcAUC(marker1, marker2, D, alpha=0.05,  plots="on")
}
