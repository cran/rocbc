\name{threerocs2}
\alias{threerocs2}

\title{
Provides visual comparison of three ROC estimation methods (two-marker version)
}

\description{
This function provides a visual comparison of the Empirical ROC, the Box-Cox ROC, and the Metz binormal semi-parametric estimator of the ROC curve. It also computes the AUC for the curve corresponding to each method.
}

\usage{
threerocs2(marker1, marker2, D, plots)
}

\arguments{
  \item{marker1}{A vector of length n that contains the biomarker scores of all individuals for the first marker.}
  \item{marker2}{A vector of length n that contains the biomarker scores of all individuals for the second marker.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", the user gets a single plot containing the estimated ROC curves using the Empirical, Box-Cox, and Metz methods for each of the two provided markers.}
}

\value{
  \item{AUC_Empirical1}{The AUC of the empirical ROC curve for the first marker.}
  \item{AUC_Empirical2}{The AUC of the empirical ROC curve for the second marker.}
  \item{AUC_Metz1}{The AUC of the Metz binormal curve (as calculated by MRMCaov package using the "binormal" option) for the first marker.}
  \item{AUC_Metz2}{The AUC of the Metz binormal curve (as calculated by MRMCaov package using the "binormal" option) for the second marker.}
  \item{AUC_BoxCox1}{The AUC of the Box-Cox based ROC curve for the first marker.}
  \item{AUC_BoxCox2}{The AUC of the Box-Cox based ROC curve for the second marker.}
}

\references{
Box GEP, Cox DR. (1964). An Analysis of Transformations. Journal of the Royal Statistical Society. 26(2):211-252. https://www.jstor.org/stable/2984418

Smith BJ, Hillis SL, Pesce LL (2023). MCMCaov: Multi-Reader Multi-Case Analysis of Variance. R package version 0.3.0, https://github.com/brian-j-smith/MRMCaov.

Smith BJ, Hillis SL (2020). “Multi-reader multi-case analysis of variance software for diagnostic performance comparison of imaging modalities.” In Samuelson F, Taylor-Phillips S (eds.), Proceedings of SPIE 11316, Medical Imaging 2020: Image Perception, Observer Performance, and Technology Assessment, 113160K. doi:10.1117/12.2549075, https://pubmed.ncbi.nlm.nih.gov/32351258.

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

out=threerocs2(marker1, marker2, D, plots = "on")
summary(out)
}
