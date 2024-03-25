\name{threerocs}
\alias{threerocs}

\title{
Provides visual comparison of three ROC estimation methods
}

\description{
This function provides a visual comparison of the Empirical ROC, the Box-Cox ROC, and the Metz binormal semi-parametric estimator of the ROC curve. It also computes the AUC for the curve corresponding to each method.
}

\usage{
threerocs(marker, D, plots)
}

\arguments{
  \item{marker}{A vector of length n that contains the biomarker scores of all individuals.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", the user gets the histograms of the biomarker for both the healthy and the diseased group before and after the Box-Cox transformation. In addition, all four corresponding qq-plots are provided.}
}

\value{
  \item{AUC_Empirical}{The AUC of the empirical ROC curve.}
  \item{AUC_Metz}{The AUC of the Metz binormal curve (as calculated by MRMCaov package using the "binormal" option).}
  \item{AUC_BoxCox}{The AUC of the Box-Cox based ROC curve.}
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
x <- rgamma(100, shape=2, rate = 8) # generates biomarker data from a gamma
                                 # distribution for the healthy group.
y <- rgamma(100, shape=2, rate = 4) # generates biomarker data from a gamma
                                 # distribution for the diseased group.
scores <- c(x,y)
D=c(pracma::zeros(1,100), pracma::ones(1,100))
out=threerocs(marker=scores, D, plots="on")
summary(out)
}
