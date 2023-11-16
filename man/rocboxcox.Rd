\name{rocboxcox}
\alias{rocboxcox}

\title{
Comprehensive Box-Cox ROC analysis for a single marker
}

\description{
This function applies the Box-Cox transformation to provide a comprehensive ROC analysis that involves the AUC (and its CI), the maximized Youden index (and its CI), the optimized Youden based cutoff (and its CI), and joint confidence regions for the optimized pair of sensitivity and specificity. Cutoffs are back-transformed so that they lie on the original scale of the data.
}

\usage{
rocboxcox(marker, D, alpha, plots, printProgress = FALSE)
}

\arguments{
  \item{marker}{A vector of length n that contains the biomarker scores of all individuals.}
  \item{D}{A vector of length n that contains the true disease status of an individual. It is a binary vector containing 0 for the healthy/control individuals and 1 for the diseased individuals.}
  \item{alpha}{Nominal level used to calculate the confidence intervals. A common choice is 0.05.}
  \item{plots}{Valid inputs are "on" and "off". When set to "on", the function returns a comprehensive figure with the ROC estimate and several point estimates: AUC, Youden index, optimized Youden cutoff, and the Youden based sensitivty and specificity along with the corresponding marginal confidence intervals and the joint confidence region of the estimated sensitivity and specificity.}
  \item{printProgress}{Boolean. When set to TRUE, messages describing the progress of the bootstrapping will be printed to the console window. When set to FALSE, these messages are suppressed. Default value is FALSE.}
}

\value{
  \item{transx}{The Box-Cox transformed scores for the healthy.}
  \item{transy}{The Box-Cox transformed scores for the diseased.}
  \item{transformation.parameter}{The estimated Box-Cox transformation parameter (lambda).}
  \item{AUC}{The estimated area under the Box-Cox based ROC curve.}
  \item{AUCCI}{The (1-a)100\% CI for the AUC. This CI is based on probit transforming the AUC estimate, finding the CI on the real line, and then back-transforming its endpoints to the ROC space.}
  \item{pvalueAUC}{The corresponding p-value for the AUC estimate.}
  \item{J}{The maximized Youden index.}
  \item{JCI}{The corresponding CI for maximized Youden index. For this CI consider we consider the probit transformation and then back-transform its endpoints to derive a 95\% CI for the Youden index itself.}
  \item{pvalueJ}{The corresponding two-tailed p-value for the maximized Youden index.}
  \item{Sens}{The sensitivity that corresponds to the Youden based optimized cutoff.}
  \item{CImarginalSens}{the marginal (1-a)100\% CI for the sensitivity that corresponds to the Youden based optimized cutoff. This is derived by first employing the probit transformation, finding a CI on the real line, and then back-transforming its endpoints to the ROC space.}
  \item{Spec}{The Specificity that corresponds to the Youden based optimized cutoff.}
  \item{CImarginalSpec}{The marginal (1-a)100\% CI for the specificity that corresponds to the Youden based optimized cutoff. This is derived by first employing the probit transformation, finding a CI on the real line, and then backtransforming its endpoints to the ROC space.}
  \item{cutoff}{The Youden-based optimized cutoff.}
  \item{CIcutoff}{The (1-a)100\% CI for the Youden-based optimized cutoff. This is based on the bootstrap. It involves using the Box-Cox transformation for every bootstrap iteration and then using the inverse Box-Cox transformation to obtain the cutoff on its original scale.}
  \item{areaegg}{The area of the (1-a)100\% egg-shaped joint confidence region that refers to the optimized pair of sensitivity and specificity. This takes into account the fact that the estimated sensitivity and specificity are correlated as opposed to the corresponding rectangular area that ignores this.}
  \item{arearect}{The area of the (1-a)100\% rectangular joint confidence region that refers to the optimized pair of sensitivity and specificity. This ignores the correlation of the optimized sensitivity and specificity and tends to yield a larger area compared to the one of the egg-shaped region.}
  \item{mxlam}{The mean of the marker scores of the healthy group after the Box-Cox transformation.}
  \item{sxlam}{The standard deviation of the marker scores of the healthy group after the Box-Cox transformation.}
  \item{mylam}{The mean of the marker scores of the diseased group after the Box-Cox transformation.}
  \item{sylam}{The standard deviation of the marker scores of the diseased group after the Box-Cox transformation.}
  \item{results}{A table that provides some indicative results: the AUC, the J (maximized Youden index), the estimated cutoff, the sensitivity, and the specificity along with their marginal CIs.}
  \item{roc}{A function of the estimated Box-Cox ROC curve. You can use this to simply request TPR values for given FPR values.}
}

\references{
Bantis LE, Nakas CT, Reiser B. (2014). Construction of confidence regions in the ROC space after the estimation of the optimal Youden index-based cut-off point. Biometrics, 70(1):212-223. doi: https://doi.org/10.1111/biom.12107

Box GEP, Cox DR. (1964). An Analysis of Transformations. Journal of the Royal Statistical Society. 26(2):211-252. https://www.jstor.org/stable/2984418
}

\author{
Leonidas Bantis
}

\examples{
set.seed(123)
x=rgamma(100, shape=2, rate = 8)
y=rgamma(100, shape=2, rate = 4)
scores=c(x,y)
D=c(pracma::zeros(1,100), pracma::ones(1,100))
out=rocboxcox(marker=scores,D, 0.05, plots="on")
}
