% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesian.R
\name{bSignedRankTest}
\alias{bSignedRankTest}
\title{Bayesian equivalent to Wilcoxon's signed-rank test}
\usage{
bSignedRankTest(
  x,
  y = NULL,
  s = 0.5,
  z0 = 0,
  rope = c(-0.01, 0.01),
  nsim = 1000,
  seed = as.numeric(Sys.time())
)
}
\arguments{
\item{x}{First sample}

\item{y}{Second sample (if not provided, x is assumed to be the difference)}

\item{s}{Scale parameter of the prior Dirichlet Process. The default value is set to 0.5}

\item{z0}{Position of the pseudo-observation associated to the prior Dirichlet Process. The default value is set to 0 (inside the rope)}

\item{rope}{Interval for the difference considered as "irrelevant"}

\item{nsim}{Number of samples used to estimate the posterior distribution}

\item{seed}{Optional parameter used to fix the random seed}
}
\value{
A list with the following elements: 
\item{\code{method}}{a string with the name of the method used}
\item{\code{posterior.probabilities}}{a vector with the left, rope and right probabilities}
\item{\code{approximated}}{a logical value, \code{TRUE} if the posterior distribution is approximated (sampled) and \code{FALSE} if it is exact}
\item{\code{parameters}}{parameters used by the method}
\item{\code{posterior}}{Sampled probabilities (see details)}
}
\description{
Implementation of the Bayesian version of the signed-rank test presented in Benavoli \emph{et al.} 2017
}
\details{
The results includes the typical information relative to the three 
areas of the posterior density (left, right and rope probabilities), but also  
the result of the simulation used to estimate the probabilities.

The posterior contains the sampled probabilities in a matrix where each colums corresponds
to the sampled (posterior) probability of the measure (z), falling in that particular area (left to the rope, inside the rope and right to the rope). Conversely, the posterior probabilities refer to the probability of each region having the highest probability

As for the prior parameters, they are set to the default values indicated in Benavoli \emph{et al.} 2017 and you should not modify the unless you know what you are doing.
}
\examples{
sample1 <- rnorm(25, 1, 1)
sample2 <- rnorm(25, 1.2, 1)
results <- bSignedRankTest (x=sample1, y=sample2, s=0.5, z0=0)
res$posterior.probabilities

}
\references{
A. Benavoli, G. Corani, J. Demsar, M. Zaffalon (2017) Time for a Change: a Tutorial
for Comparing Multiple Classifiers Through Bayesian Analysis. \emph{Journal of Machine Learning Research}, 18, 1-36.
}
