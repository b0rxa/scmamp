getNemenyiCD <- function (alpha = 0.05, num.alg, num.problems) {
  # Auxiliar function to compute the critical difference for Nemenyi test
  # Args:
  #   alpha:        Alpha for the test
  #   num.alg:      Number of algorithms tested
  #   num.problems: Number of problems where the algorithms have been tested
  #
  # Returns:
  #   Corresponding critical difference
  #
  df <- num.alg * (num.problems - 1)
  qa <- qtukey(p=1 - alpha, nmeans=num.alg, df=df)/sqrt(2)
  cd <- qa * sqrt((num.alg * (num.alg + 1)) / (6 * num.problems))
  return(cd)
}


# EXPORTED FUNCTIONS -----------------------------------------------------------

#' @title Wilcoxon signed-rank est
#'
#' @description This function implements the paired Wilcoxon signed-rank test
#' @param x First sample
#' @param y Second sample
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Demsar (2006), page 7
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @examples
#' x <- rbeta(50, 2, 20)
#' y <- x + runif(50) * 0.2
#' wilcoxonSignedTest(x, y)

wilcoxonSignedTest <- function (x, y, ...) {
  if (length(x) != length(y)) {
    stop("This is a paired test, so the two vectors have to have the same length")
  }
  N  <- length(x)
  d  <- x-y
  # Compute the statistic based on the ordering of the differences. 
  # We need to assign the highest rank to the biggest value, so we change the 
  # sign of the absolute difference.
  o  <- rank(-abs(d), ties.method="average")
  rp <- sum(o[d>0]) + sum(o[d == 0])
  rn <- sum(o[d<0]) + sum(o[d == 0])
  t  <- min(rp, rn)
  
  num    <- t - 0.25 * (N * (N + 1))
  den    <- sqrt((N * (N + 1) * (2 * N + 1)) / 24)
  z      <- num / den
  pvalue <- pnorm(z)  
  
  names(t)  <- "T"
  method    <- "Wilcoxon Signed-Rank test"
  data.name <- deparse(substitute(data))
  
  # Create and return the htest object
  res <- list(statistic=t, p.value=pvalue, method=method, data.name=data.name )
  class(res) <- "htest"
  return(res)
}

#' @title ANOVA test for multiple comparisons
#'
#' @description This function performs F-test for K populations means
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to Test 22 in Kanji (2006).
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @examples
#' data(data_gh_2008)
#' anovaTest(data.gh.2008)
#' 

anovaTest <- function (data, ...){
  
  # Implemented according to Test 22 in 100 statistical tests
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  m.x <- colMeans(data)
  
  # Variance of the observations with respect to their own sample
  var.1 <- mean(apply(data, MARGIN=2, FUN=var))
  
  # Variance ofthe sample means with respect to the grand mean
  g.m   <- mean(unlist(data))
  var.2 <- sum(N * (m.x - g.m) ^ 2) / (k - 1)
  
  # Statistic that follows an F distribution with (k - 1); k * (N - 1) degrees of freedom
  F.stat <- var.2 / var.1 
  pvalue <- 1 - pf(F.stat, k - 1, k * (N - 1))
  
  # Build and return the htest object
  parameters    <- c(df1=k - 1, df2=k * (N - 1))
  names(F.stat) <- "F statistic"
  method        <- "F-test for K population means (Analysis of Variance)"
  data.name     <- deparse(substitute(data))
  res <- list(statistic=F.stat, p.value=pvalue, parameters=parameters, 
              method=method, data.name=data.name)
  class(res) <- "htest"
  return(res)
}


#' @title Get the ranking matrix.
#'
#' @description This function returns, given a matrix, the ranking of the colums in each row.
#' @param data The matrix to rank.
#' @param decreasing Logical value indicating whether the top ranked has to be the highest value or not.
#' @param ... Not used
#' @return A matrix containing the per-row rankings. In case of ties, the mean rank is obtained (e.g, if there is a tie between the 4th and the 5th column, both are assigned a mean rank of 4.5)
#' @examples
#' data(data_gh_2008)
#' rankMatrix(data.gh.2008)
#' 
rankMatrix <- function(data, decreasing=TRUE, ...){
  # The rank function is based on an increasing ordering. In case we need to
  # get the rank of the descreasing ordering, just rank -x instead of x
  if (decreasing){
    f <- function(x){
      rank(-x, ties.method="average")
    }
  } else {
    f <- function(x){
      rank(x, ties.method="average")
    }
  }
  
  rankings <- t(apply (data, MARGIN=1, FUN=f))
  colnames(rankings) <- colnames(data)
  rownames(rankings) <- rownames(data)
  return(rankings)
}


#' @title Friedman's test
#'
#' @description This function performs Friedman's test for multiple comparisons
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{parameter}, the two degrees of freedom of the F distribution; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Demsar (2006), page 11
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' 
#' @examples
#' data(data_gh_2008)
#' friedmanTest(data.gh.2008)
#'
friedmanTest <- function (data, ...) {
  N <- dim(data)[1]
  k <- dim(data)[2]
  mr <- colMeans(rankMatrix(data))
  
  friedman.stat <- 12 * N / (k * (k + 1)) * (sum(mr^2) - (k * (k + 1)^2) / 4)
  p.value <- 1 - pchisq(friedman.stat, df=(k - 1))
  
  names(friedman.stat) <- "Friedman's chi-squared"
  parameter <- (k - 1)
  names(parameter) <- c("df")
  method <- "Friedman's rank sum test"
  data.name <- deparse(substitute(data))
  htest.result <- list(statistic=friedman.stat, parameter=parameter, 
                       p.value=p.value, method=method, data.name=data.name)
  class(htest.result) <- "htest"
  return(htest.result)
}


#' @title Friedman's Aligned Ranks test
#'
#' @description This function performs Friedman's Aligned Rank test for multiple comparisons
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{parameter}, the two degrees of freedom of the F distribution; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2008).
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' 
#' @examples
#' data(data_gh_2008)
#' friedmanTest(data.gh.2008)
#'

friedmanAlignedRanksTest <- function (data, ...) {
  N <- dim(data)[1]
  k <- dim(data)[2]
  
  # Compute the p-values for each pair
  dataset.means <- rowMeans(data)
  diff.matrix <- data - matrix(rep(dataset.means, k), ncol=k)
  # To get the rank of a decreasing ordering, change the sign of the differences
  ranks <- matrix(rank(-unlist(diff.matrix)), ncol=k, byrow=FALSE)
  colnames(ranks) <- colnames(data)
  
  r.j <- colSums(ranks)
  r.i <- rowSums(ranks)
  
  r.i.sq.sum <- sum(r.i^2)
  r.j.sq.sum <- sum(r.j^2)
  kn         <- k * N
  
  t.stat.num   <- (k - 1) * (r.j.sq.sum - (kn * N / 4) * (kn + 1)^2)
  t.stat.denom <- ((kn * (kn + 1)) * ((2 * kn) + 1)) / 6 - (1 / k) * r.i.sq.sum
  t.stat       <- t.stat.num / t.stat.denom
  
  p.value <- 1 - pchisq(t.stat, df=(k - 1))
  
  names(t.stat) <- "T"
  parameter <- (k - 1)
  names(parameter) <- c("df")
  method <- "Friedman's Aligned Rank Test for Multiple Comparisons"
  data.name <- deparse(substitute(data))
  htest.result <- list(statistic=t.stat, parameter=parameter, 
                       p.value=p.value, method=method, data.name=data.name)
  class(htest.result) <- "htest"
  return(htest.result)
}


#' @title Quade's test
#'
#' @description This function performs Quade's test for multiple comparisons
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{parameter}, the two degrees of freedom of the F distribution; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2008).
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' 
#' @examples
#' data(data_gh_2008)
#' quadeTest(data.gh.2008)
#' 
quadeTest <- function (data, ...) {
  N <- dim(data)[1]
  k <- dim(data)[2]
  
  # Compute the p-values for each pair
  range.dataset <- apply(data, MARGIN=1, 
                         FUN=function(x) {
                           rg <- range(x)
                           return(diff(rg))
                         })
  q.i  <- rank(range.dataset)
  r.ij <- rankMatrix(data)
  aux.q.i <- matrix(rep(q.i, k), ncol=k, byrow=FALSE)
  s.ij <- aux.q.i * (r.ij - ((k+1)/2))
  w.ij <- aux.q.i * r.ij
  
  s.j <- colSums(s.ij)
  w.j <- colSums(w.ij)
  
  t.j <- w.j / (N * (N + 1) / 2)
  num <- k * (k + 1) * (2 * N + 1) * (k - 1)
  den <- 18 * N * (N + 1)
  nrm <- sqrt(num / den)
  
  a.2 <- N * (N + 1) * (2 * N + 1) * k * (k + 1) * (k - 1) / 72
  b   <- 1 / N * sum(s.j^2)
  t.3 <- ((N - 1) * b) / (a.2 - b)
  
  p.value <- 1 - pf(t.3, df1=(k - 1), df2 = (k - 1) * (N - 1))
  
  names(t.3) <- "T"
  parameter <- (k - 1)
  names(parameter) <- c("df")
  method <- "Quade for Multiple Comparisons"
  data.name <- deparse(substitute(data))
  htest.result <- list(statistic=t.3, parameter=parameter, 
                       p.value=p.value, method=method, data.name=data.name)
  class(htest.result) <- "htest"
  return(htest.result)
}

#' @title Iman Davenport's modification of Friedman's test
#'
#' @description This function performs Iman-Davenport modification of Friedman's test
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{parameter}, the two degrees of freedom of the F distribution; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Demsar (2006), page 11
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' 
#' @examples
#' data(data_gh_2008)
#' imanDavenportTest(data.gh.2008)

imanDavenportTest <- function (data, ...) {
  N <- dim(data)[1]
  k <- dim(data)[2]
  mr <- colMeans(rankMatrix(data))
  
  friedman.stat <- 12 * N / (k * (k + 1)) * (sum(mr^2) - (k * (k + 1)^2) / 4)
  # Iman Davenport correction of Friedman's test
  id.stat <- (N - 1) * friedman.stat / (N * (k - 1) - friedman.stat)
  p.value <- 1 - pf(id.stat, df1=(k - 1), df2=(k - 1) * (N - 1))
  
  names(id.stat)<-"Corrected Friedman's chi-squared"
  parameter <- c((k - 1), (k - 1) * (N - 1))
  names(parameter) <- c("df1", "df2")
  method <- "Iman Davenport's correction of Friedman's rank sum test"
  data.name <- deparse(substitute(data))
  htest.result <- list(statistic=id.stat, parameter=parameter, 
                       p.value=p.value, method=method, data.name=data.name)
  class(htest.result)<-"htest"
  htest.result
}

#' @title Nemenyi test 
#'
#' @description This function performs the Nemenyi test
#' @param data Matrix or data frame where each algorithm is in a column
#' @param alpha Significance level
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{method}, a character string indicating what type of test was performed; \code{data.name}, a character string giving the name of the data and \code{diff.matirx}, a matrix with all the pairwise differences of average rankings
#' @details The test has been implemented according to the version in Demsar (2006), page 7
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @examples
#' data(data_gh_2008)
#' res <- nemenyiTest(data.gh.2008, alpha = 0.1)
#' res
#' res$diff.matrix

nemenyiTest <- function (data, alpha=0.05) {
  k <- dim(data)[2]
  N <- dim(data)[1]
  cd <- getNemenyiCD (alpha=alpha, num.alg=k, num.problems=N)
  
  mean.rank <- colMeans(rankMatrix(data))
  pairs <- generatePairs(k, control=NULL)
  
  differences <- apply(pairs, MARGIN=1, 
                       FUN=function(x) {
                         mean.rank[x[1]] - mean.rank[x[2]]
                       })
  difference.matrix <- matrix(rep(0, k^2), ncol=k)
  difference.matrix[pairs] <- differences
  difference.matrix[pairs[,c(2,1)]] <- differences
  colnames(difference.matrix) <- rownames(difference.matrix)
  colnames(difference.matrix) <- colnames(data)
    
  names(cd)<-"Critical difference"
  parameter <- c(k, (N - 1) * k)
  names(parameter) <- c("k","df")
  method <- "Nemenyi test"
  data.name <- deparse(substitute(data))
  htest.results <- list(statistic=cd, parameter=parameter, method=method, 
                        data.name=data.name, diff.matrix=difference.matrix)
  class(htest.results)<-"htest"
  return(htest.results)
}

#' @title Tukey test 
#'
#' @description This function performs the Tukey test
#' @param data Matrix or data frame where each algorithm is in a column
#' @param alpha Significance level
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{method}, a character string indicating what type of test was performed; \code{data.name}, a character string giving the name of the data and \code{diff.matirx}, a matrix with all the pairwise absolute difference of average values.
#' @details The test has been implemented according to Test 28 in Kanji (2006).
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @examples
#' data(data_gh_2008)
#' res <- tukeyTest(data.gh.2008, alpha=0.1)
#' res
#' res$diff.matrix

tukeyTest <- function (data, alpha=0.05) {
  # Implemented as Test 28 in 100 statistical tests
  k <- dim(data)[2]
  N <- dim(data)[1]
  nu <- k * (N - 1)
  
  # Sample variances
  var.vector <- apply(data, MARGIN=2, FUN=var)
  m.vector <- colMeans(data)
  
  # Total variance. Given that all the samples have the same size, 
  # it is the average variance
  var <- mean(var.vector)
  
  # Studentized range value
  q <- qtukey(alpha, nmeans=k, df=nu)
  W <- q * sqrt(var / N)
  
  # Get the absolute difference of means
  pairs <- generatePairs(k, control=NULL)
  getAbsDiff <- function(x) {
    adf <- abs(m.vector[x[1]] - m.vector[x[2]])
    return(adf)
  }
  differences <- apply(pairs, MARGIN=1, FUN=getAbsDiff)
  difference.matrix <- matrix(rep(NA, k^2), k)
  difference.matrix[pairs] <- differences
  difference.matrix[pairs[,c(2,1)]] <- differences
  colnames(difference.matrix) <- colnames(data)
  colnames(difference.matrix) <- rownames(difference.matrix)
  
  names(W)<-"Critical difference"
  parameter <- c(K=k, df2=(N - 1) * k)
  method <- "Tukey test"
  data.name <- deparse(substitute(data))
  htest.results <- list(statistic=W, parameter=parameter, method=method, 
                        data.name=data.name, diff.matrix=difference.matrix)
  class(htest.results)<-"htest"
  return(htest.results)
}