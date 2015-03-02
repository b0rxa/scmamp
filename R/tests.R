## Function that returns the order an, in case of ties, the average rank
order.with.ties <- function(x , decreasing = TRUE){
  sorted <- sort(x,decreasing = decreasing)
  sapply(x,FUN=function(i) mean(which(sorted==i)))
}

#' @title Wilcoxon signed-rank est
#'
#' @description This function implements the paired Wilcoxon signed-rank test
#' @param x First sample
#' @param y Second sample
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Demsar (2006), page 7
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.

wilcoxon.signed.test <- function (x , y , ...){
  if (length(x)!=length(y)) stop("This is a paired test, so the two vectors have to have the same length")
  N <- length(x)
  d <- x-y
  o <- order.with.ties(abs(d))
  rp <- sum(o[d>0]) + sum(o[d==0])
  rn <- sum(o[d<0]) + sum(o[d==0])
  t <- min(rp , rn)
  num <- t - 0.25*(N*(N+1))
  den <- sqrt((N*(N+1)*(2*N+1))/24)
  z <- num / den
  
  pvalue <- pnorm(z)  
  names(t)<-"T"
  method <- "Wilcoxon Signed-Rank test"
  data.name <- deparse(substitute(data))
  res <- list(statistic = t , p.value = pvalue , method = method , data.name = data.name )
  class(res)<-"htest"
  res
}

#' @title ANOVA test for multiple comparisons
#'
#' @description This function performs F-test for K populations means
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to Test 22 in Kanji (2006).
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.

anova.test <- function (data , ...){
  
  ## Implemented according to Test 22 in 100 statistical tests
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  ## Variance of the observations with respect to their own sample
  var.1 <- mean(apply(data , MARGIN = 2 , FUN = var))
  
  ## Variance ofthe sample means with respect to the grand mean
  g.m <- mean(unlist(data))
  var.2 <- sum(N*(m.x - g.m)^2) / (k-1)
  
  ## Statistic that follows an F distribution with (k-1); k*(N-1) degrees of freedom
  F = var.2 / var.1 
  
  pvalue <- 1-pf(F,k-1,k*(N-1))
  
  parameters <- c(df1=k-1 , df2=k*(N-1))
  names(F)<-"F statistic"
  method <- "F-test for K population means (Analysis of Variance)"
  data.name <- deparse(substitute(data))
  res <- list(statistic = F , p.value = pvalue , parameters = parameters, method = method , data.name = data.name )
  class(res)<-"htest"
  res
}

#' @title Get the ranking matrix
#'
#' @description This function returns, given a matrix, the ranking of the colums in each row
#' @param data The matrix to rank.
#' @param decreasing Logical value indicating whether the top ranked has to be the highest value or not
#' @return A matrix containing the per-row rankings. In case of ties, the mean rank is obtained (e.g, if there is a tie between the 4th and the 5th column, both are assigned a mean rank of 4.5)
#' @examples
#' data("garcia.herrera")
#' rank.matrix(data.garcia.herrera)

rank.matrix <- function(data , decreasing = TRUE){
  rankings <- t(apply (data , MARGIN = 1 , FUN = order.with.ties))
  colnames(rankings) <- colnames(data)
  rankings
}


#' @title Friedman's test
#'
#' @description This function performs Friedman's test
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{parameter}, the two degrees of freedom of the F distribution; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.
#' @details The test has been implemented according to the version in Demsar (2006), page 11
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' 
#' @examples
#' data("garcia.herrera")
#' iman.davenport.test(data.garcia.herrera)

friedman.test <- function (data , ...){
  N<-dim(data)[1]
  k<-dim(data)[2]
  mr <- colMeans(rank.matrix(data))
  
  friedman.stat <- 12*N/(k*(k+1))*(sum(mr^2) - (k*(k+1)^2)/4)
  p.value <- 1-pchisq(friedman.stat , df = (k-1))
  names(friedman.stat)<-"Friedman's chi-squared"
  parameter <- (k-1)
  names(parameter) <- c("df")
  method <- "Friedman's rank sum test"
  data.name <- deparse(substitute(data))
  res <- list(statistic = friedman.stat , parameter = parameter , p.value = p.value , method = method , data.name = data.name )
  class(res)<-"htest"
  res
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
#' data("garcia.herrera")
#' iman.davenport.test(data.garcia.herrera)

iman.daveport.test <- function (data , ...){
  N<-dim(data)[1]
  k<-dim(data)[2]
  mr <- colMeans(rank.matrix(data))
  
  friedman.stat <- 12*N/(k*(k+1))*(sum(mr^2) - (k*(k+1)^2)/4)
  ## Iman Davenport correction of Friedman's test
  id.stat <- (N-1)*friedman.stat / (N*(k-1) - friedman.stat)
  p.value <- 1-pf(id.stat , df1 = (k-1) , df2 = (k-1)*(N-1))
  names(id.stat)<-"Corrected Friedman's chi-squared"
  parameter <- c((k-1) , (k-1)*(N-1))
  names(parameter) <- c("df1","df2")
  method <- "Iman Davenport's correction of Friedman's rank sum test"
  data.name <- deparse(substitute(data))
  res <- list(statistic = id.stat , parameter = parameter , p.value = p.value , method = method , data.name = data.name )
  class(res)<-"htest"
  res
}


#' @title Nemenyi test critical difference
#'
#' @description This function computes the critical difference for the Nemenyi test
#' @param alpha Significance of the test
#' @param num.alg Number of algorithms (treatments, etc.)
#' @param num.problems Number of problems (samples)
#' @return Critical difference of averge rankings. When the difference is, in absolute value, greater than this value the null hypothesis (no differences) is rejected.
#' @details The test has been implemented according to the version in Demsar (2006), page 11
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @seealso \code{\link{nemenyi.test}}
#' 
nemenyi.cd <- function (alpha = 0.05, num.alg , num.problems){
  qa <- qtukey(1-alpha , num.alg , num.alg * (num.problems-1))/sqrt(2)
  qa*sqrt((num.alg*(num.alg+1))/(6*num.problems))
}

#' @title Nemenyi test 
#'
#' @description This function performs the Nemenyi test
#' @param data Matrix or data frame where each algorithm is in a column
#' @param alpha Significance level
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{method}, a character string indicating what type of test was performed; \code{data.name}, a character string giving the name of the data and \code{diff.matirx}, a matrix with all the pairwise differences of average rankings
#' @details The test has been implemented according to the version in Demsar (2006), page 7
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
nemenyi.test <- function (data , alpha = 0.05){
  k <- dim(data)[2]
  N <- dim(data)[1]
  cd <- nemenyi.cd (alpha = alpha , num.alg = k , num.problems = N)
  
  mean.rank <- colMeans(rank.matrix(data))
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  differences <- apply(pairs , MARGIN = 1 , FUN = function(x){mean.rank[x[1]] - mean.rank[x[2]]})
  difference.matrix <- matrix(rep(0 , k^2) , ncol = k)
  difference.matrix[pairs] <- differences
  difference.matrix[pairs[,c(2,1)]] <- differences
  colnames(difference.matrix) <- rownames(difference.matrix) <- colnames(data)
  
  
  names(cd)<-"Critical difference"
  parameter <- c(k , (N-1)*k)
  names(parameter) <- c("k","df")
  method <- "Nemenyi test"
  data.name <- deparse(substitute(data))
  res <- list(statistic = cd , parameter = parameter , method = method , data.name = data.name , diff.matrix = difference.matrix )
  class(res)<-"htest"
  res
}

#' @title Tukey test 
#'
#' @description This function performs the Tukey test
#' @param data Matrix or data frame where each algorithm is in a column
#' @param alpha Significance level
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{method}, a character string indicating what type of test was performed; \code{data.name}, a character string giving the name of the data and \code{diff.matirx}, a matrix with all the pairwise absolute difference of average values.
#' @details The test has been implemented according to Test 28 in Kanji (2006).
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.

tukey.test <- function (data , alpha = 0.05){
  ## Implemented as Test 28 in 100 statistical tests
  
  k <- dim(data)[2]
  N <- dim(data)[1]
  nu <- k*(N-1)
  
  ## Sample variances
  var.vector <- apply(data , MARGIN = 2 , FUN = var)
  m.vector <- colMeans(data)
  
  ## Total variance. Given that all the samples have the same size, it is the average variance
  var <- mean(var.vector)
  
  ## Studentized range value
  q <- qtukey(alpha , nmeans = k , df = nu)
  W <- q * sqrt(var / N)
  
  ## Get the absolute difference of means
  f<-function(x) abs(m.vector[x[1]] - m.vector[x[2]])
  differences <- apply(pairs , MARGIN = 1 , FUN = f)
  difference.matrix <- matrix(rep(NA , k^2) , k)
  difference.matrix[pairs] <- differences
  difference.matrix[pairs[,c(2,1)]] <- differences
  colnames(difference.matrix) <- rownames(difference.matrix) <- colnames(data)
  
  names(W)<-"Critical difference"
  parameter <- c(K = k , df2 = (N-1)*k)
  method <- "Tukey test"
  data.name <- deparse(substitute(data))
  res <- list(statistic = W , parameter = parameter , method = method , data.name = data.name , diff.matrix = difference.matrix )
  class(res)<-"htest"
  res
}
