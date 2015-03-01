#' @title ANOVA test for multiple testing
#'
#' @description This function performs ANOVA test to compare multiple samples
#' @param data Matrix where the test is performed
#' @param ... Ignored
#' @return A list with class "htest" containing the following components: \code{statistic}, the value of the statistic used in the test; \code{p.value}, the p-value for the test; \code{method}, a character string indicating what type of test was performed and \code{data.name}, a character string giving the name of the data.

anova.test <- function (data , ...){
  k <- dim(data)[2]
  N <- dim(data)[1]
  ## Prepare the data
  values <- vector()
  groups <- vector()
  rmv <- sapply(1:k , FUN = function(x){
    values <<- c(values , data[,x])
    groups <<- c(groups , rep(colnames(data)[x],N))
  })
  model <- lm(values ~ groups)
  anova <- anova(model)
  pvalue <- anova$"Pr(>F)"[1]
  Fstatistic <- anova$"F value"[1]
  
  names(Fstatistic)<-"F statistic"
  method <- "ANOVA test for multiple comparisons"
  data.name <- deparse(substitute(data))
  res <- list(statistic = Fstatistic , p.value = pvalue , method = method , data.name = data.name )
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
  order.with.ties <- function(x){
    sorted <- sort(x,decreasing = decreasing)
    sapply(x,FUN=function(i) mean(which(sorted==i)))
  }
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
nemenyi.cd <- function (alpha = 0.05, num.alg , num.problems){
  qa <- qtukey(1-alpha , num.alg , num.alg * (num.problems-1))/sqrt(2)
  qa*sqrt((num.alg*(num.alg+1))/(6*num.problems))
}

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
