#' Statistical comparison of multiple algorithms
#'
#' This package has been develop to simplify the statistical assessment of algorithms when tested in different problems. It includes statistical tests, as well as some plotting functions.
#' @author Borja Calvo
#' @docType package
#' @name scma
#' @aliases scma-package
NULL


#' Example in Garcia and Herrera (2008)
#'
#' Dataset corresponding to the accuracy of 5 classifiers in 30 datasets
#' Each algorithm is in a column
#'
#' @format A data frame with 5 columns and 30 rows
#' @source García S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for all Pairwise Comparisons. \emph{Journal of Machine Learning Research}. 9, 2677-2694.
#' @name data.garcia.herrera
NULL


#' @title Statistical pairwise comparisons of algorithms
#'
#' @description This function gives access to different alternatives to test the algorithms pair-wise
#' @param results.matrix A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows).
#' @param test Parameter that indicates the statistical test to be used. It can be either a string indicating one of the available test (\code{t-test} and \code{Wilcoxon}) or a function that, given two parameters, \code{x} and \code{y}, return the p-value associated with the comparison. For compatibility, this function has to have the ... special argument.
#' @param correction A string indicating the type of correction that has to be applied to the p-values. Valid names are \code{Shaffer}, \code{Bergmann-Hommel}, \code{Nemenyi}, \code{Tukey} or any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}
#' @param ... Special argument used to pass additional parameters to the statistical test or the correction method.
#' @return The function returns a list that contains:
#' \itemize{
#' \item {\code{raw.pvalues} - Matrix with the raw p-values}
#' \item {\code{corrected.pvalues} - Matrix with the corrected p-values}
#' \item {\code{test} - String character indicating the test used}
#' \item {\code{correction} - String character indicating the correction used}
#' }
#' @details The most powerfull method is the dynamic procedure by Bergmann and Hommel, but its computational requirements render this method only applicable for few algorithms. In the current version of the package it accepts up to 8 classifiers. Shaffer's static approach increases the power without much less cost. Regarding the Nemenyi test, it is equivalent to using Bonferroni's correction implemented in the \code{'mt.raw2adjp'} function. Note that Tukey post hoc test is designed for ANOVA and, thus, if selected the \code{test} parameter is ignores. Simalrly, the Nemenyi test is the non parametric version and thus, it also ignores the \code{test} parameter.
#' @seealso \code{plot.pvalues}, \code{plot.hypothesis}, \code{algorithm.graph} , \code{scma}

pairwise.test <- function(results.matrix ,  test="Wilcoxon" , correction="Shaffer" , ...){
  k <- dim(results.matrix)[2]
  N <- dim(results.matrix)[1]
  
  if (correction == "Nemenyi" && (is.function(test) || (!is.function(test) & test!="Wilcoxon"))){
    warning("Nemenyi test has to be coupled with Wilcoxon test. Ignoring the 'test' parameter and setting it to 'Wilcoxon'")
    test="Wilcoxon"
  }
  
  ## Define the function to perform the post-hoc test
  test.name <- "NA"
  if (is.function(test)){
    pw.test <- test
    test.name <- paste("Ad hoc function:" , deparse(substitute(test)))
    matrix.h <- NULL
    matrix.raw <- NULL
    matrix.adj <- NULL
  }else{
    pw.test <- switch(test,
                      "t-test" = {
                        test.name <- "Paired t-test"
                        function(x,y,...) t.test(x,y, paired = T , ...)$p.value
                      },
                      "Wilcoxon" = {
                        test.name <- "Paired Wilcoxon test"
                        function(x,y,...) wilcox.test(x,y, paired = T , ...)$p.value
                      },
                      stop("Unknown test. Valid options in the current version are 't-test' and 'Wilcoxon'. Alternatively, you can pass a function that performs a paired statistical texts which should have, at least, two parameters, 'x' and 'y' and returns the p-value associted to the comparison"))
  }
  
  ## Generate all the pairs to test
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  ## Run the tests
  if (correction=="Tukey"){
    warning("Tukey is the post hoc ANOVA. Ignoring the 'test' parameter")
    matrix.raw <- NA
    test.name <- "NA"
    correction.name <- "Tukey post hoc"
    matrix.corrected <- tukey.test(results.matrix)
  }else{  
    f<-function(x) pw.test(x = results.matrix[ , x[1]] , y = results.matrix[ , x[2]] , ...)
    pvalues <- apply(pairs , MARGIN = 1 , FUN = f)
    matrix.raw <- matrix(rep(NA , k^2) , k)
    matrix.raw[pairs] <- pvalues
    matrix.raw[pairs[,c(2,1)]] <- pvalues
    colnames(matrix.raw) <- colnames(results.matrix)
    rownames(matrix.raw) <- colnames(results.matrix)  
    ## Do the correction
    correction.name <- "NA"
    correct <- switch(correction ,
                      "Shaffer" = {
                        correction.name <- "Shaffer static"
                        function(raw , ...) shaffer.static(raw)
                      },
                      "Bergmann Hommel" = {
                        if (k>8) stop("Bergmann Hommel dynamic approach can only be used with 8 or less algorithms. Try instead 'Shaffer' for the static approach")
                        correction.name <- "Bermannn Hommel dynamic"
                        function(raw , ...) bergmann.hommel.dynamic(raw)
                      },
                      "Nemenyi" = {
                        correction.name <- "Nemenyi test"
                        function(raw , ...) raw*dim(raw)[2]
                      },
                      {
                        if (!(correction %in% p.adjust.methods))
                          stop(paste("Non valid method for p.adjust function. Valid options are " , paste(p.adjust.methods,collapse=";"),sep=""))
                        correction.name <- paste("Multtest's implementation of " , correction , sep = "")
                        function(raw , ...){
                          raw.vector <- raw[pairs]
                          corrected.vector <- p.adjust(p = raw.vector , method = correction)
                          corrected.matrix <- matrix(rep(NA , dim(raw)[1]^2),ncol=dim(raw)[1])
                          corrected.matrix[pairs] <- corrected.vector
                          corrected.matrix[pairs[,c(2,1)]] <- corrected.vector
                          colnames(corrected.matrix) <- rownames(corrected.matrix) <- colnames(raw)
                          corrected.matrix
                        }
                        })
    matrix.corrected <- correct(matrix.raw)
  }
  list(raw.pvalues = matrix.raw , corrected.pvalues = matrix.corrected , test = test.name ,  correction = correction.name)
}