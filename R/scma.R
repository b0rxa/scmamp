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
#' @param test Parameter that indicates the statistical test to be used. It can be either a string indicating one of the available test (\code{'t-test'} for paired t-test,  \code{'Wilcoxon'}) for Wilcoxon Signed rank test, \code{'Friedman post'}, for raw p-values in Friedman test's Nemenyi post hoc or a function that, given two parameters, \code{x} and \code{y}, return the p-value associated with the comparison. For compatibility, this function has to have the ... special argument.
#' @param correction A string indicating the type of correction that has to be applied to the p-values. Valid names are \code{Shaffer}, \code{Bergmann-Hommel}, \code{Nemenyi}, \code{Tukey} or any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}.
#' @param ... Special argument used to pass additional parameters to the statistical test or the correction method.
#' @return The function returns a list that contains:
#' \itemize{
#' \item {\code{raw.pvalues} - Matrix with the raw p-values}
#' \item {\code{corrected.pvalues} - Matrix with the corrected p-values}
#' \item {\code{test} - String character indicating the test used}
#' \item {\code{correction} - String character indicating the correction used}
#' }
#' @details The most powerfull method is the dynamic procedure by Bergmann and Hommel, but its computational requirements render this method only applicable for few algorithms. In the current version of the package it accepts up to 8 classifiers. If selected in a dataset with more than 8 columns the correction is automatically changed to \code{'Shaffer'} and a warning is displayed. Shaffer's static approach increases the power without much less cost. Regarding the Nemenyi test, when used with \code{'Friedman post'} test it is equivalent to using Bonferroni's correction implemented in the \code{'mt.adj'} function. Note that Tukey post hoc test is designed for ANOVA and, thus, if selected the \code{test} parameter has to be \code{'t-test'}; otherwise a warning is shown and the \code{test} parameter is changed. Simalrly, the Nemenyi test has to be coupled with Friedman test's post, and if not it is modified an a warning is displayed.
#' @seealso \code{plot.pvalues}, \code{plot.hypothesis}, \code{algorithm.graph}.

pairwise.test <- function(results.matrix ,  test="Friedman post" , correction="Shaffer" , ...){
  k <- dim(results.matrix)[2]
  N <- dim(results.matrix)[1]

  ## Corrections of the elections  
  if (correction == "Nemenyi" && (is.function(test) || (!is.function(test) & test!="Friedman post"))){
    warning("The Nemenyi test is Friedman test's post hoc and, thus, can only be coupled with this option. Ignoring the 'test' parameter and setting it to 'Friedman post'")
    test="Friedman post"
  }
  
  if (correction == "Tukey" && (is.function(test) || (!is.function(test) & test!="t-test"))){
    warning("The Tukey test is, essentially, a corrected version of the pairwise t-test, so it has to be coupled with this option. Ignoring the 'test' parameter and setting it to 't-test'")
    test="t-test"
  }
  
  if (correction == "Bergmann Hommel" & k>8){
    warning("Currently the package only supports Bergmann Hommel procedure for 8 or less algorithms. The correction procedure has been changed to 'Shaffer'")
  }
  
  ## Compute the raw p-value matrix
  
  test.name <- "NA"
  if (is.function(test)){
    matrix.raw <- custom.post(data, test)
    test.name <- paste("Ad hoc function:" , deparse(substitute(test)))
  }else{
    matrix.raw <- switch(test,
                      "t-test" = {
                        test.name <- "Paired t-test"
                        custom.post(data , function(x,y) t.test(x,y,paired=T)$p.value)
                      },
                      "Wilcoxon" = {
                        test.name <- "Paired Wilcoxon test"
                        custom.post(data , function(x,y) wilcoxon.signed.test(x,y)$p.value)
                      },
                      "Friedman post" = {
                        test.name <- "Friedman test's post hoc"
                        friedman.post(data)
                      },
                      stop("Unknown test. Valid options in the current version are 'Friedman post', 'Wilcoxon' ,'ANOVA post' and 't-test'. Alternatively, you can pass a function that performs a paired statistical texts which should have, at least, two parameters, 'x' and 'y' and returns the p-value associted to the comparison"))
  }
  
  matrix.adj <- switch(correction ,
                      "Shaffer" = {
                        correction.name <- "Shaffer static"
                        shaffer.static(matrix.raw)
                      },
                      "Bergmann Hommel" = {
                        correction.name <- "Bermannn Hommel dynamic"
                        bergmann.hommel.dynamic(matrix.raw)
                      },
                      "Nemenyi" = {
                        correction.name <- "Nemenyi test"
                        res <- friedman.post(data)*(k*(k-1)/2)
                        res[res>1] <- 1
                        res
                      },
                      "Tukey" = {
                        correction.name <- "Tukey test"
                        anova.post(data)
                      },
                      {
                        if (!(correction %in% p.adjust.methods))
                          stop(paste("Non valid method for p.adjust function. Valid options are " , paste(p.adjust.methods,collapse=";"),sep=""))
                        correction.name <- paste("p.adjust functions with method set at '" , correction , "'", sep = "")
                        raw.vector <- matrix.raw[pairs]
                        corrected.vector <- p.adjust(p = raw.vector , method = correction )
                        corrected.matrix <- matrix(rep(NA , dim(raw)[1]^2),ncol=dim(raw)[1])
                        corrected.matrix[pairs] <- corrected.vector
                        corrected.matrix[pairs[,c(2,1)]] <- corrected.vector
                        colnames(corrected.matrix) <- rownames(corrected.matrix) <- colnames(raw)
                        corrected.matrix
                      })
  
  list(raw.pvalues = matrix.raw , corrected.pvalues = matrix.adj , test = test.name ,  correction = correction.name)
}