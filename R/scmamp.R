#' Statistical comparison of multiple algorithms
#'
#' This package has been develop to simplify the statistical assessment of algorithms when tested in different problems. It includes statistical tests, as well as some plotting functions.
#' @author Borja Calvo, Guzman Santafe
#' @docType package
#' @name scmamp
#' @aliases scmamp-package
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
#' @export
#' @description This function gives access to different alternatives to test the algorithms pair-wise
#' @param results.matrix A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows).
#' @param test Parameter that indicates the statistical test to be used. It can be either a string indicating one of the available test (\code{'t-test'} for paired t-test,  \code{'Wilcoxon'}) for Wilcoxon Signed rank test, \code{'Friedman post'}, for raw p-values in Friedman test's Nemenyi post hoc or a test or a function with, at least, two parameters, \code{x} and \code{y}, which are the two samples to be compared. The function has to return a list that contains, at least, one element called p.value (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param correction A string indicating the type of correction that has to be applied to the p-values. Valid names are \code{Shaffer}, \code{Bergmann-Hommel}, \code{Nemenyi}, \code{Tukey} or any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}.
#' @param ... Special argument used to pass additional parameters to the statistical test or the correction method.
#' @return The function returns a list that contains:
#' \itemize{
#' \item {\code{raw.pvalues} - Matrix with the raw p-values}
#' \item {\code{corrected.pvalues} - Matrix with the corrected p-values}
#' \item {\code{test} - String character indicating the test used}
#' \item {\code{correction} - String character indicating the correction used}
#' }
#' @details The most powerfull method is the dynamic procedure by Bergmann and Hommel, but its computational requirements render this method only applicable for few algorithms. In the current version of the package it accepts up to 9 algorithms. If selected in a dataset with more than 9 columns the correction is automatically changed to \code{'Shaffer'} and a warning is displayed. Shaffer's static approach increases the power without much less cost. Regarding the Nemenyi test, when used with \code{'Friedman post'} test it is equivalent to using Bonferroni's correction implemented in the \code{'mt.adj'} function. Note that Tukey post hoc test is designed for ANOVA and, thus, if selected the \code{test} parameter has to be \code{'t-test'}; otherwise a warning is shown and the \code{test} parameter is changed. Simalrly, the Nemenyi test has to be coupled with Friedman test's post, and if not it is modified an a warning is displayed.
#' @seealso \code{plot.pvalues}, \code{plot.hypothesis}, \code{algorithm.graph}.
#' @examples
#' data(data.garcia.herrera)
#' pwcomp.bh <- pairwise.test(results.matrix = data.garcia.herrera , test = "Friedman post" , correction = "Bergmann Hommel")
#' plot.pvalues(pwcomp.bh$corrected.pvalues)

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
  
  if (correction == "Bergmann Hommel" & k>length(E)){
    warning("Currently the package only supports Bergmann Hommel procedure for 8 or less algorithms. The correction procedure has been changed to 'Shaffer'")
  }
  
  ## Compute the raw p-value matrix
  
  test.name <- "NA"
  if (is.function(test)){
    matrix.raw <- custom.post(results.matrix, test , ...)
    test.name <- paste("Ad hoc function:" , deparse(substitute(test)))
  }else{
    matrix.raw <- switch(test,
                      "t-test" = {
                        test.name <- "Paired t-test"
                        custom.post(results.matrix , function(x,y) t.test(x,y,paired=T)$p.value)
                      },
                      "Wilcoxon" = {
                        test.name <- "Paired Wilcoxon test"
                        custom.post(results.matrix , function(x,y) wilcoxon.signed.test(x,y)$p.value)
                      },
                      "Friedman post" = {
                        test.name <- "Friedman test's post hoc"
                        friedman.post(results.matrix)
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
                        res <- friedman.post(results.matrix)*(k*(k-1)/2)
                        res[res>1] <- 1
                        res
                      },
                      "Tukey" = {
                        correction.name <- "Tukey test"
                        anova.post(results.matrix)
                      },
                      {
                        if (!(correction %in% p.adjust.methods))
                          stop(paste("Non valid method for p.adjust function. Valid options are " , paste(p.adjust.methods,collapse="; "),sep=""))
                        correction.name <- paste("p.adjust functions with method set at '" , correction , "'", sep = "")
                        ## Generate all the pairs to test
                        pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
                        raw.vector <- matrix.raw[pairs]
                        corrected.vector <- p.adjust(p = raw.vector , method = correction )
                        corrected.matrix <- matrix(rep(NA , k^2),ncol=k)
                        corrected.matrix[pairs] <- corrected.vector
                        corrected.matrix[pairs[,c(2,1)]] <- corrected.vector
                        colnames(corrected.matrix) <- rownames(corrected.matrix) <- colnames(raw)
                        corrected.matrix
                      })
  
  list(raw.pvalues = matrix.raw , corrected.pvalues = matrix.adj , test = test.name ,  correction = correction.name)
}






#' @title All vs. best statistical comparisons
#'
#' @export
#' @description Given a results matrix, this function tests the differences between all the algorithms with respect to the best.
#' @param results.matrix A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows).
#' @param test Statistical test to be applied. The first two arguments should be the two vectors to compare. The result has to be a list with an element named p.value. An example of this output is the typical \code{\link{htest}} class returned by R's statistical tests.
#' @param group.by A vector of names or indices of columns that will be used to group the results in \code{results.matrix}.
#' @param alg.col A vector of names or indices of columns indicating which columns contain the results for the algorithms to compare.
#' @param best A string indicating which should be considered the best result. Valid options are \code{'min'} and \code{'max'}.
#' @param summary A function used to summarize the data when looking for the best algorithm. By default, the median is used.
#' @param correction Type of correction to be applied to the p-values. Any method accepted by the \code{\link{p.adjust}} function is a valid option.
#' @param ... Additional parameters to be passed to the test or the summarization function. 
#' @return A list with three matrices, \code{summary}, \code{raw.pvalues} and \code{adj.pvalues}. The first one contains the summarized values, the second one the raw p-values obtained in the comparison and the third one the adjusted pvalues. In the p-values matrices \code{NA} indicates the reference used in that row (i.e., the algorithms with the best value)
#' @examples
#' dir <- system.file("loading_tests",package="scmamp")
#' file <- paste(dir , "beta_complete_comparison.out" , sep="/")
#' data <- read.comparison.file (file = file , alg.cols = c('kakizawa','vitale','boundarykernel','betakernel'))
#' all.vs.best.test (data , group.by = c('size' , 'alpha' , 'beta') , alg.col=4:7 , test = t.test , best='min' , summary = mean , correction = 'hommel' , na.rm = TRUE , paired = TRUE)


all.vs.best.test <- function (results.matrix, test = wilcoxon.signed.test ,  group.by , alg.col , best='max' , summary = mean , correction='holm' ,  ...){
  
  if (length(alg.col)<2) stop("At least two algorithms are required to run the function")
  
  if (is.character(group.by)) group.by <- which(colnames(data) %in% group.by)
  if (is.character(alg.col)) alg.col <- which(colnames(data) %in% alg.col)
  
  ## Remove any index out of bounds
  group.by <- subset(group.by , subset = group.by>0 & group.by<=ncol(results.matrix))
  alg.col <- subset(alg.col , subset = alg.col>0 & alg.col<=ncol(results.matrix))
  
  
  groups <- unique(results.matrix[,group.by])
  if(length(group.by)) groups <- data.frame(groups)
  
  ##########################################################
  test.row <- function (group.values){
    if (length(group.by)==1){
      sub <- results.matrix[,group.by] == group.values  
    }else{
      sub <- apply(results.matrix[,group.by] , MARGIN = 1 , FUN = function(x) all(x ==group.values))
    }
    
    matrix <- subset(results.matrix , 
                     subset = sub, 
                     select = alg.col)
    
    if(nrow(matrix)<15) warning(paste("Running test with less than 15 samples. Combination: " , 
                                      paste(paste(colnames(groups) , "=" , group.values),collapse="; ") , sep = ""))
    
    summ <- apply(matrix , MARGIN = 2 , FUN = function(x) summary (x , ...))
    id <- switch(best,'max' = which.max(summ) , 'min' = which.min(summ) , stop("The 'best' parameter has to be either min or max"))
    res <- sapply(1:ncol(matrix) , FUN = function(i){
      if (i!=id){
        test (matrix[,i] , matrix[,id] , ...)$p.value  
      }else{
        NA
      }
    }) 
    res
  }
  ##########################################################
  
  pvalues <- t(apply(groups , MARGIN = 1 , FUN = test.row))
  pvalues.adj <- matrix(p.adjust(unlist(pvalues)) , byrow = F , ncol = ncol(pvalues))
  colnames(pvalues.adj) <- colnames(results.matrix)[alg.col]
  
  raw.matrix <- cbind(groups,pvalues)
  adj.matrix <- cbind(groups,pvalues.adj)
  colnames(raw.matrix) <- colnames(adj.matrix)
  
  ignore <- which(!colnames(results.matrix) %in% colnames(raw.matrix)) 
  
  summ <- summarize.data(results.matrix , fun = summary , group.by = group.by , ignore = ignore , ...)
  
  list(summary = summ , raw.pvalues = raw.matrix , adj.pvalues = adj.matrix)
}
