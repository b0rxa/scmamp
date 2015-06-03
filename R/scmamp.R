#' Statistical comparison of multiple algorithms
#'
#' This package has been develop to simplify the statistical assessment of algorithms when tested in different problems. It includes statistical tests, as well as some plotting functions.
#' @author Borja Calvo, Guzman Santafe
#' @docType package
#' @name scmamp
#' @aliases scmamp-package
#' @details
NULL


#' Example in Garcia and Herrera (2008)
#'
#' Dataset corresponding to the accuracy of 5 classifiers in 30 datasets
#' Each algorithm is in a column
#'
#' @format A data frame with 5 columns and 30 rows
#' @source García S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for all Pairwise Comparisons. \emph{Journal of Machine Learning Research}. 9, 2677-2694.
#' @name data.gh.2008
#' @details
NULL

#' @title Auxiliar function to perform the post hoc tests
#'
#' @param data A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows). It can contain additional columns, but in that case it is mandatory to indicate, in the next parameter, which columns contain the algorithm information, unless these columns are (all) used in the \code{group.by} parameter.
#' @param test Parameter that indicates the statistical test to be used. It can be a string indicating one of the available test (\code{'friedman'} for Friedman test, \code{'aligned ranks'}) for Friedman aligned ranks test, \code{'quade'}, for Quade test, \code{anova}, for ANOVA test. Alternatively, it can be a function that recives as first argument a matrix containing the columns to be compared, and that returns a a list with, at least, an element named \code{p.value} (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param control Colum used as control. If NULL, all the pairwise comparisons are performed
#' @details
#' @return The p-values obtained in the comparisons.
#' 
runPostHoc <- function (data, test, control, ...) {
  # Auxiliar function to conduct the post hoc test
  # Args:
  #   data:    Dataset where the test is conducted. It should contain only the values 
  #            to compare
  #   test:    Test to be performed
  #   control: Algorithm used as control
  #
  # Returns:
  #   The obtained p-values
  #
  if (is.function(test)) {
    matrix.raw <- customPost(data, test, ...)
  }else{
    matrix.raw <- switch(test,
                         "t-test"= {
                           test.name <- "Paired t-test"
                           customPost(data=data, control=control, 
                                       test=function(x, y) {
                                         return(t.test(x, y, paired=TRUE))
                                         })
                       },
                       "wilcoxon"= {
                         test.name <- "Paired Wilcoxon test"
                         customPost(data=data, control=control, 
                                    test=function(x,y) {
                                      return(wilcoxonSignedTest(x,y))
                                    })
                       },
                       "friedman"= {
                         test.name <- "Friedman test's post hoc"
                         friedmanPost(data=data, 
                                      control=control)
                       },
                       "aligned ranks"= {
                         test.name <- "Friedman Aligned Rank test's post hoc"
                         friedmanAlignedRanksPost(data=data, 
                                                  control=control)
                       },
                       "quade"= {
                         test.name <- "Quade test's post hoc"
                         quadePost(data=data, 
                                   control=control)
                       },
                       stop("Unknown test. Valid options in the current version ",
                            "are 'Wilcoxon', 'Friedman', 'Aligned ranks', 'Quade', ",
                            "and 't-test'. Alternatively, you can pass a function ",
                            "that performs a paired statistical texts which ",
                            "should have, at least, two parameters, 'x' and 'y' ",
                            "and returns the p-value associted to the comparison"))
  }

  return(matrix.raw)
}


#' @title Post hoc test for multiple comparison analises
#'
#' @export
#' @description This function is a wrapper to run the post hoc tests.
#' @param data A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows). It can contain additional columns, but in that case it is mandatory to indicate, in the next parameter, which columns contain the algorithm information, unless these columns are (all) used in the \code{group.by} parameter.
#' @param algorithms Vector with either the names or the indices of the columns that contain the values to be tested. If not provided, the function assumes that all the columns except those indicated in \code{group.by} represent the results obtained by an algorithm.
#' @param group.by Vector with either the names or the indices of the columns to be used to group the data. A test is run for each group.
#' @param test Parameter that indicates the statistical test to be used. It can be a string indicating one of the available test (\code{'friedman'} for Friedman test, \code{'aligned ranks'}) for Friedman aligned ranks test, \code{'quade'}, for Quade test, \code{anova}, for ANOVA test. Alternatively, it can be a function that recives as first argument a matrix containing the columns to be compared, and that returns a a list with, at least, an element named \code{p.value} (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param correct A string indicating the type of correction that has to be applied to the p-values, in case more than one group is tested. Valid names are \code{holland}, for Holland's method, \code{finner}, for Finner's method, \code{Rom}, for Rom's method, \code{li}, for Li's method, or any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}. Alternatively, a function that performs the correction can be passed. This function has to recieve, as first argument, a vector of pvalues to be corrected and has to return a verctor with the corrected p-values {\emph in the same order} as the input vector.
#' @param ... Special argument used to pass additional parameters to the statistical test or the correction method.
#' @return In case the \code{group.by} argument is not provided (or it is \code{NULL}), the function return an object of class \code{htest}. If columns for grouping are provided, then the function returns a matrix that includes, for each group, the values of the \code{group.by} columns, the raw p-value and the corrected p-value.
#' @details 
#' @examples
#' # Grouped data
#' data(data_blum_2015)
#' multipleComparisonTest (data=data.blum.2015, 
#'                         algorithms=c("FrogCOL", "FrogMIS", "FruitFly"), 
#'                         group.by=c("Size", "Radius"), 
#'                         test="quade", correct="li")
#' # Not grouped data
#' data(data_gh_2008)
#' multipleComparisonTest (data=data.gh.2008, test="aligned ranks", 
#'                         correct="hochberg")
#' 
postHocTest <- function (data, algorithms=NULL, group.by=NULL, test="friedman", 
                         control=NULL, use.rank=FALSE, sum.fun=mean, 
                         correct="finner", alpha=0.05, ... ) {
 
  # Convert string columns to their corresponding ID
  if (!is.null(group.by) & is.character(group.by)) {
    if (!all(group.by %in% colnames(data))) {
      warning("Not all the columns indicated in the 'group.by' argument are in ",
              "the dataset. None existing columns will be ignored")
    }
    group.by <- which(colnames(data) %in% group.by)
  }

  # Remove any index out of bounds
  sbt <- group.by > 0 & group.by <= ncol(data)
  if (!all(sbt)) {
    warning("Not all the columns indicated in the 'group.by' argument are in ",
            "the dataset. Out of range values will be ignored.")
  }
  group.by <- subset(group.by, subset=sbt)
 
  # In case there is not a list of algorithms, then all the columns except those
  # used to group the data are regarded as algorithm results
  if (is.null(algorithms)) {
    algorithms <- which(!(1:ncol(data) %in% group.by))
  } else {
    if (is.character(algorithms)) {
      if (!all(algorithms %in% colnames(data))) {
        warning("Not all the columns indicated in the 'algorithms' argument are in ",
                "the dataset. None existing columns will be ignored")
      }
      algorithms <- which(colnames(data) %in% algorithms)
    }
  
    sbt <- algorithms > 0 & algorithms <= ncol(data)
    if (!all(sbt)) {
      warning("Not all the columns indicated in the 'group.by' argument are in ",
              "the dataset. Out of range values will be ignored.")
    }
    algorithms <- subset(algorithms, subset=sbt)
  }
  
  # Just in case ...
  if (length(algorithms) < 2) {
    stop("At least two algorithms are required to run the function")
  }
  
  # Use name for the control to avoid problems when filtering
  if (!is.null(control)) {
    if (is.character(correct) & (correct == "shaffer" | correct == "bergmannHommel")) {
      warning("Shaffer's and Bergman and Hommel's correction can only be used ",
              "when all the pairs are compared. For comparisons with a control ",
              "use any of the other corrections that do not take into account ",
              "the particular nature of all pairwise comparisons.")
    }
    if (is.numeric(control)) {
      control <- names(data)[control]
    }
  }  
  
  # Prepare the correction
  if(is.character(correct)) {
    correct.name <- correct
    correct <- switch (correct,
                       "shaffer"=adjustShaffer,
                       "bergmann"=adjustBergmannHommel,
                       "holland"=adjustHolland,
                       "finner"=adjustFinner,
                       "rom"={
                         function(pvalues) {
                           return(adjustRom(pvalues=pvalues, alpha=alpha))
                         }
                       },
                       "li"=adjustLi,
                       {
                         if (!(correct %in% p.adjust.methods)){
                           stop("Non valid method for p.adjust function. ",
                                "Valid options are ",
                                paste(p.adjust.methods, collapse="; "),
                                ". Additionally, 'holland', 'finner', 'rom' and ",
                                "'li' are also valid options")
                         }
                         function(pvalues, ...) {
                           p.adjust(p=pvalues, method=correct.name)
                         }
                       })
  } else {
    correct.name <- deparse(substitute(test))  
  }

  # Build the summary matrix
  if (is.null(group.by)){
   if (use.rank) {
     aux <- rankMatrix(data=data[, algorithms], ...)
   } else {
     aux <- data[, algorithms]
   }
   # Note that in aux the group.by columns are at the begining
   # Some operations may change the name of the algorithms (special characters)
   sum.matrix <- summarizeData(data=aux, fun=sum.fun, group.by=NULL, ...)
   names(sum.matrix) <- colnames(data)[algorithms]
  } else {
    if (use.rank) {
      aux <- cbind(data[, group.by], rankMatrix(data=data[, algorithms]), ...)
    } else {
      aux <- cbind(data[, group.by], data[, algorithms])
    }
    # Note that in aux the group.by columns are at the begining
    sum.matrix <- summarizeData(data=aux, fun=sum.fun, group.by=1:length(group.by), ...)
    # Some operations may change the name of the algorithms (special characters)
    colnames(sum.matrix) < names(data)[c(group.by, algorithms)]
  }
  
  if (!is.null(group.by)){
    # Generate all the groups as a data.frame
    groups <- unique(data[, group.by])
    # In case the result is a vector, convert it into a data.frame
    if(length(group.by)==1) {
      groups <- data.frame(groups)
    }
    names(groups) <- names(data)[group.by]
        
    getRawPvalues <- function (i) {
      # Filter the data
      rows <- rep(TRUE, nrow(data))
      for (j in seq(along.with=group.by)) {
        g <- group.by[j]
        rows <- rows & data[, g]==groups[i, j]
      }
      data.sub <- subset (data, rows)[, algorithms]
      # Check the control algorithm
      if (use.rank){
        aux <- rankMatrix(data.sub)
      } else {
        aux <- data.sub
      }
      ref <- apply(aux, MARGIN=2, FUN=sum.fun)
      if (!is.null(control) & is.character(control)) {
        if (control=="max") {
          control <- which.max(ref)
        } else if (control=="min") {
          control <- which.min(ref)
        }
      }
      group.result <- runPostHoc (data.sub, test=test, control=control)
      return(group.result)
    }
    
    group.raw.pval <- unlist(lapply (1:nrow(groups), FUN=getRawPvalues))
    group.corrected.pval <- correct(group.raw.pval)
    
    # Now we create either the arrays or the matrix, depending on whether we have
    # a control or not
    k <- length(algorithms)
    p <- nrow(groups)
    if (is.null(control)){
      dim(group.raw.pval) <- c(k, k, p)
      dim(group.corrected.pval) <- c(k, k, p)
      
      # Name the dimensions
      group.names <- paste(names(groups)[1], groups[, 1], sep=": ")
      if(ncol(groups) > 1) {
        for (j in 2:ncol(groups)) {
          group.names <- cbind(group.names, 
                               paste(names(groups)[j], groups[, j], sep=": "))
        }
        
        group.names <- apply(group.names, MARGIN=1, 
                             FUN=function(x) {
                               return(paste(x, collapse="; "))
                             })
      }
      
      dimnames(group.raw.pval) <- list(names(data)[algorithms],
                                       names(data)[algorithms], 
                                       group.names)
      
      dimnames(group.corrected.pval) <- list(names(data)[algorithms],
                                             names(data)[algorithms],
                                             group.names)
    } else {
      dim(group.raw.pval) <- c(k,p)
      dim(group.corrected.pval) <- c(k,p)
      group.raw.pval <- t(group.raw.pval)
      group.corrected.pval <- t(group.corrected.pval)
      colnames(group.raw.pval) <- names(data)[algorithms]
      colnames(group.corrected.pval) <- names(data)[algorithms]
      group.raw.pval <- cbind(groups, group.raw.pval)
      group.corrected.pval <- cbind(groups, group.corrected.pval)
    }
    raw.pval <- group.raw.pval
    corrected.pval <- group.corrected.pval
  } else {
    if (is.null(control)){
      # Check the control algorithm
      aux <- data[, algorithms]
      if (use.rank){
        aux <- rankMatrix(aux, ...)
      }
      ref <- apply(aux, MARGIN=2, FUN=sum.fun)
      if (is.character(control)) {
        if (control=="max") {
          control <- which.max(ref)
        } else if (control=="min") {
          control <- which.min(ref)
        }
      }
    }
    raw.pval <- runPostHoc(data[, algorithms], test=test, control=control, ...)
    corrected.pval <- correct(raw.pval)
  }
  
  results <- list(summary=sum.matrix, raw.pval=raw.pval, 
                  corrected.pval=corrected.pval)
  return(results)
}



#' @title Tests for multiple comparisons
#'
#' @export
#' @description This function is a wrapper to multiple comparison tests.
#' @param data A matrix or data.frame containing the results obtained by the algorithms (columns) in each problem (rows). It can contain additional columns, but in that case it is mandatory to indicate, in the next parameter, which columns contain the algorithm information, unless these columns are (all) used in the \code{group.by} parameter.
#' @param algorithms Vector with either the names or the indices of the columns that contain the values to be tested. If not provided, the function assumes that all the columns except those indicated in \code{group.by} represent the results obtained by an algorithm.
#' @param group.by Vector with either the names or the indices of the columns to be used to group the data. A test is run for each group.
#' @param test Parameter that indicates the statistical test to be used. It can be a string indicating one of the available test (\code{'friedman'} for Friedman test, \code{'aligned ranks'}) for Friedman aligned ranks test, \code{'quade'}, for Quade test, \code{anova}, for ANOVA test. Alternatively, it can be a function that recives as first argument a matrix containing the columns to be compared, and that returns a a list with, at least, an element named \code{p.value} (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param correct A string indicating the type of correction that has to be applied to the p-values, in case more than one group is tested. Valid names are \code{holland}, for Holland's method, \code{finner}, for Finner's method, \code{Rom}, for Rom's method, \code{li}, for Li's method, or any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}. Alternatively, a function that performs the correction can be passed. This function has to recieve, as first argument, a vector of pvalues to be corrected and has to return a verctor with the corrected p-values {\emph in the same order} as the input vector.
#' @param ... Special argument used to pass additional parameters to the statistical test or the correction method.
#' @return In case the \code{group.by} argument is not provided (or it is \code{NULL}), the function return an object of class \code{htest}. If columns for grouping are provided, then the function returns a matrix that includes, for each group, the values of the \code{group.by} columns, the raw p-value and the corrected p-value.
#' @details 
#' @examples
#' # Grouped data
#' data(data_blum_2015)
#' multipleComparisonTest (data=data.blum.2015, 
#'                         algorithms=c("FrogCOL", "FrogMIS", "FruitFly"), 
#'                         group.by=c("Size", "Radius"), 
#'                         test="quade", correct="li")
#' # Not grouped data
#' data(data_gh_2008)
#' multipleComparisonTest (data=data.gh.2008, test="aligned ranks", 
#'                         correct="hochberg")
#' 
multipleComparisonTest <- function (data, algorithms=NULL, group.by=NULL, 
                                    test="aligned ranks", correct="finner", 
                                    alpha=0.05, ...){

  # Convert string columns to their corresponding ID
  if (!is.null(group.by) & is.character(group.by)) {
    if (!all(group.by %in% colnames(data))) {
      warning("Not all the columns indicated in the 'group.by' argument are in ",
              "the dataset. None existing columns will be ignored")
    }
    group.by <- which(colnames(data) %in% group.by)
  }
  
  # Remove any index out of bounds
  sbt <- group.by > 0 & group.by <= ncol(data)
  if (!all(sbt)) {
    warning("Not all the columns indicated in the 'group.by' argument are in ",
            "the dataset. Out of range values will be ignored.")
  }
  group.by <- subset(group.by, subset=sbt)
  
  # In case there is not a list of algorithms, then all the columns except those
  # used to group the data are regarded as algorithm results
  if (is.null(algorithms)) {
    algorithms <- which(!(1:ncol(data) %in% group.by))
  } else {
    # Same processing as with 'group.by'
    if (is.character(algorithms)) {
      if (!all(algorithms %in% colnames(data))) {
        warning("Not all the columns indicated in the 'algorithms' argument are in ",
                "the dataset. None existing columns will be ignored")
      }
      algorithms <- which(colnames(data) %in% algorithms)
    }
  
    sbt <- algorithms > 0 & algorithms <= ncol(data)
    if (!all(sbt)) {
      warning("Not all the columns indicated in the 'group.by' argument are in ",
              "the dataset. Out of range values will be ignored.")
    }
    algorithms <- subset(algorithms, subset=sbt)
  }
  
  # Just in case ...
  if (length(algorithms) < 2) {
    stop("At least two algorithms are required to run the function")
  }

  # Prepare the test function 
  
  if(is.character(test)) {
    test.name <- test
    test <- switch (test,
                    "friedman"=friedmanTest,
                    "iman"=imanDavenportTest,
                    "aligned ranks"=friedmanAlignedRanksTest,
                    "quade"=quadeTest,
                    "anova"=anovaTest,
                    stop("Unknown test. Valid options are 'friedman', ",
                         "'aligned ranks', 'quade', 'anova' or a function that ",
                         "gets as input a data.frame or matrix where each ",
                         "algorithm is in a column"))
  } else {
    test.name <- deparse(substitute(test))  
  }
  
  # Prepare the correction function
  if(is.character(correct)) {
    correct.name <- correct
    correct <- switch (correct,
                       "holland"=adjustHolland,
                       "finner"=adjustFinner,,
                       "rom"={
                         function(pvalues) {
                           return(adjustRom(pvalues=pvalues, alpha=alpha))
                         }
                       },
                       "li"=adjustLi,
                       {
                         if (!(correct %in% p.adjust.methods)){
                           stop("Non valid method for p.adjust function. ",
                                "Valid options are ",
                                paste(p.adjust.methods, collapse="; "),
                                ". Additionally, 'holland', 'finner', 'rom' and ",
                                "'li' are also valid options")
                         }
                         function(pvalues, ...) {
                           p.adjust(p=pvalues, method=correct.name)
                         }
                       })
  } else {
    correct.name <- deparse(substitute(test))  
  }

  if (!is.null(group.by)){
    # Generate all the groups as a data.frame
    groups <- unique(data[, group.by])
    # In case the result is a vector, convert it into a data.frame
    if(length(group.by)==1) {
      groups <- data.frame(groups)
    }
    names(groups) <- names(data)[group.by]
  
    getRawPvalues <- function (i) {
      # Filter the data
      rows <- rep(TRUE, nrow(data))
      for (j in seq(along.with=group.by)) {
        g <- group.by[j]
        rows <- rows & data[, g]==groups[i, j]
      }
      data.sub <- subset (data, subset=rows, select=algorithms)
      group.result <- test(data.sub, ...)$p.value
      return(group.result)
    }
    group.raw.pval <- sapply (1:nrow(groups), FUN=getRawPvalues)
  
    # Correct the p-values and generate the final matrix
    group.corrected.pval <- correct(group.raw.pval, ...)
    res.matrix <- cbind(groups, group.raw.pval, group.corrected.pval)
    colnames(res.matrix) <- c(colnames(groups), 
                              paste("Raw_p-val_", test.name, sep=""), 
                              paste("Corrected_p-val_", correct.name, sep=""))
    result <- res.matrix
  } else {
    data.multipleComparisonTest <- data[, algorithms]
    result <- test(data.multipleComparisonTest, ...)
  }
  return(result)
}