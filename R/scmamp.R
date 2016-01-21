# PACKAGE AND DATASET INFORMATION ----------------------------------------------

#' Statistical comparison of multiple algorithms
#'
#' This package has been develop to simplify the statistical assessment of algorithms when tested in different problems. It includes statistical tests, as well as some plotting functions.
#' @name scmamp
#' @author 
#' Borja Calvo \email{borja.calvo@@ehu.eus},
#' Guzman Santafe \email{guzman.santafe@@unavarra.es}
#' 
#' Maintainer: Borja Calvo \email{borja.calvo@@ehu.es}
#' @docType package
#' @aliases scmamp-package
#' @seealso For an overview of the use see #' \code{vignette(topic=
#' "Statistical_comparison_of_multiple_algorithms_in_multiple_problems", 
#' package="scmamp")} and \code{vignette(topic="Data_loading_and_manipulation", 
#' package="scmamp")}
NULL

#' Example in Garcia and Herrera (2008)
#'
#' Dataset corresponding to the accuracy of 5 classifiers in 30 datasets.
#' Each algorithm is in a column. This is the dataset used as example in 
#' Garcia and Herrera (2008).
#'
#' @format A data frame with 5 columns and 30 rows
#' @source S. Garcia and F. Herrera (2008) An Extension on "Statistical 
#' Comparisons of Classifiers over Multiple Data Sets" for all Pairwise 
#' Comparisons. \emph{Journal of Machine Learning Research}. 9, 2677-2694.
#' @name data.gh.2008
NULL

#' Example in Garcia and Herrera (2010)
#'
#' Dataset corresponding to the accuracy of 4 classifiers in 24 datasets.
#' Each algorithm is in a column. This is the dataset used as example in 
#' Garcia and Herrera (2010).
#'
#' @format A data frame with 4 columns and 24 rows
#' @source S. Garcia and F. Herrera (2010) Advanced Nonparametric Tests for 
#' Multiple Comparison in the Design of Experiments in Computational Intelligence 
#' and Data Mining: Experimental Analysis of Power. \emph{Information Sciences}, 
#' 180, 2044-2064.
#' @name data.gh.2010
NULL


#' Comparison of optimization algorithms in Blum \emph{et al.} (2015)
#'
#' This dataset contains part of the results obtained in the comparison of decentralyzed
#' optimization algorithms presented in Blum \emph{et al.} (2015). The dataset contains 
#' 900 rows and 10 colums. Each row reprsents an instance of the maximum independent
#' set problem (a graph). The first two are descriptors of the problem in each
#' row (size and radius used to create random geometric graphs) and the other 8
#' contain the results obtained by 8 algorithms for the MIS problem instance.
#' 
#'
#' @format A data frame with 10 columns and 900 rows
#' @source C. Blum, B. Calvo and M.J. Blesa (2015)  FrogCOL and FrogMIS: New Decentralized Algorithms for Finding Large Independent Sets in Graphs. \emph{Swarm Intelligence}. In press. 
#' @name data.blum.2015
NULL


# AUXILIAR FUNCTIONS -----------------------------------------------------------

correctPValues <- function(pvalues, correct) {
  # Auxiliar function to correct pvalues in either a vector (or a one row matrix) or 
  # a symmetric matrix where pvalues are repeated in both the upper and the lower half 
  # of the matrix
  # Args:
  #   pvalues:    Vector or matrix of the pvalues to correct
  #   correct:    Function to perform the correction
  #
  # Returns:
  #   The corrected p-values matrix
  #
  
  if (is.vector(pvalues) || nrow(pvalues)==1){
    res <- correct(pvalues)
  }else{
    k <- nrow(pvalues)
    pairs <- generatePairs(k, NULL)
    pval.vector <- pvalues[pairs]
    pval.vector.corrected <- correct(pval.vector)
    corrected.pvalues <- pvalues
    corrected.pvalues[pairs] <- pval.vector.corrected
    corrected.pvalues[pairs[,2:1]] <- pval.vector.corrected
    colnames(corrected.pvalues) <- colnames(pvalues)
    rownames(corrected.pvalues) <- rownames(pvalues)
    res <- corrected.pvalues
  }
  return(res)
}


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
    matrix.raw <- customPost(data=data, control=control, test=test, ...)
  }else{
    matrix.raw <- switch(test,
                         "t-test"= {
                           customPost(data=data, control=control, 
                                       test=function(x, y) {
                                         return(t.test(x, y, paired=TRUE))
                                         })
                       },
                       "wilcoxon"= {
                         customPost(data=data, control=control, 
                                    test=function(x,y) {
                                      return(wilcoxonSignedTest(x,y))
                                    })
                       },
                       "friedman"= {
                         friedmanPost(data=data, 
                                      control=control)
                       },
                       "aligned ranks"= {
                         friedmanAlignedRanksPost(data=data, 
                                                  control=control)
                       },
                       "quade"= {
                         quadePost(data=data, 
                                   control=control)
                       },
                       "tukey"= {
                         tukeyPost(data=data, 
                                   control=control)
                       },
                       stop("Unknown test. Valid options in the current version ",
                            "are 'wilcoxon', 'friedman', 'aligned ranks', 'quade', ",
                            "and 't-test'. Alternatively, you can pass a function ",
                            "that performs a paired statistical texts which ",
                            "should have, at least, two parameters, 'x' and 'y' ",
                            "and returns the p-value associted to the comparison"))
  }

  return(matrix.raw)
}
# EXPORTED FUNCTIONS -----------------------------------------------------------

#' @title Post hoc tests for multiple comparison analises
#'
#' @description This function is a wrapper to run the post hoc tests. It can run both all vs. control and all vs. all post hoc tests.
#' @param data A matrix or data frame containing the results obtained by the algorithms (columns) in each problem (rows). It can contain additional columns, but if any of the column has to be discarderd (not used neither to group the problems nor to be part of the comparison), then it is mandatory to indicate, in the \code{algorithms} parameter, which columns contain the algorithm information.
#' @param algorithms Vector with either the names or the indices of the columns that contain the values to be tested. If not provided, the function assumes that all the columns except those indicated in \code{group.by} represent the results obtained by an algorithm.
#' @param group.by Vector with either the names or the indices of the columns to be used to group the data. Each group is tested independently. If \code{NULL}, all the data is used for a single comparison.
#' @param test Parameter that indicates the statistical test to be used. It can be either a string indicating one of the available test or a function. As a string, it can take the following values:
#'  \itemize{
#'    \item {\code{'wilcoxon'} - Wilcoxon Signed Rank test, as in Demsar (2006)}
#'    \item {\code{'t-test'} - t-test (R's t.test function with paired option set at \code{TRUE})}
#'    \item {\code{'friedman'} - Friedman post hoc test, as in Demsar (2006)}
#'    \item {\code{'aligned ranks'} Friedman's Aligned Ranks post hoc test, as in Garcia and Herrera (2010)}
#'    \item {\code{'quade'} - Quade post hoc test, as in Garcia and Herrera (2010)} 
#'    \item {\code{'tukey'} - Tukey's ANOVA post hoc test, as in Test 28 in Kanji (2006).}
#'  }
#'  
#'  If a function is provided, then it has to have as first argument a matrix containing the columns to be compared. The function has to return a  list with, at least, an element named \code{p.value} (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param control Either the name or the index of a column in the dataset (one of those in the \code{algorithms} vector), to be used as control. Alternatively, this argument can be \code{'min'}, to select the algorithm with the minimum value, \code{'max'}, to select the algorithm with the maximum value as control. If the argument is not provided (or is \code{NULL}), all the pairwise comparisons are performed instead of all vs. control comparisons.
#' @param use.rank If \code{TRUE}, then the summarization of the data is based on the ranks, rather than on the actual values. The selecion of the algorithm with the maximum or minimum value is also done in terms of the summarized ranking.
#' @param sum.fun Function to be used to summarize the data. By default, average is used.
#' @param correct Either string indicating the type of correction that has to be applied or a function to correct the p-values for multiple testing; This parameter is only need in case the data is grouped. As a string, the valid values are:
#' \itemize{
#'   \item{\code{shaffer} - Shaffer's (static) procedure, as in Garcia and Herrera (2008)}
#'   \item{\code{bergmann} - Bergman and Hommel's  procedure (similar to Shaffer dynamic), as in Garcia and Herrera (2008)}
#'   \item{\code{holland} - Holland's procedure, as in Garcia and Herrera (2010)}
#'   \item{\code{finner} - Finner's procedure, as in Garcia and Herrera (2010)}
#'   \item{\code{rom} - Rom's procedure, as in Garcia and Herrera (2010)} 
#'   \item{\code{li} - Li's procedure, as in Garcia and Herrera (2010)}
#'   \item{Any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}}
#' }. 
#' If a function is provided, the it has to recieve, as first argument, a vector of pvalues to be corrected and has to return a verctor with the corrected p-values \emph{in the same order} as the input vector.
#' @param alpha Alpha value used in Rom's correction. By default, it is set at 0.05.
#' @param ... Special argument used to pass additional parameters to the statistical test and the correction method.
#' @return In all cases the function returns a list with three elements, the summarization of the data (a row per group), the raw p-values and the corrected p-values. When the data is grouped and all the pairwise comparisons are performed (no control is provided), the p-values are in three dimensional arrays where the last dimension is corresponds to the group. In any other cases the result is a matrix with one or more rows.
#' 
#' Note that Shaffer and Bergmann and Hommel's correction can only be applied when all the pairwise tests are conducted, due to their assumptions. Moreover, its use when the data is grouped (multiple pairwise comparsions) is not trivial and, thus, it is not possible to use it when the data is grouped.
#' 
#' @seealso \code{\link{friedmanPost}}, \code{\link{friedmanAlignedRanksPost}}, \code{\link{quadePost}}, \code{\link{tukeyPost}}, \code{\link{adjustShaffer}}, \code{\link{adjustBergmannHommel}}, \code{\link{adjustHolland}}, \code{\link{adjustFinner}}, \code{\link{adjustRom}}, \code{\link{adjustLi}}
#' 
#' @references S. Garcia and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references Garcia S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for All Pairwise Comparisons. \emph{Journal of Machine Learning Research}, 9, 2677-2694.
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' 
#' @examples
#' # Grouped data, all pairwise
#' data(data_blum_2015)
#' res <- postHocTest (data=data.blum.2015, algorithms=c("FrogCOL", "FrogMIS", "FruitFly"), 
#'                     use.rank=TRUE, group.by=c("Size"), test="quade", correct="finner")
#'                    
#' # Data summarization
#' res$summary
#' 
#' # Corrected pvalues for the first group
#' res$corrected.pval[, , 1]
#' 
#' # Grouped data, all vs. control
#' res <- postHocTest (data=data.blum.2015, control="max", use.rank=FALSE, 
#'                     group.by=c("Size","Radius"), test="wilcoxon", correct="finner")
#'                    
#' # Data summarization
#' res$summary
#' 
#' # Corrected pvalues
#' res$corrected.pval
#'                                        
#' # Not grouped data
#' data(data_gh_2008)
#' postHocTest (data=data.gh.2008, test="aligned ranks", correct="bergmann")
#' 
postHocTest <- function (data, algorithms=NULL, group.by=NULL, test="friedman", 
                         control=NULL, use.rank=FALSE, sum.fun=mean, 
                         correct="finner", alpha=0.05, ... ) {
 
  #If data is a matrix, convert it to a data.frame to avoid problems with the use of further methods
  if(is.matrix(data)){
    data <- data.frame(data)
  }
  
  # If there are only two algorithms all vs. all approach is, actually equalt to
  # all vs control. To avoid problems trying to generate matrix with only one row
  # we will convert take the second algorithm as the control
  if (length(algorithms) == 2 & is.null(control)) {
    control <- algorithms[2]
  }
  
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
    if (is.character(correct) & (correct == "shaffer" | correct == "bergmann")) {
      stop("Shaffer's and Bergman and Hommel's correction can only be used ",
              "when all the pairs are compared. For comparisons with a control ",
              "use any of the other corrections that do not take into account ",
              "the particular nature of all pairwise comparisons.")
    }
    if (is.numeric(control)) {
      control <- names(data)[control]
    }
  }  
  
  if (!is.null(group.by) & (correct == "shaffer" | correct == "bergmann")) {
    stop("Shaffer's and Bergmann and Hommel's corrections cannot be used with grouped data.",
         " Please select another correction method")
  }
  
  # Prepare the correction
  if(is.character(correct)) {
    correct.name <- correct
    correct <- switch (correct,
                       "shaffer"=adjustShaffer,
                       "bergmann"=adjustBergmannHommel,
                       "holland"={
                         function(pvalues) {
                           fun <- adjustHolland
                           return(correctPValues(pvalues=pvalues,
                                                         correct=fun))
                         }
                       },
                       "finner"={
                         function(pvalues) {
                           fun <- adjustFinner
                           return(correctPValues(pvalues=pvalues,
                                                         correct=fun))
                         }
                       },
                       "rom"={
                         function(pvalues) {
                           
                           fun <- function(pvalues) {
                             adjustRom(pvalues=pvalues, alpha=alpha)
                           }
                           return(correctPValues(pvalues=pvalues,
                                                         correct=fun))
                         }
                       },
                       "li"={
                         function(pvalues) {
                           
                           fun <- adjustLi
                           return(correctPValues(pvalues=pvalues,
                                                         correct=fun))
                         }
                       },
                       {
                         if (!(correct %in% p.adjust.methods)){
                           stop("Non valid method for p.adjust function. ",
                                "Valid options are ",
                                paste(p.adjust.methods, collapse="; "),
                                ". Additionally, 'holland', 'finner', 'rom' and ",
                                "'li' are also valid options")
                         }
                         function(pvalues) {
                           fun <- function(pvalues) {
                             p.adjust(p=pvalues, method=correct.name)
                           }
                           return(correctPValues(pvalues=pvalues,
                                                         correct=fun))
                         }
                       })
  } else {
    correct.name <- deparse(substitute(test))  
  }

  # Build the summary matrix
  if (is.null(group.by)){
   if (use.rank){
     aux <- rankMatrix(data=data[, algorithms], ...)
   }else{
     aux <- data[, algorithms]
   }
   # Note that in aux the group.by columns are at the begining
   # Some operations may change the name of the algorithms (special characters)
   sum.matrix <- summarizeData(data=aux, fun=sum.fun, group.by=NULL, ...)
   ## when the data is not summaryced by groups summarizeData returns a verctor of numeric values.
   ## As this may produce problems in other method, we force the sum.matrix to be a matrix
   if(is.numeric(sum.matrix)){
     sum.matrix <- t(as.matrix(sum.matrix))
   }
   colnames(sum.matrix) <- colnames(data)[algorithms]
  } else {
    if (use.rank) {
      aux <- cbind(data[, group.by], data.frame(rankMatrix(data=data[, algorithms])), ...)
    } else {
      aux <- cbind(data[, group.by], data[, algorithms])
    }
    # Note that in aux the group.by columns are at the begining
    sum.matrix <- summarizeData(data=aux, fun=sum.fun, group.by=1:length(group.by), ...)
    # Some operations may change the name of the algorithms (special characters)
    colnames(sum.matrix) <- names(data)[c(group.by, algorithms)]
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
      group.result <- runPostHoc (data.sub, test=test, control=control, ...)
      return(group.result)
    }
    
    group.raw.pval <- unlist(lapply (1:nrow(groups), FUN=getRawPvalues))
    
    # Now we create either the arrays or the matrix, depending on whether we have
    # a control or not
    k <- length(algorithms)
    p <- nrow(groups)
    if (is.null(control)){
      dim(group.raw.pval) <- c(k, k, p)
      # Now correct the pvalues. Note that we have to take only the upper part of each matrix
      # We will create the index triplets we need for selecting the correct pvalues
      pairs <- generatePairs(k, NULL)
      triplets <- cbind(pairs[rep(1:nrow(pairs), p), ], 
                        unlist(lapply(1:p,
                                      FUN=function(i) {
                                        return (rep(i, nrow(pairs)))
                                      })))
      corrected.pvalues <- correct(group.raw.pval[triplets])
      group.corrected.pval <- group.raw.pval
      group.corrected.pval[triplets] <- corrected.pvalues
      group.corrected.pval[triplets[,c(2,1,3)]] <- corrected.pvalues
      
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
      ## Correct with all the values
      group.corrected.pval <- correct(group.raw.pval)
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
#' @description This function is a wrapper to multiple comparison tests.
#' @param data A matrix or data frame containing the results obtained by the algorithms (columns) in each problem (rows). It can contain additional columns, but if any of the column has to be discarderd (not used neither to group the problems nor to be part of the comparison), then it is mandatory to indicate, in the \code{algorithms} parameter, which columns contain the algorithm information.
#' @param algorithms Vector with either the names or the indices of the columns that contain the values to be tested. If not provided, the function assumes that all the columns except those indicated in \code{group.by} represent the results obtained by an algorithm.
#' @param group.by Vector with either the names or the indices of the columns to be used to group the data. Each group is tested independently. If \code{NULL}, all the data is used for a single comparison.
#' @param test Parameter that indicates the statistical test to be used. It can be either a string indicating one of the available test or a function. As a string, it can take the following values:
#'  \itemize{
#'    \item {\code{'friedman'} - Friedman test, as in Garcia and Herrera (2010)}
#'    \item {\code{'aligned ranks'} Friedman's Aligned Ranks test, as in Garcia and Herrera (2010)}
#'    \item {\code{'quade'} - Quade test, as in Garcia and Herrera (2010)} 
#'    \item {\code{'anova'} - ANOVA test, as in Test 22 in Kanji (2006).}
#'  }
#'  
#'  If a function is provided, then it has to have as first argument a matrix containing the columns to be compared. The function has to return a  list with, at least, an element named \code{p.value} (as the \code{htest} objects that are usually returned by R's test implementations).
#' @param correct Either string indicating the type of correction that has to be applied or a function to correct the p-values for multiple testing; This parameter is only need in case the data is grouped. As a string, the valid values are:
#' \itemize{
#'   \item{\code{holland} - Holland's procedure, as in Garcia and Herrera (2010)}
#'   \item{\code{finner} - Finner's procedure, as in Garcia and Herrera (2010)}
#'   \item{\code{rom} - Rom's procedure, as in Garcia and Herrera (2010)} 
#'   \item{\code{li} - Li's procedure, as in Garcia and Herrera (2010)}
#'   \item{Any of the methods implemented in the \code{p.adjust} function. For a list of options, type \code{p.adjust.methods}}
#' }. 
#' If a function is provided, the it has to recieve, as first argument, a vector of pvalues to be corrected and has to return a verctor with the corrected p-values \emph{in the same order} as the input vector.
#' @param alpha Alpha value used in Rom's correction. By default, set at 0.05.
#' @param ... Special argument used to pass additional parameters to the statistical test and the correction method.
#' @return In case the \code{group.by} argument is not provided (or it is \code{NULL}), the function return an object of class \code{htest}. If columns for grouping are provided, then the function returns a matrix that includes, for each group, the values of the \code{group.by} columns, the raw p-value and the corrected p-value.
#' 
#' #' @seealso \code{\link{friedmanTest}}, \code{\link{friedmanAlignedRanksTest}}, \code{\link{quadeTest}}, \code{\link{anovaTest}}, \code{\link{adjustShaffer}}, \code{\link{adjustBergmannHommel}}, \code{\link{adjustHolland}}, \code{\link{adjustFinner}}, \code{\link{adjustRom}}, \code{\link{adjustLi}}
#' 
#' @references S. Garcia and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' 
#' @examples
#' # Grouped data
#' data(data_blum_2015)
#' multipleComparisonTest (data=data.blum.2015, 
#'                         algorithms=c("FrogCOL", "FrogMIS", "FruitFly"), 
#'                         group.by=c("Size", "Radius"), 
#'                         test="quade", correct="finner")
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
  if (length(algorithms) < 3) {
    stop("At least three algorithms are required to run the function")
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
                    {
                      stop("Unknown test. Valid options are 'friedman', ",
                           "'aligned ranks', 'quade', 'anova' or a function that ",
                           "gets as input a data.frame or matrix where each ",
                           "algorithm is in a column")
                    })
  } else {
    test.name <- deparse(substitute(test))  
  }
  
  # Prepare the correction function
  if(is.character(correct)) {
    correct.name <- correct
    correct <- switch (correct,
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


#' @title Contrast estimation based on medians
#'
#' @description This function performs estimates the contrast between algorithms through the medians
#' @param data Matrix or data frame with the data to compare
#' @return A matrix where the estimation of all the pairs of differences are output. 
#' The differences correspond to row-column. 
#' @details The test has been implemented according to Garcia \emph{et al.} (2010), Section 3.3.
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @examples
#' data(data_gh_2008)
#' contrastEstimationMatrix(data.gh.2008)

contrastEstimationMatrix <- function (data) {
  k <- dim(data)[2]
  pairs <- generatePairs(k=k, control=NULL)
  medians <- apply(pairs, MARGIN=1,
                   FUN=function(x) {
                     differences <- data[, x[1]] - data[, x[2]]
                     return(median(differences))
                   })
  median.matrix <- matrix(rep(0, k^2), ncol=k)
  median.matrix[pairs] <- medians
  median.matrix[pairs[,c(2,1)]] <- -1*medians
  adjusted.m <- rowMeans(median.matrix)
  
  adjusted.differences <- apply(pairs, MARGIN=1,
                                FUN=function(x) {
                                  differences <- adjusted.m[x[1]] - adjusted.m[x[2]]
                                  return(median(differences))
                                })
  adjusted.matrix <- matrix(rep(0, k^2), ncol=k) 
  adjusted.matrix[pairs] <- adjusted.differences
  adjusted.matrix[pairs[,c(2,1)]] <- -1 * adjusted.differences
  colnames(adjusted.matrix) <- colnames(data)
  rownames(adjusted.matrix) <- colnames(data)
  
  return(adjusted.matrix)
}