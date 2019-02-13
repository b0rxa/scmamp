# NON-EXPORTED, AUXILIAR FUNCTIONS --------------------------------------------

processControlColumn <- function(data, control){
  # Function to process the control argument used in some functions in this script
  #
  # Args:
  #   data:    Dataset where the column is
  #   control: Name or ID of the column
  #
  # Returns:
  #   The ID of the column. If there is any problem, an error is rised.
  #
  if(!is.null(control)) {
    if (is.character(control)) {
      control <- which(names(data) %in% control)
      if (length(control)==0) {
        stop("The name of the column to be used as control does not exist ",
             "in the data matrix")
      }
    } else {
      if (control > ncol(data) | control < 1) {
        stop("Non-valid value for the control parameter. It has to be either ",
             "the name of a column or a number between 1 and ", ncol(data))
      }
    }
  }
  return(control)
}

generatePairs <- function(k, control) {
  # Non-exported function to create the pairs needed for the comparisons
  #
  # Args:
  #  k:       Number of algorithms.
  #  control: Id of the control. It can be NULL to indicate that all pairs have
  #           to be created.
  #
  # Returns:
  #   A data.frame with all the pairs to be compared
  #
  
  if(is.null(control)) { 
    # All vs all comparison
    aux.pairs <- sapply(1:(k - 1), 
                        FUN=function(x) {
                          cbind((x),(x+1):k)
                        })
    pairs <- do.call(what=rbind, args=aux.pairs)
  }  else{ 
    # All vs control comparison
    pairs <- cbind(control, (1:k)[-control])
  }
  return(pairs)
}

buildPvalMatrix <- function(pvalues, k, pairs, cnames, control){
  # Function to create the final p-value matrix
  #
  # Args:
  #   pvalues: Vector of p-values.
  #   k:       Number of algorithms.
  #   pairs:   Pair of algorithms to which each p-value correspond.
  #   cnames:  Character vector with the names of the algorithms.
  #   control: Id of the algorithm used as control (if any). Only usedt to know
  #            the type of comparison.
  #
  # Returns:
  #   A correctly formated matrix with the results.
  #
  if(is.null(control)){ 
    # All vs all comparison
    matrix.raw <- matrix(rep(NA, times=k^2), nrow=k)
    # The matrix has to be symetric
    matrix.raw[pairs]            <- pvalues
    matrix.raw[pairs[, c(2, 1)]] <- pvalues
    colnames(matrix.raw) <- cnames
    rownames(matrix.raw) <- cnames
  } else { 
    # All vs control comparison
    matrix.raw <- matrix(rep(NA, k), ncol=k)
    matrix.raw[(1:k)[-control]] <- pvalues
    colnames(matrix.raw) <- cnames    
  }  
  return(matrix.raw)
}


correctForMonotocity <- function (pvalues){
  # Function to ensure the monotocity of the p-values after being corrected
  # Args:
  #   pvalues: Corrected p-values, in DECRESING ORDER ACCORDING TO THE RAW PVALUE
  # Returns:
  #   The corrected p-values such as there is not a p-value smaller than the any 
  #   of the previous ones.
  #
  pvalues <- sapply(1:length(pvalues), 
                    FUN=function(x) {
                      return(max(pvalues[1:x]))
                    })
  return(pvalues)
}

setUpperBound <- function (vector, bound){
  # Simple function to limit the maximum value of a vector to a given value
  #
  # Args:
  #   vector: Vector whose values will be bounded
  #   bound:  Maximum value in the result
  #
  # Returns:
  #   A vector equal to 'vector' except for those values above 'bound', that are
  #   set to this upper bound
  #
  res <- sapply(vector, 
                FUN=function(x) {
                  return(min(bound, x))
                })
  return(res)
}


countRecursively <- function(k) {
  # Auxiliar function to get the S(k) decribed in Shaffer (1985)
  # Args:
  #   k: Number of algorithms
  #
  # Returns:
  #   List of the maximum number of true hypothesis in a pairwise comparsison of k classifiers
  #
  res <- c(0)
  if (k > 1) {
    res <- c(res, countRecursively(k - 1))
    for (j in 2:k) {
      res <- c(res, countRecursively(k - j) + (factorial(j) / (2 * factorial(j - 2))))
    }
  }
  return(sort(unique(res)))
}


computeSubdivisions <- function (set){
  # Auxiliar function to get all the possible subdivisions of a set
  # Args:
  #   set: The set of element whose subdivisions will be obtained
  #
  # Returns:
  #   All the posible subdivisions of the set as a list. Each element is a list 
  #   of two vectors, corresponding to a subdivision.
  #
  if (length(set) == 1) {  # The trivial case
    res <- list(list(s1=set, s2 = vector()))
  } else {  # In the general case, subdivide the set without the last element and then add ...
    n    <- length(set)
    last <- set[n]
    sb   <- computeSubdivisions (set[-n])
    # ... the last element in all the first sets ...
    res <- lapply(sb, 
                  FUN=function (x) {
                    return(list(s1=c(x$s1, last), s2=x$s2))
                  })
    # ... the last element in all the second sets ...
    res <- c(res, 
             lapply(sb, 
                    FUN=function (x) {
                      return(list(s1=x$s1, s2=c(x$s2, last)))
                    }))
    # ... and finally the subdivision that has the last element alone in s2
    res <- c(res, list(list(s1=last, s2=set[-n])))
  }
  return(res)
}


partition <- function(set) {
  # Auxiliar function to compute the complete set of non-empty partitions
  # Args:
  #   set: Set to get the partitions
  # 
  # Returns:
  #   A list with the complete set of partitions
  #
  # Details:
  #   This implementation has to do with the 6th step in the algorithm shown in
  #   Figure 1 in Garcia and Herrera (2008).
  #
  n <- length(set)
  
  if (n == 1){ 
    # In the trivial case the function returns the only possible subdivision. 
    # This case clearly violates the idea of having the last element in the last 
    # set (used to avoid empty sets), but it is the only exception and it does not 
    # pose a problema as the repetition of the set is under control.
    res <- computeSubdivisions(set)
  } else {
    # In the general case, get all the subdivisions (which cannot have empty 
    # sets in the first subset) and add the last element at the end of each 
    # second subset.
    last <- set[n]
    sb   <- computeSubdivisions(set[-n])
    # Add the last element only to the second subset 
    # (see Figure 1 in Garcia and Herrera, 2008)
    res <- lapply(sb, 
                  FUN=function (x) {
                    return(list(s1=x$s1 , s2=c(x$s2 , last)))
                  })
  }
  return(res)
}

## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title Remove repetitions in a list of subsets
#
# @description This function removes the repetitions in the result returned by the function \code{\link{exhaustive.sets}}
# @param E list of exhaustive sets
# @return The list without repetitions

megeSets <- function (e1, e2){
  # Auxiliar function to avoid repetitions in the results returned by
  # computeExhaustivSets
  #
  # Args:
  #   e1: First set
  #   e2: Second set
  # 
  # Returns:
  #   The union of the sets e1 and e2
  #  
  isIn <- function(e) {
    # Function to check wheter e is in e2 or not
    compareWithE <- function (en) {
      # Function to compare en with e
      eq <- FALSE
      if (all(dim(en) == dim(e))) {
        eq <- all(en == e)
      }
      return(eq)
    }    
    comparisons <- unlist(lapply(e2, FUN=compareWithE))
    return(any(comparisons))
  }
  
  # Sequentially add partitions if they are not already in the solution
  bk <- lapply (e1, 
                FUN=function (x) {
                  if (!isIn (x)) {
                    e2 <<- c(e2 , list(x))
                  }
                })
  return(e2)
}


# EXPORTED FUNCTIONS -----------------------------------------------------------

#' @title Friedman's post hoc raw p-values
#'
#' @description This function computes the raw p-values for the post hoc based on Friedman's test.
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples.
#' @param control Either the number or the name of the column for the control algorithm. If this parameter is not provided, the all vs all comparison is performed.
#' @param ... Not used. 
#' @return A matrix with all the pairwise raw p-values (all vs. all or all vs. control).
#' @details The test has been implemented according to the version in Demsar (2006), page 12.
#' @references J. Demsar (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @examples
#' data(data_gh_2008)
#' friedmanPost(data.gh.2008)
#' friedmanPost(data.gh.2008, control=1)

friedmanPost <- function (data, control=NULL, ...) {
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  # Parameter checking
  control <- processControlColumn(data, control)
  
  # Pairs to be tested
  pairs <- generatePairs(k, control)
  
  # Compute the p-values for each pair
  mean.rank <- colMeans(rankMatrix(data))
  sd <- sqrt((k * (k + 1)) / (6 * N))
  computePvalue  <- function(x) {
    stat <- abs(mean.rank[x[1]] - mean.rank[x[2]]) / sd
    # Two tailed test
    (1 - pnorm(stat))*2 
  }
  pvalues <- apply(pairs, MARGIN=1, FUN=computePvalue)
  
  # Create the result matrix, depending on the comparison performed
  matrix.raw <- buildPvalMatrix(pvalues=pvalues, k=k, pairs=pairs, 
                                cnames=colnames(data), control=control)
  
  return(matrix.raw)
}



#' @title Friedman's Aligned Ranks post hoc raw p-values
#'
#' @description This function computes the raw p-values for the post hoc based on Friedman's Aligned Ranks test.
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples.
#' @param control Either the number or the name of the column for the control algorithm. If this parameter is not provided, the all vs all comparison is performed.
#' @param ... Not used. 
#' @return A matrix with all the pairwise raw p-values (all vs. all or all vs. control).
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2010), pages 2051,2054
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @examples
#' data(data_gh_2008)
#' friedmanAlignedRanksPost(data.gh.2008)
#' friedmanAlignedRanksPost(data.gh.2008, control=1)

friedmanAlignedRanksPost <- function (data, control=NULL, ...) {
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  # Parameter checking
  control <- processControlColumn(data, control)
  
  # Pairs to be tested
  pairs <- generatePairs(k, control)
  
  # Compute the p-values for each pair
  dataset.means <- rowMeans(data)
  diff.matrix <- data - matrix(rep(dataset.means, k), ncol=k)
  # To get the rank of a decreasing ordering, change the sign of the differences
  ranks <- matrix(rank(-unlist(diff.matrix)), ncol=k, byrow=FALSE)
  colnames(ranks) <- colnames(data)
  mean.rank <- colMeans(ranks)
  
  sd <- sqrt((k * (N*k + 1)) / 6)
  computePvalue  <- function(x) {
    stat <- abs(mean.rank[x[1]] - mean.rank[x[2]]) / sd
    # Two tailed test
    (1 - pnorm(stat))*2 
  }
  pvalues <- apply(pairs, MARGIN=1, FUN=computePvalue)
  
  # Create the result matrix, depending on the comparison performed
  matrix.raw <- buildPvalMatrix(pvalues=pvalues, k=k, pairs=pairs, 
                                cnames=colnames(data), control=control)
  
  return(matrix.raw)
}


#' @title Quade post hoc raw p-values
#'
#' @description This function computes the raw p-values for the post hoc based on Quade's test.
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples.
#' @param control Either the number or the name of the column for the control algorithm. If this parameter is not provided, the all vs all comparison is performed.
#' @param ... Not used. 
#' @return A matrix with all the pairwise raw p-values (all vs. all or all vs. control).
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2010), pages 2052,2054
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparisons in the design of experiments in computational intelligence and ata mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @examples
#' data(data_gh_2008)
#' quadePost(data.gh.2008)
#' quadePost(data.gh.2008, control=1)
#' 

quadePost <- function (data, control=NULL, ...){
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  # Parameter checking
  control <- processControlColumn(data, control)
  
  #Pairs to be tested
  pairs <- generatePairs(k, control)
  
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
  
  computePvalue  <- function(x) {
    z <- abs(t.j[x[1]] - t.j[x[2]]) / nrm
    # Two tailed test
    (1 - pnorm(z)) * 2 
  }
  pvalues <- apply(pairs, MARGIN=1, FUN=computePvalue)
  
  # Create the result matrix, depending on the comparison performed
  matrix.raw <- buildPvalMatrix(pvalues=pvalues, k=k, pairs=pairs, 
                                cnames=colnames(data), control=control)
  
  return(matrix.raw)
}


#' @title Function to use custom tests to perform post hoc comparisons.
#'
#' @description This function computes the raw p-values for all vs. all or all vs. control comparisons using a custom function.
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples.
#' @param control Either the number or the name of the column for the control algorithm. If this parameter is not provided, the all vs all comparison is performed.
#' @param test Function to perform the test. It requires two parameters, \code{x} and \code{y}, the two samples to be compared, and it has to return a list that contains, at least, one element called p.value (as the \code{htest} objects that are usually returned by R's statistical test implementations).
#' @param ... Additional parameters for the test function.
#' @return A matrix with all the pairwise raw p-values.
#' @examples
#' data(data_gh_2008)
#' test <- function(x, y, ...) {
#'   t.test(x, y, paired=TRUE)
#' }
#' customPost(data.gh.2008, control=1, test=test)
#' customPost(data.gh.2008, test=test)

customPost <- function(data, control=NULL, test, ...){
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  # Parameter checking
  control <- processControlColumn(data, control)
  if (!is.function(test)) {
    stop("The 'test' argument has to be a function with at least two parameters",
         "x and y, corresponding to two vectors. It should return a list with, ",
         "at least, one element named p.value")
  }
  
  #Pairs to be tested
  pairs <- generatePairs(k, control)
  
  # Compute the p-values for each pair
  computePvalue  <- function(x) {
    res <- test(x=data[, x[1]], 
                y=data[, x[2]], ...)
    return(res$p.value)
  }
  pvalues <- apply(pairs, MARGIN=1, FUN=computePvalue)
  
  # Create the result matrix, depending on the comparison performed
  matrix.raw <- buildPvalMatrix(pvalues = pvalues, pairs = pairs, k=k,
                                cnames = colnames(data), control = control)
  
  return(matrix.raw)
}

#' @title Tukey post hoc test for ANOVA.
#'
#' @description This function computes all the pairwise p-values corrected using Tukey post hoc test.
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples.
#' @param control Either the number or the name of the column for the control algorithm. If this parameter is not provided, the all vs all comparison is performed.
#' @param ... Not used.
#' @return A matrix with all the pairwise corrected p-values.
#' @details The test has been implemented according to Test 22 in Kanji (2006).
#' @references G. K. Kanji (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @examples
#' data(data_gh_2008)
#' tukeyPost(data.gh.2008)
#' tukeyPost(data.gh.2008, control=1)

tukeyPost <- function (data, control=NULL, ...){
  
  # Implemented as Test 28 in 100 statistical tests
  k <- dim(data)[2]
  N <- dim(data)[1]
  nu <- k * (N - 1)
  
  
  # Generate all the pairs to test
  pairs <- generatePairs(k=k, control)
  
  # Sample variances
  var.vector <- apply(data, MARGIN=2, FUN=var)
  m.vector <- colMeans(data)
  
  # Total variance. Given that all the samples have the same size, it is the 
  # average variance
  var <- mean(var.vector)
  
  # Compute the p-value
  computePvalue <- function(x) {
    d <- abs(m.vector[x[1]] - m.vector[x[2]])
    q <- d * sqrt(N / var)
    1 - ptukey(q, k, nu)
  }
  pvalues <- apply(pairs, MARGIN=1, FUN=computePvalue)
  
  # Generate the final matrix with the p-values
  matrix.raw <- buildPvalMatrix(pvalues=pvalues, k=k, pairs=pairs, 
                                cnames=colnames(data), control=control)
  return(matrix.raw)
}

#' @title Shaffer's correction of p-values in pairwise comparisons.
#'.
#' @description This function implements the Shaffer's (static) multiple testing correction when the p-values correspond with pairwise comparisons.
#' @param raw.matrix A matrix with the pairwise p-values. The p-values have to be, at least, in the upper part of the matrix.
#' @return A symetric matrix with the corrected p-values.
#' @details The test has been implemented according to the version in Garcia and Herrera (2008), page 2680.
#' @references S. Garcia and F. Herrera (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for All Pairwise Comparisons. \emph{Journal of Machine Learning Research}, 9, 2677-2694.
#' @references J.P. Shaffer (1986) Modified sequentially rejective multiple test procedures. \emph{Journal of the American Statistical Association}, 81(395), 826-831.
#' 
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanPost(data.gh.2008)
#' raw.pvalues
#' adjustShaffer(raw.pvalues)

adjustShaffer <- function (raw.matrix){  
  if (!(is.matrix(raw.matrix) | is.data.frame(raw.matrix))) {
    stop("This correction method requires a square matrix or data.frame with ",
         "the p-values of all the pairwise comparisons.")
  }
  
  if (diff(dim(raw.matrix)) != 0) {
    stop("This correction method requires a square matrix or data.frame with ",
         "the p-values of all the pairwise comparisons.")
  }
  
  k <- dim(raw.matrix)[1]
  
  # Transform the matrix into a vector
  pairs <- generatePairs(k, NULL)
  raw.pvalues <- raw.matrix[pairs]
  
  sk <- countRecursively(k)[-1]
  # Get the number of hypothesis that can be simultaneously true for each
  # number of rejected hypothesis
  t.i <- c(rep(sk[-length(sk)], diff(sk)), sk[length(sk)])
  t.i <- rev(t.i)
  
  # Order the p-values to apply the correction and correct them with the number
  # of hypothesis that can be simultaneously true
  o <- order(raw.pvalues)
  adj.pvalues <- raw.pvalues[o] * t.i 
  adj.pvalues <- setUpperBound(adj.pvalues, 1)
  adj.pvalues <- correctForMonotocity(adj.pvalues)
  
  # Put the adjusted pvalues in the correct order and regenerate the matrix
  adj.pvalues <- adj.pvalues[order(o)]
  adj.matrix <- raw.matrix
  adj.matrix[pairs] <- adj.pvalues
  adj.matrix[pairs[,2:1]] <- adj.pvalues  
  
  return(adj.matrix)
}

#' @title Complete set of exhaustive sets.
#'
#' @description This function implements the algorithm in Figure 1, Garcia and Herrera (2008) to create, given a set, the complete set of exhaustive sets E.
#' @param set Set to create the exhaustive sets. The complexity of this algorithm is huge, so use with caution for sets of more than 7-8 elements. Indeed, the implementation, as it is, can be hardly used from sizes beyond 9.
#' @return A list with all the possible exhaustive sets, without repetitions.
#' @details The algorithm makes use of `exhaustive.sets`, a structure provided with the pacakge that contains the precomputed sets for size up to 9. With this structure the exhaustive sets are generated inmediately, but if the data is, for some reason, not loaded, the computation may take several hours (or even days, depending on the size of the set).
#' @examples
#' exhaustiveSets(c("A","B","C","D"))
#' 
exhaustiveSets <- function (set){
  k <- length(set)
  # Reuse the computed sets stored in the variable 'exhaustive.sets'
  if (!is.null(exhaustive.sets) & k <= length(exhaustive.sets)) { 
    es.l <- lapply(exhaustive.sets[[k]], 
                   FUN=function(x) {
                     return(matrix(set[x], nrow=2))
                   })
  } else {
    # All possible pairwise comparisons, Garcia and Herrera, Table 1 steps 1-5
    if (k == 1) {
      es.l <- NULL
    } else if (k==2) { ## Base case, no lists
      es.l <- list(matrix(c(set[1], set[2]), ncol=1))
    } else {
      pairs <- generatePairs(k, control=NULL)
      es.l <- list(apply (pairs, MARGIN=1, 
                          FUN=function(x) {
                            return(c(set[x[1]], set[x[2]]))
                          }))
    }
    # Main loop
    if (length(set) > 2) {
      partitions <- partition(set)  # Sets in step 6
      processPartition <- function (x) { # Function to perform the main loop
        es1 <- exhaustiveSets(x$s1)
        es2 <- exhaustiveSets(x$s2)
        if (!is.null(es1)) {
          es.l <<- megeSets(es1, es.l)
          if (!is.null(es2)) {
            lapply (es1, 
                    FUN=function(e1) {
                      es.aux <- lapply(es2,
                                       FUN=function(e2) {
                                         return(cbind(e1,e2))
                                       })
                      es.l <<- megeSets(es.aux, es.l)
                    })
          }
        }
        if (!is.null(es2)) es.l <<- megeSets(es2, es.l)
      }
      # Do the loop (bk is just to avoid printing trash in the screen ...)
      bk <- lapply(partitions, FUN=processPartition)
    }
  }
  gc()
  return(es.l)
}

#' @title Bergmann and Hommel dynamic correction of p-values.
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Bergmann and Hommel (1994).
#' @param raw.matrix Raw p-values in a matrix.
#' @return A matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia and Herrera (2008), page 2680-2682.
#' @references S. Garcia and F. Herrera (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for All Pairwise Comparisons. \emph{Journal of Machine Learning Research}, 9, 2677-2694.
#' @references G. Bergmann and G. Hommel (1988) Improvements of general multiple test procedures for redundant systems of hypogheses. In P. Bauer, G. Hommel and E. Sonnemann, editors, \emph{Multiple Hypotheses Testing}, 100-115, Springer, Berlin.
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanAlignedRanksPost(data.gh.2008)
#' raw.pvalues
#' adjustBergmannHommel (raw.pvalues)
#' 
adjustBergmannHommel <- function (raw.matrix){
  if (!(is.matrix(raw.matrix) | is.data.frame(raw.matrix))) {
    stop("This correction method requires a square matrix or data.frame with ",
         "the p-values of all the pairwise comparisons.")
  }
  
  if (diff(dim(raw.matrix)) != 0) {
    stop("This correction method requires a square matrix or data.frame with ",
         "the p-values of all the pairwise comparisons.")
  }
  
  k <- dim(raw.matrix)[1]
  # The exhaustive.sets is a pre-computed global variable. Computing it every time the function is called is
  # unaffordable
  if(k > length(exhaustive.sets)) {
    stop ("Sorry, this method is only available for", 
          length(exhaustive.sets),
          "or less algorithms.")
  }
  
  if (k==9) {
    message("Applying Bergmann Hommel correction to the p-values computed in ",
            "pairwise comparisions of 9 algorithms. This requires checking", 
            length(exhaustive.sets$k9), "sets of hypothesis. It may take a few seconds.")
  }
  es.k <- exhaustive.sets[[k]]
  
  pairs <- generatePairs(k, control=NULL)
  raw.pvalues <- raw.matrix[pairs]
  correct <- function(i){
    aux <- lapply(es.k, 
                  FUN = function(s) {
                    # check that i is in s
                    if (any(colSums(pairs[i,] == s) == 2)) {
                      nu = dim(s)[2] * min(raw.matrix[t(s)])
                      return(min(nu, 1))
                    }
                  })
    return(max(unlist(aux)))
  }
  
  adj.pvalues <- sapply(1:dim(pairs)[1], FUN=correct)
  
  # Correct any possible inversion
  o <- order(raw.pvalues)
  aux <- setUpperBound(adj.pvalues[o], 1)
  aux <- correctForMonotocity(aux)
  adj.pvalues <- aux[order(o)]
  
  # Build the final matrix
  adj.matrix <- raw.matrix
  adj.matrix[pairs] <- adj.pvalues
  adj.matrix[pairs[,2:1]] <- adj.pvalues  
  
  return(adj.matrix)
}


#' @title Holland correction of p-values.
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Holland and Copenhaver (1987)
#' @param pvalues Raw p-values in either a vector or a matrix. Note that in case the p-values are in a matrix, all the values are used for the correction. Therefore, if the matrix contains repeated values (as those output by some methods in this package), the repetitions have to be removed.
#' @return A vector or matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2010), page 2680-2682.
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera, F. (2010) Advanced nonparametric tests for multiple comparison in the design of experiments in computational intelligence and data mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references B. S. Holland and M. D. Copenhaver (1987) An improved sequentially rejective Bonferroni test procedure \emph{Biometrics}, 43, 417-423.
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanPost(data.gh.2008)
#' raw.pvalues
#' adjustHolland (raw.pvalues)
#' 
adjustHolland <- function(pvalues) {
  ord <- order(pvalues, na.last=NA)
  pvalues.sorted <- pvalues[ord]
  k <- length(pvalues.sorted) + 1
  
  p.val.aux <- sapply(1:(k - 1), 
                      FUN=function(j, p.val, k){
                        r <- 1 - (1 - p.val[j])^(k - j)
                        return(r)
                      },
                      p.val=pvalues.sorted, k=k)
  
  p.adj.aux <- setUpperBound(p.val.aux, 1)
  p.adj.aux <- correctForMonotocity(p.adj.aux)
  
  p.adj <- rep(NA, length(pvalues))
  suppressWarnings(expr = {
                     p.adj[ord] <- p.adj.aux
                   })
  
  if(is.matrix(pvalues)){
    p.adj <- matrix(p.adj, ncol=ncol(pvalues))
    colnames(p.adj) <- colnames(pvalues)
    rownames(p.adj) <- rownames(pvalues)
  }
  return(p.adj)
}

#' @title Finner correction of p-values
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Finner (1993)
#' @param pvalues Raw p-values in either a vector or a matrix. Note that in case the p-values are in a matrix, all the values are used for the correction. Therefore, if the matrix contains repeated values (as those output by some methods in this package), the repetitions have to be removed.
#' @return A vector or matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2010), page 2680-2682.
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparison in the design of experiments in computational intelligence and data mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references H. Finner (1993) On a monotocity problem in ste-down mulitple test procedures. \emph{Journal of the American Statistical Association}, 88, 920-923.
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanPost(data.gh.2008)
#' raw.pvalues
#' adjustFinner (raw.pvalues)
#' 
adjustFinner <- function(pvalues) {
  ord <- order(pvalues, na.last=NA)
  pvalues.sorted <- pvalues[ord]
  k <- length(pvalues.sorted) + 1
  
  p.val.aux <- sapply(1:(k - 1),
                      FUN=function(j, p.val, k) {
                        r <- 1 - (1 - p.val[j])^((k - 1) / j)
                        return(r)
                      },
                      p.val=pvalues.sorted, k=k)
  
  p.adj.aux <- setUpperBound(p.val.aux, 1)
  p.adj.aux <- correctForMonotocity(p.adj.aux) 
  
  p.adj <- rep(NA, length(pvalues))
  suppressWarnings(expr = {
                     p.adj[ord] <- p.adj.aux
                   })
  
  if(is.matrix(pvalues)) {
    p.adj <- matrix(p.adj, ncol=ncol(pvalues))
    colnames(p.adj) <- colnames(pvalues)
    rownames(p.adj) <- rownames(pvalues)
  }
  return(p.adj)
}


#' @title Rom correction of p-values
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Rom (1990)
#' @param pvalues Raw p-values in either a vector or a matrix. Note that in case the p-values are in a matrix, all the values are used for the correction. Therefore, if the matrix contains repeated values (as those output by some methods in this package), the repetitions have to be removed.
#' @param alpha value for the averall test
#' @return A vector or matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia\emph{et al.} (2010), page 2680-2682.
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparison in the design of experiments in computational intelligence and data mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references D. M. Rom (1990) A sequentially rejective test procedure based on a modified Bonferroni inequality. \emph{Biometrika}, 77, 663-665.
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanPost(data.gh.2008)
#' raw.pvalues
#' adjustRom(raw.pvalues, alpha=0.05)
adjustRom <- function(pvalues, alpha=0.05){
  ord <- order(pvalues, na.last=NA)
  pvalues.sorted <- pvalues[ord]
  k <- length(pvalues.sorted) + 1  
  
  alpha.adj <- numeric(k - 1)
  alpha.adj[k - 1] <- alpha
  alpha.adj[k - 2] <- alpha/2  
  for(i in 3:(k - 1)){
    j1 <- 1:(i - 1)
    j2 <- 1:(i - 2)
    aux <- choose(i,j2) * (alpha.adj[k - 1 - j2])^(i - j2)
    alpha.adj[k-i] <- (sum(alpha^j1) - sum(aux)) / i      
  }
  
  r <- alpha / alpha.adj
  
  p.val.aux <- r * pvalues.sorted
  p.val.aux <- setUpperBound(p.val.aux, 1)
  p.adj.aux <- correctForMonotocity(pvalues=p.val.aux)
  
  p.adj.aux <- sapply((k - 1):1,
                      FUN=function(i) {
                        min(p.adj.aux[(k-1):i])
                      })
  
  # Note that the code above ensures that any value is not bigger than any of the
  # values bellow, but it reverse the order of the corrected p-values.
  
  p.adj.aux <- rev(p.adj.aux)
  
  p.adj <- rep(NA,length(pvalues))
  suppressWarnings(expr = {
                     p.adj[ord] <- p.adj.aux
                   })
  
  if(is.matrix(pvalues)){
    p.adj <- matrix(p.adj, ncol=ncol(pvalues))
    colnames(p.adj) <- colnames(pvalues)
    rownames(p.adj) <- rownames(pvalues)
  }
  return(p.adj)
}

#' @title Li correction of p-values.
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Li (2008).
#' @param pvalues Raw p-values in either a vector or a matrix. Note that in case the p-values are in a matrix, all the values are used for the correction. Therefore, if the matrix contains repeated values (as those output by some methods in this package), the repetitions have to be removed.
#' @return A vector or matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia \emph{et al.} (2010), page 2680-2682. This is a simple procedure that provides good results when the highest p-value corrected is below 0.5. However, When the highest p-value is close to 1 the correction is extremely conservative. Actually, when the highest p-value is 1, all the corrected p-values are set at 1. Therefore, it is not advisable to be used under these circumstances.
#' @references S. Garcia, A. Fernandez, J. Luengo and F. Herrera (2010) Advanced nonparametric tests for multiple comparison in the design of experiments in computational intelligence and data mining: Experimental analysis of power. \emph{Information Sciences}, 180, 2044-2064.
#' @references J. Li (2008) A two-step rejection prcedure for testing mulitple hypotheses. \emph{Journal of Statistical Planning and Inference}, 138, 1521-1527.
#' @examples
#' data(data_gh_2008)
#' raw.pvalues <- friedmanPost(data.gh.2008)
#' adjustLi(raw.pvalues)
adjustLi <- function(pvalues){
  ord <- order(pvalues, na.last=NA)
  pvalues.sorted <- pvalues[ord]
  k <- length(pvalues.sorted) + 1   
  if (max(pvalues, na.rm=TRUE) > 0.5)
    warning("The highest p-value is above 0.05. In such a situation the method is ",
            "far too conservative, so consider using another method.")
  
  
  p.adj.aux <- sapply(1:(k - 1),
                      FUN=function(i, ps) {
                        l <- length(ps)
                        return(min(ps[l], ps[i] / (ps[i] + 1 - ps[l])))
                      },
                  ps=pvalues.sorted)    
  
  p.adj <- rep(NA, length(pvalues))
  suppressWarnings(expr = {
                     p.adj[ord] <- p.adj.aux
                   })
  if(is.matrix(pvalues)) {
    p.adj <- matrix(p.adj, ncol=ncol(pvalues))  
    colnames(p.adj) <- colnames(pvalues)
    rownames(p.adj) <- rownames(pvalues)
  }
  
  return(p.adj)
}
