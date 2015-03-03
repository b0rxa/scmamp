#' @title Friedman's post hoc raw p-values
#'
#' @description This function computes the raw p-values for the post hoc based on Friedman's test
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples
#' @return A matrix with all the pairwise raw p-values.
#' @details The test has been implemented according to the version in Demsar (2006), page 12
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @examples
#' data(data.garcia.herrera)
#' friedman.post(data.garcia.herrera)

friedman.post <- function (data , ...){
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  meanrank <- colMeans(rank.matrix(data))
  
  ## Generate all the pairs to test
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  ## Compute the p-value
  sd <- sqrt((k*(k+1))/(6*N))
  f<-function(x) (1 - pnorm(abs(meanrank[x[1]] - meanrank[x[2]])/sd))*2 ## Two tailed ...
  pvalues <- apply(pairs , MARGIN = 1 , FUN = f)
  matrix.raw <- matrix(rep(NA , k^2) , k)
  matrix.raw[pairs] <- pvalues
  matrix.raw[pairs[,c(2,1)]] <- pvalues
  colnames(matrix.raw) <- rownames(matrix.raw) <- colnames(data)
  matrix.raw
}

#' @title Function to use custom test to get the pairwise raw p-values
#'
#' @description This function computes the raw p-values using a custom function
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples
#' @param test Function to perform the test. It requires two parameters, \code{x} and \code{y}, the two samples to be compared, and it has to return the associeted p-value.
#' @return A matrix with all the pairwise raw p-values.
#' @examples
#' data(data.garcia.herrera)
#' test <- function(x , y , ...) t.test(x,y,paired=T)$p.value
#' custom.post(data.garcia.herrera , test)

custom.post <- function(data , test , ...){
  k <- dim(data)[2]
  N <- dim(data)[1]
  
  meanrank <- colMeans(rank.matrix(data))
  
  ## Generate all the pairs to test
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  ## Compute the p-value
  sd <- sqrt((k*(k+1))/(6*N))
  f<-function(x) test(x = data[ , x[1]] , y = data[ , x[2]])
  pvalues <- apply(pairs , MARGIN = 1 , FUN = f)
  matrix.raw <- matrix(rep(NA , k^2) , k)
  matrix.raw[pairs] <- pvalues
  matrix.raw[pairs[,c(2,1)]] <- pvalues
  colnames(matrix.raw) <- rownames(matrix.raw) <- colnames(data)
  matrix.raw
}

#' @title Tukey post hoc test for ANOVA
#'
#' @description This function computes all the pairwise p-values corrected using Tukey post hoc test
#' @param data Data set (matrix or data.frame) to apply the test. The column names are taken as the groups and the values in the matrix are the samples
#' @return A matrix with all the pairwise corrected p-values.
#' @details The test has been implemented according to Test 22 in Kanji (2006).
#' @references Kanji, G. K. (2006) \emph{100 Statistical Tests}. SAGE Publications Ltd, 3rd edition.
#' @examples
#' data(data.garcia.herrera)
#' anova.post(data.garcia.herrera)

anova.post <- function (data , ...){
  
  ## Implemented as Test 28 in 100 statistical tests
  k <- dim(data)[2]
  N <- dim(data)[1]
  nu <- k*(N-1)
  
  ## Sample variances
  var.vector <- apply(data , MARGIN = 2 , FUN = var)
  m.vector <- colMeans(data)
  
  ## Total variance. Given that all the samples have the same size, it is the average variance
  var <- mean(var.vector)
    
  ## Generate all the pairs to test
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  ## Compute the p-value
  f<-function(x) {
    d <- abs(m.vector[x[1]] - m.vector[x[2]])
    q <- d*sqrt(N/var)
    1 - ptukey(q,k,nu)
  }
  pvalues <- apply(pairs , MARGIN = 1 , FUN = f)
  matrix.raw <- matrix(rep(NA , k^2) , k)
  matrix.raw[pairs] <- pvalues
  matrix.raw[pairs[,c(2,1)]] <- pvalues
  colnames(matrix.raw) <- rownames(matrix.raw) <- colnames(data)
  matrix.raw
}

## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title Enforce the monotocity of a sequence
#
# @description This function ensures that the i-th element is not smaller than any previous one
# @param pvalues Vector where the correction has to be performed
# @return A vector with the corrected values
# @details In some corrections it can happen that after the correction the order of the resulting p-values is not te same as in the raw pvalues. A very simple example is the Holm correction where the ordered p-values are multiplied by m, m-1, ... , 2, 1. In case the last two p-values are, for instance, 0.7 and 0.8, the corrected p-values would be ... , 0.7*2, 0.8*1, i.e., ... , 1 , 0.8. This situation needs a correction to set the last value at 1.
# @examples
# pv <- c(4.4870e-07,1.0416e-06,0.00116,0.0288,0.032, 0.0303 , 0.0384 , 0.0384 , 1 , 1)
# correct.for.monotocity(pv)
correct.for.monotocity <- function (pvalues){
  sapply(1:length(pvalues) , function(x) max(pvalues[1:x]))
}


## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title Maximum number of true hypothesis
#
# @description This function gets the S(k) set as described in Shaffer (1985)
# @param k Number of algorithms
# @return List of maximum number of true hypothesis in a pair-wise comparsion of \code{k} classifiers
# @examples
# recursive.count(5)
# 
recursive.count <- function (k){
  res <- c(0)
  if (k>1){
    res <- c(res , recursive.count(k-1))
    for (j in 2:k){
      res <- c(res , recursive.count(k-j) + (factorial(j) / (2*factorial(j-2))))
    }
  }
  sort(unique(res))
}

#' @title Shaffer's correction of p-values in pair-wise comparisons
#'
#' @description This function implements the Shaffer's multiple testing correction when the p-values correspond with pair-wise comparisons
#' @param raw.matrix A matrix with the pair-wise p-values. The p-values have to be, at least, in the upper part of the matrix.
#' @return A symetric matrix with the corrected p-values.
#' @details The test has been implemented according to the version in Garcia and Herrera (2008), page 2680.
#' @references Garcia S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for All Pairwise Comparisons. \emph{Journal of Machine Learning Research}, 9, 2677-2694.
#' 
#' @examples
#' data(data.garcia.herrera)
#' raw.pvalues <- friedman.post(data.garcia.herrera)
#' shaffer.static (raw.pvalues)

shaffer.static <- function (raw.matrix){
  k <- dim(raw.matrix)[1]
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  
  raw.pvalues <- raw.matrix[pairs]
  sk <- recursive.count(k)[-1]
  t_i <- c(rep(sk[-length(sk)] , diff(sk)) , sk[length(sk)])
  t_i <- rev(t_i)
  
  o <- order(raw.pvalues)
  ## Order the p-values to apply the correction
  adj.pvalues <- raw.pvalues[o] * t_i  
  adj.pvalues <- sapply(adj.pvalues, function(x){min(x,1)})
  adj.pvalues <- correct.for.monotocity(adj.pvalues)
  
  adj.pvalues <- adj.pvalues[order(o)]
  adj.matrix <- raw.matrix
  adj.matrix[pairs] <- adj.pvalues
  adj.matrix[pairs[,2:1]] <- adj.pvalues  
  
  adj.matrix
}

## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title All the ordered subdivisions of a set
#
# @description This function is the base to crate all the partitions needed in the algorithm to create all the exhaustive sets in Figure 1, Garcia and Herrera (2008).
# @param set Set to be subdivided
# @return All the possible subdbisions of the set, including those where there are empty sets, but without repetitions. The format is a list of lists. Each element in the main list is a list containing two vectors, \code{s1} and \code{s2}, the two subsets of the set passed in the argument
 
subdivisions <- function (set){
  if (length(set)==1){ ## trivial case, only one posibility
    res <- list(list(s1=set, s2 = vector()))
  }else{## In the general case, subdivide the set without the last element and then add ...
    n <- length(set)
    last <- set[n]
    sb <- subdivisions (set[-n])
    ## ... the last element in all the first sets ...
    res <- lapply(sb , FUN = function (x) {list(s1=c(x$s1 , last) , s2=x$s2)})
    ## ... the last element in all the second sets ...
    res <- c(res , lapply(sb , FUN = function (x) {list(s1=x$s1 , s2=c(x$s2 , last))}))
    ## ... and finally the subdivision that has the last element alone in s2
    res <- c(res , list(list(s1=last , s2=set[-n])))
  }
  res
}

## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title Complete set of (non-empty) partitions of a set
#
# @description This function creates all the possible partitions indicated in the 6th step of the algorithm shown in Figure 1, Garcia and Herrera (2008).
# @param set Set to be partitioned
# @return All the possible partitions of the set where the last element is in the second subset and the first subset is not empty.

partition <- function (set){
  n <- length(set)
  if (n == 1){ ##In the trivial case it returns the only possible subdivision. This case violates the idea of having the last element in the last set, but it is the only exception and it does not pose a problema as the repetition of the set is under control.
    res <- subdivisions(set)
  }else{ ## In the general case, get all the subdivisions (which cannot have empty sets in the first subset) and add the last element at the end of each second subset.
    last <- set[n]
    sb <- subdivisions(set[-n])
    ## Add the last element only to the second subset (see Figur 1 in Garcia and Herrera, 2008)
    res <- lapply(sb , FUN = function (x) {list(s1=x$s1 , s2=c(x$s2 , last))})
  }
  res
}

## DO NOT USE ROXYGENIZE, IT IS NOT A PUBLIC FUNCTION!
# @title Remove repetitions in a list of subsets
#
# @description This function removes the repetitions in the result returned by the function \code{\link{exhaustive.sets}}
# @param E list of exhaustive sets
# @return The list without repetitions

unique.exhaustive.sets <- function (E){
  E.new <- list(E[[1]])
  
  already.in <- function (e){
    compare.with.e <- function (en){
      if (all(dim(en)==dim(e))) {
        return (all(en==e))
      }else{
        return (FALSE)
      }
    }    
    comparisons <- unlist(lapply (E.new , FUN = compare.with.e))
    any(comparisons)
  }
  ## Sequentially add partitions if they are not already in the solution
  bk<-lapply (E , FUN = function (x) if (!already.in (x)) E.new <<- c(E.new , list(x)))
  E.new
}


#' @title Create the complet set of exhaustive sets
#'
#' @description This function implements the algorithm in Figure 1, Garcia and Herrera (2008) to create, given a set, the complete set of exhaustive sets E
#' @param set Set to create the exhaustive sets. The complexity of this algorithm is huge, so use with caution for sets of more than 7-8 elements.
#' @return A list with all the possible exhaustive sets, without repetitions
#' @examples
#' exhaustive.sets(c("A","B","C","D"))
exhaustive.sets <- function (set){
  k <- length(set)
  if (k<=length(E)){ ## Reuse the computed sets stored in the variable E
    E_l <- lapply(E[[k]] , FUN=function(x) matrix(set[x],nrow=2))
  }else{
    ## All possible pairwise comparisons, Garcia and Herrera, Table 1 steps 1-5
    if (k==1) return(NULL) ## No sense the set with one classifier
    if (k==2){ ## Base case, no lists
      E_l <- list(matrix(c(set[1],set[2]),ncol=1))
    }else{
      pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
      E_l <- list(apply (pairs , MARGIN = 1 , FUN = function(x) c(set[x[1]] , set[x[2]])))
    }
    ## Main loop
    if (length(set)>2){
      partitions <- partition(set)  ## Sets in step 6
      process.partition <- function (x){ ## Function to perform the main loop
        E1 <- exhaustive.sets (x$s1)
        E2 <- exhaustive.sets (x$s2)
        if (!is.null(E1)){
          E_l <<- c(E_l , E1)
          if (!is.null(E2)){
            lapply (E1 , FUN = function(e1){
              E_l <<- c(E_l , lapply(E2 , FUN = function(e2) cbind(e1,e2)))
            })
          }
        }      
        if (!is.null(E2)) E_l <<- c(E_l , E2)
      }
      ## Do the loop (bk is just to avoid printing trash in the screen ...)
      bk <- lapply(partitions , FUN = process.partition)
    }
    E_l <- unique.exhaustive.sets(E_l)
  }
  gc()
}

#' @title Bergmann and Hommel dynamic correction of p-values
#'
#' @description This function takes the particular list of possible hypthesis to correct for multiple testing, as defined in Bergmann and Hommel (1994)
#' @param raw.matrix Raw p-values in a matrix
#' @return A matrix with the corrected p-values
#' @details The test has been implemented according to the version in Garcia and Herrera (2008), page 2680-2682.
#' @references Garcia S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for All Pairwise Comparisons. \emph{Journal of Machine Learning Research}, 9, 2677-2694.
#' @examples
#' data(data.garcia.herrera)
#' raw.pvalues <- friedman.post(data.garcia.herrera)
#' bergmann.hommel.dynamic (raw.pvalues)
bergmann.hommel.dynamic <- function (raw.matrix){
  ## Load the exhaustive sets
  data("exhaustive.sets")
  k <- dim(raw.matrix)[1]
  if(k>8) stop ("Sorry, this method is only available for 8 or less algorithms.")
  pairs <- do.call(rbind,sapply(1:(k-1), FUN=function(x) cbind((x),(x+1):k)))
  Ek <- E[[k]]
  raw.pvalues <- raw.matrix[pairs]
  get.corrected <- function(i){
    aux <- lapply(Ek , FUN = function(I){
      ## check that i is in I
      if (any(colSums(pairs[i,] == I)==2)){
        nu = dim(I)[2] * min(raw.matrix[t(I)])
        min(nu,1)
      }
    })
    max(unlist(aux))
  }
  adj.pvalues <- sapply(1:dim(pairs)[1] , FUN = get.corrected)
  ## Correct any possible inversion
  o <- order(raw.pvalues)
  aux <- correct.for.monotocity(adj.pvalues[o])
  adj.pvalues <- aux[order(o)]
    
  adj.matrix <- raw.matrix
  adj.matrix[pairs] <- adj.pvalues
  adj.matrix[pairs[,2:1]] <- adj.pvalues  
  
  adj.matrix
}