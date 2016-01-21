# NORMALITY CHECK PLOTS ---------------------------------------------------

#' @title Gaussian distribution quantile-quantile plot
#'
#' @description This function creates a quantile-quantile plot to assess the goodness of fit of a Gaussian distribution to a given sample.
#' @param data List of data points to check
#' @param ... The plot is created using \code{ggplot2}. This special parameter can be used to pass additional parameters to the \code{\link{geom_point}} function used to plot the sample points.
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{\link{plotDensities}}
#' @examples
#' ## Skewed distribution
#' sample <- rbeta(100 , 2 , 50)
#' qqplotGaussian(sample)
#' ## Symmetric distribution
#' sample <- rbeta(100 , 5 , 5)
#' qqplotGaussian(sample)

qqplotGaussian <- function (data, ...) {
#   if(!require(ggplot2)) {
#     stop("This function requires the package ggplot2, which is not ",
#          "installed. You can install it typing install.packages('ggplot2')")
#   }
  
  processVector <- function (sample) {
    # Auxiliar function to get the points for a single qqplot
    sample <- sort(sample)
    m      <- mean(sample)
    sd     <- sd(sample)
    n      <- length(sample)
    
    q.gaussian <- qnorm(seq(from=(1 / n), to=((n - 1) / n), by=(1/n)), 
                        mean=m, sd=sd)
    q.sample <- sample[-n]
    
    df <- data.frame(Empirical=q.sample, Gaussian=q.gaussian)  
    return (df)
  }
  
  if (is.vector(data)) {
    df <- processVector(data)
    facet <- NULL
  } else {
    aux <- lapply(1:ncol(data), 
                  FUN=function(i) {
                    alg.pts  <- processVector(data[, i])
                    alg.name <- names(data)[i]
                    return(cbind(alg.pts, Algorithm=alg.name))
                  })
    df <- do.call(rbind, aux)
    facet <- facet_wrap(~Algorithm, scales="free")
  }
  
  gplot <- ggplot(df, aes_string(x="Empirical", y="Gaussian")) +
           geom_abline(slope=1, intercept=0, col="darkgray", size=1.1) +  
           geom_point(...) + facet
  
  return(gplot)
}


#' @title Kernel based density estimation of the samples
#'
#' @description This function estimates and plots the densities of the results of each algorithm
#' @param data A matrix where columns represent the algorithms
#' @param ... The plot is created using \code{ggplot2}. This special parameter can be used to pass additional parameters to the \code{\link{geom_line}} function used to plot the sample points. It can also be used to pass additional arguments to the \code{density} function, which is used to eastimate the densities.
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{\link{qqplotGaussian}}
#' @examples
#' data(data_gh_2010)
#' plotDensities(data.gh.2010)
#' 

plotDensities <- function (data, ...) {
#   if(!require(ggplot2)) {
#     stop("This function requires the package ggplot2, which is not installed. ",
#          "You can install it typing install.packages('ggplot2')")
#   }
  
  if (is.vector(data)){
    d  <- density(data)
    df <- data.frame(Value=d$x, Density=d$y)
    mapping <- aes_string(x="Value", y="Density")
  } else {
    k <- dim(data)[2]
    aux <- lapply(1:k, 
                  FUN=function (x) {
                    d  <- density(data[, x], ...)
                    df <- data.frame(Algorithm=colnames(data)[x], 
                                     Value=d$x, Density=d$y)
                    return(df)
    })
    df <- do.call(rbind, aux)
    mapping <- aes_string(x="Value", y="Density", col="Algorithm")
  }
  gplot <- ggplot(df, mapping) + geom_line(...)
  return(gplot)
}

# PLOTTING THE RESULTS ----------------------------------------------------

#' @title Plotting the p-value matrix
#'
#' @description This function plots the p-value matrix as a tile plot.
#' @param pvalue.matrix Matrix with the p-values to plot
#' @param alg.order A permutation indicating the reordering of the algorithms
#' @param show.pvalue Logical value indicating whether the numerical values have to be printed
#' @param font.size Size of the p-values, if printed
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{\link{drawAlgorithmGraph}}, \code{\link{plotCD}}
#' @examples
#' data(data_gh_2008)
#' pvalues <- friedmanPost(data.gh.2008)
#' ordering <- order(summarizeData(data.gh.2008))
#' plotPvalues(pvalues, alg.order=ordering)
#' 
 
plotPvalues <- function(pvalue.matrix, alg.order=NULL, show.pvalue=TRUE, font.size=5) {
#   if(!require(reshape2)) {
#     stop("This function requires the package reshape2, which is not installed. ",
#          "You can install it typing install.packages('reshape2')")
#   }
  
#   if(!require(ggplot2)) {
#     stop("This function requires the package ggplot2, which is not installed. ",
#          "You can install it typing install.packages('ggplot2')")
#   }
  
  # Convert the matrix into a data frame and order the algorithms according to 
  # the desired order.
  df <- melt(pvalue.matrix)
  colnames(df) <- c("X", "Y", "p.value")
  if (!is.null(alg.order)) {
    l <- colnames(pvalue.matrix)[alg.order]
    df$X <- factor(df$X, levels=l)
    df$Y <- factor(df$Y, levels=l)
  }
  
  gplot <- ggplot(df, aes_string(x="X", y="Y", fill="p.value")) + geom_tile(col="white") +
          scale_fill_continuous("p-value") + labs(x="Algorithm" , y="Algorithm")
  
  if (show.pvalue) {
    p.value <- df$p.value
    gplot <- gplot + geom_text(aes_string(label="round(p.value, 2)"), 
                               size=font.size, col="white")
  }
  return(gplot)
} 


#' @title Critical difference plot
#'
#' @description This function plots the critical difference plots shown in Demsar (2006)
#' @param results.matrix Matrix or data frame with the results for each algorithm
#' @param alpha Significance level to get the critical difference. By default this value is 0.05
#' @param cex Numeric value to control the size of the font. By default it is set at 0.75.
#' @param ... Additional arguments for \code{\link{rankMatrix}}
#' @seealso \code{\link{drawAlgorithmGraph}},   \code{\link{plotRanking}}, \code{\link{plotPvalues}}
#' @references Demsar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. \emph{Journal of Machine Learning Research}, 7, 1-30.
#' @examples
#' data(data_gh_2008)
#' plotCD(data.gh.2008, alpha=0.01)
#' 
plotCD <- function (results.matrix, alpha=0.05, cex=0.75, ...) {
  k <- dim(results.matrix)[2]
  N <- dim(results.matrix)[1]
  cd <- getNemenyiCD(alpha=alpha, num.alg=k, num.problems=N)
  
  mean.rank <- sort(colMeans(rankMatrix(results.matrix, ...)))
  
  # Separate the algorithms in left and right parts
  lp <- round(k/2)
  left.algs <- mean.rank[1:lp]
  right.algs <- mean.rank[(lp+1):k]  
  max.rows <- ceiling(k/2)
  
  # Basic dimensions and definitions
  char.size    <- 0.001  # Character size
  line.spacing <- 0.25   # Line spacing for the algorithm name
  m            <- floor(min(mean.rank))
  M            <- ceiling(max(mean.rank))
  max.char     <- max(sapply(colnames(results.matrix), FUN = nchar))  # Longest length of a label
  text.width   <- (max.char + 4) * char.size
  w            <- (M-m) + 2 * text.width
  h.up         <- 2.5 * line.spacing  # The upper part is fixed. Extra space is for the CD
  h.down       <- (max.rows + 2.25) * line.spacing # The lower part depends on the no. of algorithms. 
  # The 2 extra spaces are for the lines that join algorithms
  tick.h       <- 0.25 * line.spacing
  
  label.displacement <- 0.1    # Displacement of the label with respect to the axis
  line.displacement  <- 0.025  # Displacement for the lines that join algorithms
  
  # Background of the plot
  plot(0, 0, type="n", xlim=c(m - w / (M - m), M + w / (M - m)), 
       ylim=c(-h.down, h.up), xaxt="n", yaxt="n", xlab= "", ylab="", bty="n")
  
  # Draw the axis
  lines (c(m,M), c(0,0))
  dk <- sapply(m:M, 
               FUN=function(x) {
                 lines(c(x,x), c(0, tick.h))
                 text(x, 3*tick.h, labels=x, cex=cex)
               })
  
  # Draw the critical difference
  lines(c(m, m + cd), c(1.75 * line.spacing, 1.75 * line.spacing))
  text(m + cd / 2, 2.25 * line.spacing, "CD", cex=cex)
  lines(c(m, m), c(1.75 * line.spacing - tick.h / 4, 
                   1.75 * line.spacing + tick.h / 4))
  lines(c(m + cd, m + cd), c(1.75 * line.spacing - tick.h / 4, 
                             1.75 * line.spacing + tick.h / 4))
  
  # Left part, labels
  dk <- sapply (1:length(left.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=m - label.displacement, y=line.h, 
                       labels=names(left.algs)[x], cex=cex, adj=1)
                  lines(c(m, left.algs[x]), c(line.h, line.h))
                  lines(c(left.algs[x], left.algs[x]), c(line.h, 0))
                })
  
  # Right part, labels
  dk <- sapply (1:length(right.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=M + label.displacement, y=line.h, 
                       labels=names(right.algs)[x], cex=cex, adj=0)
                  lines(c(M, right.algs[x]), c(line.h, line.h))
                  lines(c(right.algs[x], right.algs[x]), c(line.h, 0))
                })
  
  # Draw the lines to join algorithms
  getInterval <- function (x) {
    from <- mean.rank[x]
    diff <- mean.rank - from
    ls <- which(diff > 0 & diff < cd)
    if (length(ls) > 0) {
      c(from, mean.rank[max(ls)])
    }
  }
  
  intervals <- mapply (1:k, FUN=getInterval)
  aux <- do.call(rbind, intervals)
  
  # With this strategy, there can be intervals included into bigger ones
  # We remove them in a sequential way
  to.join <- aux[1,]
  if(nrow(aux) > 1) {  
    for (r in 2:nrow(aux)) {
      if (aux[r - 1, 2] < aux[r, 2]) {
        to.join <- rbind(to.join, aux[r, ])
      }
    }
  }
  
  row <- c(1)
  # Determine each line in which row will be displayed
  if (!is.matrix(to.join)) {  # To avoid treating vector separately
    to.join <- t(as.matrix(to.join))
  }
  nlines <- dim(to.join)[1]
  
  for(r in 1:nlines) {
    id <- which(to.join[r, 1] > to.join[, 2])
    if(length(id) == 0) {
      row <- c(row, tail(row, 1) + 1)
    } else {
      row <- c(row, min(row[id]))
    }
  }
  
  step <- max(row) / 2
  
  # Draw the line
  dk <- sapply (1:nlines, 
                FUN = function(x) {
                  y <- -line.spacing * (0.5 + row[x] / step)
                  lines(c(to.join[x, 1] - line.displacement, 
                          to.join[x, 2] + line.displacement), 
                        c(y, y), lwd=3)
                })
}



#' @title Ranking Plots
#'
#' @description This function creates a plot similar to the critical difference plot, but applicable to any corrected pvalue.
#' @param pvalues Matrix or data frame with the p-values used to determine the differences
#' @param summary Summary values used to place the algorithms. Typically it will be the average ranking, but it can be any other value 
#' @param alpha Significance level to determine the significativity of the differences. By default this value is 0.05
#' @param cex Numeric value to control the size of the font. By default it is set at 0.75.
#' @param decreasing A logical value to determine whether the values have to be plotter from smaller to larger or the other way round.
#' @seealso \code{\link{drawAlgorithmGraph}}, \code{\link{plotCD}}, \code{\link{plotPvalues}}
#' @examples
#' data(data_gh_2008)
#' test <- postHocTest(data.gh.2008, test="friedman", correct="bergmann", use.rank=TRUE)
#' plotRanking(pvalues=test$corrected.pval, summary=test$summary, alpha=0.05)
#' 
plotRanking <- function (pvalues, summary, alpha=0.05, cex=0.75, decreasing=FALSE) {
  k <- length(summary)

  if (is.matrix(summary)) {
    if (ncol(summary) == 1) {
      summary <- summary[, 1]    
    } else {
      summary <- summary[1, ]
    }
  }
  
  # Check the names
  if (!all(sort(colnames(pvalues)) %in% sort(names(summary))) |
      !all(sort(names(summary)) %in% sort(colnames(pvalues)))) {
    stop("The column names of 'pvalues' and the names of 'summary' have to contain the same names")
  }
  
  # Reorder the pvalues and mean ranks
  o <- order(summary, decreasing=decreasing)
  summary <- summary[o]
  pvalues <- pvalues[names(summary),names(summary)]
  
  
  # Separate the algorithms in left and right parts
  lp <- round(k/2)
  left.algs <- summary[1:lp]
  right.algs <- summary[(lp+1):k]  
  max.rows <- ceiling(k/2)
  
  # Basic dimensions and definitions
  char.size    <- 0.001  # Character size
  line.spacing <- 0.25   # Line spacing for the algorithm name
  m            <- floor(min(summary))
  M            <- ceiling(max(summary))
  max.char     <- max(sapply(colnames(pvalues), FUN = nchar))  # Longest length of a label
  text.width   <- (max.char + 4) * char.size
  w            <- (M-m) + 2 * text.width
  h.up         <- 2.5 * line.spacing  # The upper part is fixed. Extra space is for the CD
  h.down       <- (max.rows + 2.25) * line.spacing # The lower part depends on the no. of algorithms. 
  # The 2 extra spaces are for the lines that join algorithms
  tick.h       <- 0.25 * line.spacing
  
  label.displacement <- 0.1    # Displacement of the label with respect to the axis
  line.displacement  <- 0.025  # Displacement for the lines that join algorithms
  
  # Background of the plot
  plot(0, 0, type="n", xlim=c(m - w / (M - m), M + w / (M - m)), 
       ylim=c(-h.down, h.up), xaxt="n", yaxt="n", xlab= "", ylab="", bty="n")
  
  # Draw the axis
  lines (c(m,M), c(0,0))
  dk <- sapply(m:M, 
               FUN=function(x) {
                 lines(c(x,x), c(0, tick.h))
                 text(x, 3*tick.h, labels=x, cex=cex)
               })
  
  
  # Left part, labels
  dk <- sapply (1:length(left.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=m - label.displacement, y=line.h, 
                       labels=names(left.algs)[x], cex=cex, adj=1)
                  lines(c(m, left.algs[x]), c(line.h, line.h))
                  lines(c(left.algs[x], left.algs[x]), c(line.h, 0))
                })
  
  # Right part, labels
  dk <- sapply (1:length(right.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=M + label.displacement, y=line.h, 
                       labels=names(right.algs)[x], cex=cex, adj=0)
                  lines(c(M, right.algs[x]), c(line.h, line.h))
                  lines(c(right.algs[x], right.algs[x]), c(line.h, 0))
                })
  
  # Draw the lines to join algorithms
  getInterval <- function (x) {
    ls <- which(pvalues[x, ] > alpha)
    if (length(ls) > 0) {
      c(summary[x], summary[max(ls)])
    }
  }
  
  intervals <- mapply (1:k, FUN=getInterval)
  
  # Under some circumstances the function does not return a matrix ...
  if (is.matrix(intervals)) {
    aux <- t(intervals)
  } else {
    aux <- do.call(rbind, intervals)
  }
  
  # With this strategy, there can be intervals included into bigger ones
  # We remove them in a sequential way
  to.join <- aux[1,]
  if(nrow(aux) > 1) {  
    for (r in 2:nrow(aux)) {
      if (aux[r - 1, 2] < aux[r, 2]) {
        to.join <- rbind(to.join, aux[r, ])
      }
    }
  }
  
  row <- c(1)
  # Determine each line in which row will be displayed
  if (!is.matrix(to.join)) {  # To avoid treating vector separately
    to.join <- t(as.matrix(to.join))
  }
  nlines <- dim(to.join)[1]
  
  for(r in 1:nlines) {
    id <- which(to.join[r, 1] > to.join[, 2])
    if(length(id) == 0) {
      row <- c(row, tail(row, 1) + 1)
    } else {
      row <- c(row, min(row[id]))
    }
  }
  
  step <- max(row) / 2
  
  # Draw the line
  dk <- sapply (1:nlines, 
                FUN = function(x) {
                  y <- -line.spacing * (0.5 + row[x] / step)
                  lines(c(to.join[x, 1] - line.displacement, 
                          to.join[x, 2] + line.displacement), 
                        c(y, y), lwd=3)
                })
}



#' @title Hypotheses represented as a graph
#'
#' @description This function can be used to plot a graph where algorithms are nodes and  algorithms that cannot be regarded as different are joined by an edge.
#' @param pvalue.matrix Matrix with the p-values
#' @param mean.value Vector of values to be written together with the name of the algorithm
#' @param ... Additional parameters to the Rgraphviz function. This is mainly to change the layout of the graph
#' @param alpha Significance level to determine which hypotheses are rejected.
#' @param font.size Size of the font for the node labels.
#' @param highlight A character indicating which node has to be highlighted. It can be the one with the maximum value (\code{'max'}), the minimum value (\code{'min'}) or none (\code{'none'}).
#' @param highlight.color Any R valid color for the highlighted node.
#' @param node.color Any R valid color for the non-highlighted nodes.
#' @param font.color Any R valid color for the node labels.
#' @param digits Number of digits to display the value associated to each node
#' @param node.width Numeric value for the width of the node
#' @param node.height Numeric value for the height of the node
#' @seealso \code{\link{plotPvalues}}, \code{\link{plotRanking}}, \code{\link{plotCD}}
#' @examples
#' data(data_blum_2015)
#' data <- filterData(data.blum.2015, condition="Size == 1000", remove.cols=1:2)
#' res <- postHocTest(data, test = "friedman", use.rank=TRUE, correct="bergmann")
#' drawAlgorithmGraph(res$corrected.pval, res$summary)
#' 

drawAlgorithmGraph <- function (pvalue.matrix, mean.value, ..., 
                                alpha=0.05, font.size=15, highlight="min", 
                                highlight.color="chartreuse3", node.color="gray30", 
                                font.color="white", digits=2, 
                                node.width=5, node.height=2) {
#   if(!require(Rgraphviz)) {
#     stop("This function requires the package Rgraphviz, which is not installed. ",
#          "You can install it typing\n\n ",
#          "source('http://www.bioconductor.org/biocLite.R')\n\n and then\n\n ",
#          "biocLite('Rgraphviz')\n\n")
#   }
  # Just in case we have a matrix ...
  if (is.matrix(mean.value)) {
    mean.value <- mean.value[1, ]
  }
  
  if (!all(colnames(pvalue.matrix) %in% names(mean.value))) {
    stop ("The names of the algorithms in the matrix and the mean.value vector ",
          "do not match")
  }
  
  hypothesis.matrix <- pvalue.matrix > alpha
  
  nc <- rep(node.color, length(mean.value))
  if(highlight == "min") {
    nc[which.min(mean.value)] <- highlight.color
  } else if(highlight == "max") {
    nc[which.max(mean.value)] <- highlight.color
  }
    
  hypothesis.matrix[is.na(hypothesis.matrix)]<-FALSE
  adj.matrix <- hypothesis.matrix
  colnames(adj.matrix) <- names(mean.value)
  rownames(adj.matrix) <- names(mean.value)
  
  nl <- paste(names(mean.value), "\\\n", round(mean.value, digits), sep="")
  
  am.graph <- new("graphAM", adjMat=adj.matrix, edgemode="undirected")
  names(nc) <- graph::nodes(am.graph)
  names(nl) <- graph::nodes(am.graph)
  
  nAttrs <- list()
  nAttrs$label      <- nl
  nAttrs$fillcolor  <- nc
  
  attrs <- list(node=list(shape="rectangle", width=node.width, height=node.height,
                            fontcolor=font.color, fontsize=font.size))
  plot(am.graph, ... , nodeAttrs=nAttrs, attrs=attrs)
}