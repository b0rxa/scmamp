# NORMALITY CHECK PLOTS ---------------------------------------------------

#' @title Gaussian distribution quantile-quantile plot
#'
#' @description This function creates a quantile-quantile plot to assess the goodness of fit of a Gaussian distribution to a given sample.
#' @param data List of data points to check
#' @param ... The plot is created using \code{ggplot2}. This special parameter can be used to pass additional parameters to the \code{\link{geom_point}} function used to plot the sample points.
#' @return A \code{\linkS4class{ggplot}} object.
#' @seealso \code{\link{plotDensities}}
#' @examples
#' ## Skewed distribution
#' sample <- rbeta(100 , 2 , 50)
#' qqplotGaussian(sample)
#' ## Symmetric distribution
#' sample <- rbeta(100 , 5 , 5)
#' qqplotGaussian(sample)

qqplotGaussian <- function (data, ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  
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
    facet <- ggplot2::facet_wrap(~Algorithm, scales="free")
  }
  
  gplot <- ggplot2::ggplot(df, ggplot2::aes(x=Empirical, y=Gaussian)) +
    ggplot2::geom_abline(slope=1, intercept=0, col="darkgray", size=1.1) +  
    ggplot2::geom_point(...) + facet
  
  return(gplot)
}


#' @title Kernel based density estimation of the samples
#'
#' @description This function estimates and plots the densities of the results of each algorithm
#' @param data A matrix where columns represent the algorithms
#' @param ... The plot is created using \code{ggplot2}. This special parameter can be used to pass additional parameters to the \code{\link{geom_line}} function used to plot the sample points. It can also be used to pass additional arguments to the \code{density} function, which is used to eastimate the densities.
#' @return A \code{\linkS4class{ggplot}} object.
#' @seealso \code{\link{qqplotGaussian}}
#' @examples
#' data(data_gh_2010)
#' plotDensities(data.gh.2010)
#' 

plotDensities <- function (data, ...) {
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  
  if (is.vector(data)){
    d  <- density(data)
    df <- data.frame(Value=d$x, Density=d$y)
    mapping <- ggplot2::aes(x=Value, y=Density)
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
    mapping <- ggplot2::aes(x=Value, y=Density, col=Algorithm)
  }
  gplot <- ggplot2::ggplot(df, mapping) + ggplot2::geom_line(...)
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
#' @return A \code{\linkS4class{ggplot}} object.
#' @seealso \code{\link{drawAlgorithmGraph}}, \code{\link{plotCD}}
#' @examples
#' data(data_gh_2008)
#' pvalues <- friedmanPost(data.gh.2008)
#' ordering <- order(summarizeData(data.gh.2008))
#' plotPvalues(pvalues, alg.order=ordering)
#' 

plotPvalues <- function(pvalue.matrix, alg.order=NULL, show.pvalue=TRUE, font.size=5) {
  
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  
  # Convert the matrix into a data frame and order the algorithms according to 
  # the desired order.
  df <- melt(pvalue.matrix)
  colnames(df) <- c("X", "Y", "p.value")
  if (!is.null(alg.order)) {
    l <- colnames(pvalue.matrix)[alg.order]
    df$X <- factor(df$X, levels=l)
    df$Y <- factor(df$Y, levels=l)
  }
  
  gplot <- ggplot2::ggplot(df, ggplot2::aes(x=X, y=Y, fill=p.value)) + ggplot2::geom_tile(col="white") +
    ggplot2::scale_fill_continuous("p-value") + ggplot2::labs(x="Algorithm" , y="Algorithm")
  
  if (show.pvalue) {
    p.value <- df$p.value
    gplot <- gplot + ggplot2::geom_text(ggplot2::aes(label=round(p.value, 2)),
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
  
  opar <- par(mai = c(0,0,0,0))
  on.exit(par(opar))
  
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
  
  label.displacement <- 0.25    # Displacement of the label with respect to the axis
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
                  lines(c(m - label.displacement*0.75, left.algs[x]), c(line.h, line.h))
                  lines(c(left.algs[x], left.algs[x]), c(line.h, 0))
                })
  
  # Right part, labels
  dk <- sapply (1:length(right.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  text(x=M + label.displacement, y=line.h, 
                       labels=names(right.algs)[x], cex=cex, adj=0)
                  lines(c(M + label.displacement*0.75, right.algs[x]), c(line.h, line.h))
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
  if(NROW(aux) > 0) {
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
  
  opar <- par(mai=c(0, 0, 0, 0), mgp=c(0, 0, 0))
  on.exit(par(opar))
  
  k <- length(summary)
  
  # dirty patch to for decreasing=TRUE case
  if(decreasing)
  {
   summary <- k - summary + 1
   decreasing = FALSE
  }
  
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
  
  # In case we have a vector of p-values (all vs. control), identify the control
  control <- NULL
  if (nrow(pvalues)==1) {
    control <- colnames(pvalues)[is.na(pvalues)]
  }
  
  # Reorder the pvalues and mean ranks
  o <- order(summary, decreasing=decreasing)
  summary <- summary[o]
  if (nrow(pvalues)>1) {
    pvalues <- pvalues[names(summary), names(summary)]
  }else{
    pvalues <- matrix(pvalues[1,names(summary)], nrow=1)
    colnames(pvalues) <- names(summary)
  }
  
  
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
  max.char     <- max(sapply(colnames(pvalues), FUN=nchar))  # Longest length of a label
  text.width   <- (max.char + 4) * char.size
  w            <- (M-m) + 2 * text.width
  h.up         <- line.spacing  # The upper part is fixed.
  h.down       <- (max.rows + 2.25) * line.spacing # The lower part depends on the no. of algorithms. 
  # The 2 extra spaces are for the lines that join algorithms
  tick.h       <- 0.25 * line.spacing
  
  label.displacement <- 0.25   # Displacement of the label with respect to the axis
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
                  name <- names(left.algs)[x]
                  if (!is.null(control) && name==control) {
                    font = 4
                  }else{
                    font = 1
                  }
                  text(x=m - label.displacement, y=line.h, 
                       labels=name, cex=cex, adj=1, font=font)
                  lines(c(m - label.displacement*0.75, left.algs[x]), c(line.h, line.h))
                  lines(c(left.algs[x], left.algs[x]), c(line.h, 0))
                })
  
  # Right part, labels
  dk <- sapply (1:length(right.algs), 
                FUN=function(x) {
                  line.h <- -line.spacing * (x + 2)
                  name <- names(right.algs)[x]
                  if (!is.null(control) && name==control) {
                    font = 4
                  }else{
                    font = 1
                  }
                  text(x=M + label.displacement, y=line.h, 
                       labels=name, cex=cex, adj=0, font=font)
                  lines(c(M + label.displacement*0.75, right.algs[x]), c(line.h, line.h))
                  lines(c(right.algs[x], right.algs[x]), c(line.h, 0))
                })
  
  # Draw the lines to join algorithms
  if (nrow(pvalues)==1) {
    to.join <- summary[which(is.na(pvalues) |  pvalues > alpha)]
    if (length(to.join)>1) {
      lines (c(min(to.join), max(to.join)) , c(0,0), lwd=3)
    }
    
  }else{
    getInterval <- function (x) {
      ls <- which(pvalues[x, ] > alpha)
      # Only retail those to the right in the matrix
      ls <- ls[ls > x]
      res <- NULL
      if (length(ls) > 0) {
        res <- c(as.numeric(summary[x]), as.numeric(summary[max(ls)]))
      }
      return(res)
    }
    
    intervals <- mapply (1:(k-1), FUN=getInterval)
    
    # Under some circumstances the function does not return a matrix ...
    if (is.matrix(intervals)) {
      aux <- t(intervals)
    } else {
      aux <- do.call(rbind, intervals)
    }
    
    # First, chech that there are lines to draw!
    if (length(aux) > 0) {
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
  }
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
#' ## This function requieres the package Rgraphviz
#' # drawAlgorithmGraph(res$corrected.pval, res$summary)
#' 

drawAlgorithmGraph <- function (pvalue.matrix, mean.value, ..., 
                                alpha=0.05, font.size=15, highlight="min", 
                                highlight.color="chartreuse3", node.color="gray30", 
                                font.color="white", digits=2, 
                                node.width=5, node.height=2) {
  
  if (!requireNamespace("Rgraphviz", quietly=TRUE)) {
    stop("This function requires the Rgraphviz package. Please install it. Note that the packages is currently available at Bioconductor", call.=FALSE)
  }
  
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



#' @title Plotting the (marginal) posterior densities in Bayesian analyses
#'
#' @description This function plots, univariately, the posterior densities obtained
#' in Bayesian analyses. 
#' @param results A list containing, at least three elements, one named \code{approximate}, which is
#' a logical value indicating whether the posterior is a function or a sample, \code{rope}, a two dimensional
#' vector with the minimum and maximum values of the rope and \code{posterior}, either a one parameter function or
#' a matrix (or data.frame) where each row is a sample and each column a sampled parameter
#' @param parameter Either a string or a number indicating, in case the posterior is approximated, the parameter
#' to be ploted (i.e., the name or the index of a column in the sample matrix)
#' @param ... Additional parameters to the Rgraphviz function. This is mainly to change the layout of the graph
#' @param plot.rope  A logical value indicating whether the rope has to be plotted or not. Note that not for all
#' parameter the rope makes sense
#' @param num.points Number of points used to plot the functions
#' @param plot.samples A logical value. If true, the samples are plotted (only when the posterior is approximate)
#' @param alpha Numeric value for the transparency of the points, only applicable if \code{plot.samples} is true
#' @return An object of class \linkS4class{ggplot} with the plot
#' @details 
#' Note that if the methods are exact (not simulated), the true density can be plotted but,
#' for those cases where the posterior is approximated through sampling, the function will
#' plot a kernel density estimation of the posterior and, thus, the probabilities computed
#' by other functions are not directly the areas under the densities.
#' @examples
#' x <- rnorm(25, 1, 2)
#' y <- rnorm(25, 1.1, 2)
#' results <- bCorrelatedTtest(x=x, y=y, rho=0, rope=c(-0.05, 0.05))
#' plotPosterior(results)
#' 

plotPosterior <- function(results, parameter=1, num.points=1000, plot.rope=TRUE, plot.samples=TRUE, alpha=NULL, ...) {
  
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  
  if (results$approximate){
    if (is.character(parameter) & !any(names(results$posterior)==parameter)){
      stop("The parameter ", parameter, " is not contained in the sample of the posterior distribution")
    } else if(parameter <= 0 | parameter > ncol(results$posterior)){
      stop("The parameter argument has to be either a valid column name ", 
           "or a valid column index for the posterior sample matrix")
    }
    sample <- data.frame (Sample=results$posterior[, parameter])
    aux <- density(results$posterior[, parameter], ...)
    dens <- data.frame(Difference=aux$x, Posterior=aux$y)
    
    if (is.null(alpha)) {
      alpha <- min(1, 250/nrow(sample))
    }
    
    gplot <- ggplot(dens, aes(x=Difference, y=Posterior)) + geom_line() 
    
    if(plot.samples) {
      gplot <- gplot + geom_point(data=sample, aes(x=Sample, y=0), 
                                  position=position_jitter(height=0.05*max(dens$Posterior)), alpha=alpha) 
    }
  } else {
    qpos <- results$additional$qposterior
    dpos <- results$posterior
    rope <- results$parameters$rope
    
    # Get the data ready
    x <- seq(min(qpos(0.0005),rope[1]), max(qpos(0.9995), rope[2]), length.out=num.points)
    df <- data.frame(Difference=x, Posterior=dpos(x))
    
    gplot <- ggplot(df, aes(x=Difference, y=Posterior)) + geom_line() 
  }
  
  if (plot.rope) {
    gplot <- gplot + geom_vline(xintercept=rope[1], linetype=2, col="darkgreen") +
      geom_vline(xintercept=rope[2], linetype=2, col="darkgreen")
  }
  
  gplot <- gplot + labs(x="Sample", y="Posterior density")
  return(gplot)
}




#' @title Plotting the (marginal) posterior densities in Bayesian analyses
#'
#' @description This function plots, univariately, the posterior densities obtained
#' in Bayesian analyses. 
#' @param results A list containing, at least three elements, one named \code{approximate}, which is
#' a logical value indicating whether the posterior is a function or a sample, \code{rope}, a two dimensional
#' vector with the minimum and maximum values of the rope and \code{posterior}, either a one parameter function or
#' a matrix (or data.frame) where each row is a sample and each column a sampled parameter
#' @param parameter Either a string or a number indicating, in case the posterior is approximated, the parameter
#' to be ploted (i.e., the name or the index of a column in the sample matrix)
#' @param ... Additional parameters to the Rgraphviz function. This is mainly to change the layout of the graph
#' @param plot.rope  A logical value indicating whether the rope has to be plotted or not. Note that not for all
#' parameter the rope makes sense
#' @param num.points Number of points used to plot the functions
#' @param plot.samples A logical value. If true, the samples are plotted (only when the posterior is approximate)
#' @param alpha Numeric value for the transparency of the points, only applicable if \code{plot.samples} is true
#' @return An object of class \linkS4class{ggplot} with the plot
#' @details 
#' Note that if the methods are exact (not simulated), the true density can be plotted but,
#' for those cases where the posterior is approximated through sampling, the function will
#' plot a kernel density estimation of the posterior and, thus, the probabilities computed
#' by other functions are not directly the areas under the densities.
#' @examples
#' x <- rnorm(150, 1, 0.5)
#' y <- rnorm(150, 1, 0.05)
#' results <- bSignedRankTest(x=x, y=y, rope=c(-0.15, 0.15))
#' plotSimplex(results, posterior.label=TRUE)
#' 

plotSimplex <- function(results, A="Alg. A", B="Alg. B", plot.density=TRUE, plot.points=TRUE,
                        palette=c("green", "darkgray", "red"), point.size=1,
                        font.size=5, alpha=NULL, posterior.label=FALSE) {
  
  if (!requireNamespace("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it.", call.=FALSE)
  }
  if (!require("geometry", quietly=TRUE)) {
    stop("This function requires the geometry package. Please install it and try again.")
  }
  
  post.sample <- results$posterior
  post.sample <- post.sample[, c("Left","Rope","Right")]
  
  # Get the winner for the color of the points
  aux <- apply(post.sample, MARGIN=1, FUN=which.max)
  aux[aux==1] <- names(post.sample)[1]
  aux[aux==2] <- names(post.sample)[2]
  aux[aux==3] <- names(post.sample)[3]
  colors <- factor(aux, names(post.sample))
  
  # Coordinates of the eges of the Simplex
  simplex.coords <- rbind(c(2, 0), c(1,1), c(0, 0))
  
  # Auxiliar info to draw the triangle
  center <- c(1, 0.3333333)
  ab <- c(0.5, 0.5)
  bc <- c(1.5, 0.5)
  ac <- c(1, 0)
  
  #Convert from barycentric coords to cartesians and add the winner
  points <- data.frame(Color=colors, bary2cart(simplex.coords, as.matrix(post.sample)))
  names(points) <- c("Colors", "X", "Y")
  
  # Additional info for the plot (triangle and lines)
  triangle <- data.frame(simplex.coords[c(1, 2, 3), ])
  names(triangle) <- c("X", "Y")
  
  divisors1 <- data.frame(rbind(center, ab))
  names(divisors1) <- c("X", "Y")
  divisors2 <- data.frame(rbind(center, bc))
  names(divisors2) <- c("X", "Y")
  divisors3 <- data.frame(rbind(center, ac))
  names(divisors3) <- c("X", "Y")

  p.right <- results$posterior.probabilities["Right"]
  p.rope <- results$posterior.probabilities["Rope"]
  p.left <- results$posterior.probabilities["Left"]
  
  if (is.null(alpha)) {
    alpha <- min(1, 250/nrow(points))
  }
  
  # Create the plot, layer by layer
  g <- ggplot(points, aes(x=X, y=Y))
  
  # Optionally, add the points
  if (plot.points) {
    g <- g + geom_point(aes(color=Colors), alpha=alpha, shape=19, size=point.size) +
      scale_color_manual(values=c("Left"=palette[3], "Rope"=palette[2], "Right"=palette[1]), guide="none")
  }
  
  # Optionally, add the density
  if (plot.density) {
    g <- g + stat_density2d(aes(fill=..level..,alpha=..level..),geom="polygon", show.legend=FALSE) +
      scale_fill_gradient2(low="#ffffff", mid="#78c679", high="#005a32", guide="none") + 
      geom_polygon(data=data.frame(X=c(-0.1,1.1,-0.1), Y=c(-0.1,1.1,1.1)), fill="white") + 
      geom_polygon(data=data.frame(X=c(0.9,2.1,2.1), Y=c(1.1,1.1,-0.1)), fill="white") + 
      geom_polygon(data=data.frame(X=c(-0.1,2.1,2.1,-0.1), Y=c(0,0,-0.1,-0.1)), fill="white")
  }
  
  # Add the triagle and annotations
  g <- g + geom_polygon(data=triangle, color="black", fill=NA) + 
    geom_line(data=divisors1, color="black", linetype=2) + 
    geom_line(data=divisors2, color="black", linetype=2) + 
    geom_line(data=divisors3, color="black", linetype=2) + 
    annotate("text", x=0, y=-0.05, label=A, hjust=0, size=font.size) + 
    annotate("text", x=2, y=-0.05, label=B, hjust=1, size=font.size) + 
    annotate("text", x=1, y=1.05, label="Rope", hjust=0.5, size=font.size) +
    scale_y_continuous(limits=c(-0.1, 1.1)) + theme_void() + theme(panel.background=element_rect(fill="white"))
  
  
  # And, optionally, add the probabilities
  if (posterior.label) {
    g <- g + annotate("text", x=0, y=1, vjust=1, hjust=0, size=font.size,
                      label=paste0("P(", A, " Win)= ", round(p.right, 3), "\n",
                                   "P(", B, " Win)= ", round(p.left, 3), "\n",
                                   "P(Rope Win)= ", round(p.rope, 3), "\n"))
  } 
  
  return(g)
}




#' @title Plotting barycentric plots
#'
#' @description This function plots tuples of vectors that sum 1 in barycentric coordinates 
#' @param data.matrix a matrix or dataframe where each row is an n-dimension data point (whose sum is 1)
#' @return An object of class \linkS4class{ggplot} with the plot
#' @details 
#' @examples
#' data <- matrix(runif(100*5), ncol=5)
#' data <- data / rowSums(data)
#' plotBarycentric(data)
#' 


plotBarycentric <- function (data.matrix) {
  
  if (!require("ggplot2", quietly=TRUE)) {
    stop("This function requires the ggplot2 package. Please install it running\n\n     install.packages('ggplot2')")
  }
  if (!require("geometry", quietly=TRUE)) {
    stop("This function requires the geometry package. Please install it running\n\n     install.packages('ggplot2')")
  }
  
  
  if (is.null(colnames(data.matrix))) {
    colnames(data.matrix) <- paste("A", 1:ncol(data.matrix), sep="")
  }
  
  ## Beyond three dimensions the plot may be hard to interpret. To ease the interpretation we order the dimensions 
  ## according to the average value
  expected.val <- colMeans(data.matrix)
  o <- order(expected.val, decreasing=TRUE)
  
  data.matrix <- data.matrix[, o]
  
  ## For the coloring we collect the winners
  winner <- factor(apply(data.matrix, MARGIN=1, FUN=which.max), levels=1:ncol(data.matrix))
  levels(data.matrix) <- colnames(data.matrix)
  
  ## Build the plot
  card <- ncol(data.matrix)
  alpha.step <- 2*pi/card
  angle.seq <- seq(0, 2*pi-alpha.step, alpha.step)
  vertices.coords <- data.frame(X=sin(angle.seq), Y=cos(angle.seq), Angle=angle.seq, Algorithm=colnames(data.matrix))
  
  coords <- data.frame(Winner=winner, bary2cart(as.matrix(vertices.coords[, 1:2]), as.matrix(data.matrix)))
  
  g <- ggplot(vertices.coords, aes(x=X, y=Y))+ geom_polygon(fill="gray90", color="orange") +
    geom_point(data=coords, aes(col=Winner, fill=Winner), alpha=0.5, shape=16) + theme_void() +
    geom_label(data=vertices.coords, aes(label=Algorithm))
  return(g)
}
