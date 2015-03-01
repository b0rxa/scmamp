# NORMALITY CHECK PLOTS ---------------------------------------------------

#' @title Gaussian distribution quantile-quantile plot
#'
#' @description This function creates a quantile-quantile plot to assess the goodness of fit of a Gaussian distribution to a given sample.
#' @param sample List of data points to check
#' @param ... The plot is created using \code{\link{ggplot2}}. This special parameter can be used to pass additional parameters to the \code{\link{geom_point}} function used to plot the sample points.
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{scma} \code{plot.densities}

gaussian.qqplot <- function (sample , ...){
  if(!require(ggplot2)) stop("This function requires the package ggplot2, which is not installed. You can install it typing install.packages('ggplot2')")
  sample <- sort(sample)
  m <- mean(sample)
  sd <- sd(sample)
  n <- length(sample)
  q.gaussian <- qnorm(seq(1/n , (n-1)/n , 1/n) , m , sd)
  q.sample <- sample[-n]
  df <- data.frame(Empirical = q.sample , Gaussian = q.gaussian)  
  m <- min(df)
  M <- max(df)
  
  ggplot(df , aes(x=Empirical , y=Gaussian)) + scale_x_continuous(limits = c(m,M))+ scale_y_continuous(limits = c(m,M)) + geom_abline(slope=1 , intercept = 0 , col="darkgray" , size=1.1) +  geom_point(...) 
}


#' @title Kernel based density estimation of the samples
#'
#' @description This function estimates and plots the densities of the results of each algorithm
#' @param results.matrix A matrix where columns represent the algorithms
#' @param ... The plot is created using \code{\link{ggplot2}}. This special parameter can be used to pass additional parameters to the \code{\link{geom_line}} function used to plot the sample points.
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{scma} \code{gaussian.qqplot}

plot.densities <- function (results.matrix , ...){
  if(!require(ggplot2)) stop("This function requires the package ggplot2, which is not installed. You can install it typing install.packages('ggplot2')")
  k <- dim(results.matrix)[2]
  aux <- lapply(1:k , FUN = function (x) {
    d <- density(results.matrix[,x])
    data.frame (Algorithm = colnames(results.matrix)[x] , Value = d$x , Density = d$y)
  })
  df <- do.call(rbind , aux)
  ggplot(df , aes(x = Value , y = Density , col = Algorithm)) + geom_line(...)
  
}

# PLOTTING THE RESULTS ----------------------------------------------------

#' @title Plotting the p-value matrix
#'
#' @description This function plots the p-value matrix as a \code{\link{ggplot2}}'s tile plot.
#' @param pvalue.matrix Matrix with the p-values to plot
#' @param alg.order A permutation indicating the reordering of the algorithms
#' @param show.pvalue Logical value indicating whether the numerical values have to be printed
#' @param font.size Size of the p-values, if printed
#' @return A \code{\link{ggplot}} object.
#' @seealso \code{scma} \code{algorithm.graph}

plot.pvalues <- function(pvalue.matrix , alg.order = NULL , show.pvalue = TRUE , font.size = 5){
  if(!require(reshape2)) stop("This function requires the package ggplot2, which is not installed. You can install it typing install.packages('reshape2')")
  if(!require(ggplot2)) stop("This function requires the package ggplot2, which is not installed. You can install it typing install.packages('ggplot2')")
  
  df <- melt(pvalue.matrix)
  colnames(df) <- c("X" , "Y" , "p.value")
  if (!is.null(alg.order)){
    l <- colnames(pvalue.matrix)[alg.order]
    df$X <- factor(df$X , levels = l)
    df$Y <- factor(df$Y , levels = l)
  }
  g <- ggplot(df, aes(x=X,y=Y,fill=p.value)) + geom_tile(col="white") +
          scale_fill_continuous("p-value") +
          labs(x="Algorithm" , y="Algorithm")
  if (show.pvalue) g <- g + geom_text(aes(label = round(p.value,2)) , size = font.size , col="white")
  g
} 

critical.difference.plot <- function (results.matrix , alpha = 0.05 , cex=0.75){
  k <- dim(results.matrix)[2]
  N <- dim(results.matrix)[1]
  nem <- nemenyi.test (data = results.matrix , alpha = alpha)
  cd <- nem$statistic
  mean.rank <- sort(colMeans(rank.matrix(results.matrix)))
  
  ## Separate the algorithms in left and right parts
  lp <- round(k/2)
  left.algs <- mean.rank[1:lp]
  right.algs <- mean.rank[(lp+1):k]  
  max.rows <- ceiling(k/2)
  
  ## Basic dimensions and definitions
  char.size <- 0.001 ## Character size
  line.spacing <- 0.25 ## Line spacing for the algorithm name
  m <- floor(min(mean.rank))
  M <- ceiling(max(mean.rank))
  max.char <- max(sapply(colnames(results.matrix), FUN = nchar)) ## Longest length of a label
  text.width <- (max.char+4) *char.size
  w <- (M-m) + 2*text.width
  h.up <- 2.5 * line.spacing                ## The upper part is fixed
  h.down <- (max.rows + 2.25) * line.spacing ## The lower part depends on the number of algorithms. The 2 extra spaces are for the lines that join algorithms
  tick.h <- 0.25 * line.spacing
  label.displacement <- 0.1 ## Displacement of the label with respect to the axis
  line.displacement <- 0.025 ## Displacement for the lines that join algorithms
  
  ## Background
  plot(0 , 0 , type='n' , xlim = c(m - w/(M-m) , M + w/(M-m)) , ylim=c(-h.down , h.up) , 
       xaxt='n' , yaxt = 'n' , xlab = '' , ylab = '' , bty = 'n')
  
  ## Draw axis
  lines (c(m,M) , c(0,0))
  dk <- sapply(m:M , function(x) {
    lines(c(x,x) , c(0 , tick.h))
    text(x , 3*tick.h , labels = x , cex = cex)
    })
  
  ## Draw CD
  lines(c(m , m + cd) , c(1.75*line.spacing , 1.75*line.spacing))
  text(m+cd/2 , 2.25*line.spacing , "CD" , cex = cex)
  lines(c(m,m) ,  c(1.75*line.spacing-tick.h/4 , 1.75*line.spacing+tick.h/4))
  lines(c(m + cd,m + cd) , c(1.75*line.spacing-tick.h/4 , 1.75*line.spacing+tick.h/4))
  
  ## Left part, labels
  dk <- sapply (1:length(left.algs) , function(x) {
    line.h <- -line.spacing*(x+2)
    text(m-label.displacement , line.h , names(left.algs)[x] , cex=cex , adj = 1)
    lines(c(m , left.algs[x]) , c(line.h , line.h))
    lines(c(left.algs[x] , left.algs[x]) , c(line.h , 0))
  })
  
  ## Right part, labels
  dk <- sapply (1:length(right.algs) , function(x) {
    line.h <- -line.spacing*(x+2)
    text(M+label.displacement , line.h , names(right.algs)[x] , cex=cex , adj = 0)
    lines(c(M , right.algs[x]) , c(line.h , line.h))
    lines(c(right.algs[x] , right.algs[x]) , c(line.h , 0))
  })
  
  ## Draw the lines to join algorithms
  get.interval <- function (x){
    from <- mean.rank[x]
    diff <- mean.rank - from
    ls <- which(diff>0 & diff<cd)
    if (length(ls)>0)
    {
      c(from , mean.rank[max(ls)])
    }
  }
  to.join <- do.call(rbind,mapply (1:k , FUN = get.interval))
  row <- c(1)
  ## Determine each lin in which row will be displayed
  nlines <- dim(to.join)[1]
  for(r in 2:nlines){
    id <- which(to.join[r,1]>to.join[,2])
    if(length(id)==0){
      row <- c(row , tail(row,1)+1)
    }else{
      row <- c(row , min(row[id]))
    }
  }
  
  step <- max(row)/2
  
  ## Draw the line
  dk<-sapply (1:nlines , FUN = function(x){
    y <- -line.spacing*(0.5 + row[x]/step)
    lines(c(to.join[x,1]-line.displacement , to.join[x,2]+line.displacement) , c(y,y) , lwd=3)
  })
  
}


algorithm.graph <- function (hypothesis.matrix , mean.value , ... , font.size = 15 , highlight="min" , highlight.color = "chartreuse3" , node.color="gray30" , font.color="white", digits = 2 , node.width = 5 , node.height=2){
  if(!require(Rgraphviz)) stop("This function requires the package ggplot2, which is not installed. You can install it typing source('http://www.bioconductor.org/biocLite.R') and then biocLite('Rgraphviz')")
  
  if (!all(colnames(hypothesis.matrix) %in% names(mean.value))) stop ("The names of the algorithms in the matrix and the mean.valu vector do not match")
  
  nc <- rep(node.color , length(mean.value))
  if(highlight=="min"){
    nc [which.min(mean.value)] <- highlight.color
  }else if(highlight=="max"){
    nc [which.max(mean.value)] <- highlight.color
  }
    
  hypothesis.matrix[is.na(hypothesis.matrix)]<-FALSE
  adj.matrix <- hypothesis.matrix
  colnames(adj.matrix) <- rownames(adj.matrix) <- names(mean.value)
  
  nl <- paste(names(mean.value) , "\\\n" , round(mean.value,digits) , sep="")
  
  am.graph <- new("graphAM" , adjMat = adj.matrix , edgemode="undirected")
  names(nc) <- names(nl) <- nodes(am.graph)
  nAttrs <- list()
  nAttrs$label <- nl
  nAttrs$fillcolor <- nc
  attrs <- list(node = list(shape = "rectangle" , width = node.width , height = node.height ,
                            fontcolor = font.color , fontsize = font.size))
  plot(am.graph , ... , nodeAttrs = nAttrs , attrs = attrs)
}