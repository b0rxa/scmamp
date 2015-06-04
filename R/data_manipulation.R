#' @title Expression based row filtering
#'
#' @export
#' @description This is a simple function to filter data based on an expression defined using the colum names
#' @param data A NAMED matrix or data frame to be filtered (column names are required).
#' @param condition A string indicating the condition that the row have to fulfill to be retained. The column names are used as variables in the condition (see examples bellow).
#' @param remove.cols Either a vector of column names or a vector of column indices to be removed from the result
#' @return The original data where the rows for which the condition is \code{FALSE} and the columns in the vector \code{remove.cols}  have been removed
#' @seealso \code{\link{summarizeData}}, \code{\link{writeTabular}} and the vignette \code{vignette(topic="Data_loading_and_manipulation", 
#' package="scmamp")}
#' @examples
#' data(data_gh_2008)
#' names(data.gh.2008)
#' filterData(data.gh.2008, condition="CN2 > 0.7 & Kernel < 0.7", remove.cols=1:2)
#' 
filterData <- function (data, condition="TRUE", remove.cols=NULL) {
  checkRow <- function (row) {
    # Extract columns as variables
    for (i in seq(along.with=row)) {
      assign(names(row)[i], row[i])
    }
    # Evaluate the condition
    cond <- eval(parse(text=condition))
    return(cond)
  }
  
  # Generate the subset of rows
  sub <- apply(data, MARGIN=1, FUN=checkRow)
  
  ## Generate the colums to select
  if (is.character(remove.cols)) {
    id.retain <- which(!(colnames(data) %in% remove.cols))
  } else {
    id.retain <- which(!(1:ncol(data) %in% remove.cols))
  }
  
  # In case there are indices out of range, remove them
  id.retain <- subset(id.retain, 
                      subset=id.retain > 0 & id.retain <= ncol(data))
  
  # Get the subset
  sbst <- subset(data, subset=sub, select=id.retain)
  return(sbst)
}


#' @title Summarization of data
#'
#' @export
#' @description This is a simple function to apply a summarization function to a matrix or data frame.
#' @param data A matrix or data frame to be summarized.
#' @param fun Function to be used in the summarization. It can be any function that, taking as first argument a numeric vector, otuputs a numeric value. Typical examples are \code{\link{mean}}, \code{\link{median}}, \code{\link{min}}, \code{\link{max}} or \code{\link{sd}}.
#' @param group.by A vector of either column names or column indices according to which the data will be grouped to be summarized.
#' @param ignore A vector of either column names or column indices of the columns that have to be removed from the output.
#' @param ... Additional parameters to the summarization function (\code{fun}). For example, \code{na.rm=TRUE} to indicate that the missing values should be ignored.
#' @return A data frame where, for each combination of the values in the columns indicated by \code{group.by}, each column (except those in \code{ignore}) contains the summarization of the values in the original matrix that have that combination of values.
#' #' @seealso \code{\link{filterData}}, \code{\link{writeTabular}} and the vignette \code{vignette(topic="Data_loading_and_manipulation", 
#' package="scmamp")}
#' @examples
#' data(data_blum_2015)
#' # Group by size and radius. Get the mean and variance of only the last two
#' # columns.
#' summarizeData (data.blum.2015, group.by=c("Radius","Size"), ignore=3:8, 
#'                fun=mean, na.rm=TRUE)
#' summarizeData (data.blum.2015, group.by=c("Radius","Size"), ignore=3:8, 
#'                fun=sd, na.rm=TRUE)
#' 
summarizeData <- function (data, fun=mean, group.by=NULL, ignore=NULL, ... ) {
  if (!is.data.frame(data)) {
    data <- data.frame(data)
  }
  # Convert character definitions to colum id
  if (is.character(group.by)) {
    group.by <- which(colnames(data) %in% group.by)
  }
  
  if (is.character(ignore)) {
    ignore <- which(colnames(data) %in% ignore)
  }
  
  ## Only numeric columns can be summarized
  non.numeric <- which(!unlist(lapply(data, is.numeric)))
  if (!all(non.numeric %in% c(group.by, ignore))) {
    warning ("Only numeric columns can be summarized. Character and factor ",
             "columns should be either in the 'group.by' or the 'ignore' list. ",
             "Non numeric columns will be ignored")
    ignore <- unique(c(ignore, non.numeric))
  }
  
  # Remove any index out of bounds
  group.by <- subset(group.by, subset=group.by > 0 & group.by <= ncol(data))
  ignore <- subset(ignore, subset=ignore > 0 & ignore <= ncol(data))
  
  if (length(intersect(group.by,ignore)) > 0) {
    stop("The same column cannot be simultaneously in the 'group.by' and the ",
         "'ignore' list")
  }
  
  if (is.null(group.by)) {
    if (!is.null(ignore)) {
      data <- data[, -ignore]
    }
    summ <- apply(data, MARGIN=2, 
                  FUN=function(x) {
                    fun(x, ...)
                  }) 
  }else{  
    groups <- unique(data[, group.by])
    if(length(group.by)) groups <- data.frame(groups)
    to.summarize <- (1:ncol(data))[-c(ignore, group.by)]
    summGroup <- function (i) {
      sub <- rep(TRUE, nrow(data))
      for (j in seq(along.with=group.by)) {
        sub <- sub & data[, group.by[j]] == groups[i,j]
      }
      m <- subset(data, subset=sub)
      m <- m[, to.summarize]
      if (length(to.summarize) == 1) {
        m <- matrix(m, ncol=1)
      }
      apply(m, MARGIN=2, 
            FUN=function(x) {
                  fun(x , ...)
                })
    }
    aux <- lapply(1:nrow(groups), FUN=summGroup)
    summ <- cbind(groups, do.call(rbind, aux))
  }
  return(summ)
}