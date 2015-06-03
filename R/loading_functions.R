# NON-EXPORTED, AUXILIAR FUNCTIONS --------------------------------------------

processExperimentMatrix <- function (data, alg.col, value.col) {
  # Auxiliar function to process an experiment matrix
  # Args:
  #   data: Data to be processed
  #   alg.col:   Names or id of the columns containing the results
  #   value.col: Name or id of the column containing the values
  # Returns:
  #   Processed matrix
  #
  
  # Verify the column names/id's  
  if (is.character(alg.col)) {
    if (!alg.col %in% colnames(data)) {
      stop (alg.col, "column not found. Column names in the file are (", 
                  paste(colnames(data), collapse=","),  ")", sep="")
    } else {
      alg.col <- which(colnames(data) %in% alg.col)
    }
  } else if(is.numeric(alg.col)) {
    if (alg.col > dim(data)[2] | alg.col <= 0) {
      stop ("The alg.col parameter has to be a valid value (between 1 and ", 
            dim(data)[2], ")", sep="")
    }
  } else {
    stop ("The alg.col parameter has to be either a number between 1 and ",
          dim(data)[2], " or the name of a column (",
          paste(colnames(data), collapse=","), ")", sep="")
  }
  
  if (is.character(value.col)) {
    if (!value.col %in% colnames(data)) {
      stop (paste(value.col, " column not found. Column names in the file are (", 
                  paste(colnames(data), collapse=","), ")",sep=""))
    } else {
      value.col <- which(colnames(data) %in% value.col)
    }
  } else if (is.numeric(value.col)) { 
    if (value.col > dim(data)[2] | value.col <= 0) {
      stop ("The value.col parameter has to be a valid value (between 1 and " , 
            dim(data)[2], ")", sep="")
    }
  } else {
    stop ("The value.col parameter has to be either a number between 1 and ", 
          dim(data)[2], " or the name of a column (", 
          paste(colnames(data), collapse = ","), ")",sep="")
  }
  
  # Process the file to build the final matrix
  
  grouping <- which(!1:ncol(data) %in% c(alg.col, value.col))
  groups <- unique(data[, grouping])
  
  algorithms <- as.character(unique(data[, alg.col]))
  processCombination <- function (i) {
    # Subset of the whole dataset
    
    # Using the apply function here is not efficient. A loop in the 'grouping'
    # variables is far more computationally efficient.
    
    rows <- rep(TRUE, nrow(data))
    for (j in seq(along.with=grouping)) {
      g <- grouping[j]
      rows <- rows & data[, g]==groups[i, j]
    }
    sub <- subset (data, rows)
    
    # Process all the algorithms
    processAlg <- function(alg.name) {
      sb <- subset(sub, sub[, alg.col] == alg.name)[, value.col]
      return(sb)
    }
    
    aux <- lapply(algorithms, processAlg)
    
    # Check that all the vectors have the same length
    l.aux <- sapply(aux, length)
    if (length(unique(l.aux)) > 1) {
      comb <- paste("(", paste(colnames(groups), collapse=","), 
                    ") = (", paste(groups[i, ], collapse=","), 
                    ")",sep="")
      lengths <- paste("(", paste(algorithms, collapse=","), 
                       ") = (", paste(l.aux , collapse=","), ")", sep="")
      mssg <- paste("Problems while parsing the file. For every combination ",
                    "of the parameters the algorithms should have the same ","
                    number of values. In the combination ", comb, " the lengths ",
                    "associated to each algoithm are ", lengths, sep="")
      stop(mssg)
    }
    res <- do.call(cbind, aux)
    colnames(res) <- algorithms
    suppressWarnings (expr={
                        cbind(groups[i,], res)
                      })
  }
  
  aux <- lapply(1:nrow(groups), FUN=processCombination)
  return(do.call(rbind, aux))
}



processExpFile <- function(file, fname.pattern, names, alg.var.name , value.col, col.names=NULL, ...){
  # Auxiliar function to process individual experiment files in a directory
  # Args:
  #   file:          Path of the file to process
  #   fname.pattern: Pattern to extract information from the file name
  #   names:         Vector of names for the values extracted from the file name
  #   alg.var.name:  Name of the variable (either a column or an extracted value) 
  #                  containing the information about the algorithm used
  #   value.col:     Column containing the results
  #   col.names:     Names for the columns in the file. If NULL, the first row 
  #                  in the file is used as name
  #   
  #   Returns:
  #     Matrix with the information read from the file
  #
  
  
  # Process the name of the file
  fname <- basename(file)
  symbs <- c('_', ':', ',', ';', '-', '+', '*', '&', '%', '#')
  chars.in.name <- unlist(strsplit(fname, split=vector()))
  splt.id <- which(!(symbs %in% chars.in.name))[[1]]
  splt <- symbs[splt.id]
  
  # Transform the name into something easy to split and the separate it into the elements
  replacement <- paste(paste("\\", 1:length(names), sep=""), collapse=splt)
  params <- strsplit(gsub(fname.pattern, replacement, fname), splt)[[1]]
  names(params)<-names
  
  if (!is.null(col.names)) {
    data <- read.csv(file, header=FALSE, ...)
    if (length(col.names) != ncol(data)) {
      stop("The number of columns (", ncol(data), ") in the file does not match ",
           "the length of 'col.names' (", length(col.names), ")")
    }
    colnames(data) <- col.names
  } else {
    data <- read.csv(file, header=TRUE, ...)
  }
  
  # Merge the info get from the file name with that inside it
  output <- cbind(matrix(rep(params, nrow(data)), ncol=length(params), byrow=TRUE),
                  data)
  colnames(output) <- c(names, colnames(data))
  return(output)
}


processCompFile <- function(file, fname.pattern, names, alg.cols, col.names, ...){
  # Auxiliar function to process individual comparison files in a directory
  # Args:
  #   file:          Path of the file to process
  #   fname.pattern: Pattern to extract information from the file name
  #   names:         Vector of names for the values extracted from the file name
  #   alg.cols:      Columns containing the results of the algorithms
  #   col.names:     Names for the columns in the file. If NULL, the first row 
  #                  in the file is used as name
  #   
  #   Returns:
  #     Matrix with the information read from the file
  #
  
  rcsv.args <- list(...)
  if (!is.null(rcsv.args$header)) {
    stop("The argument header cannot be set by hand. It depends on whether a ",
         "col.names argument is passed or not")
  }
  
  fname <- basename(file)
  symbs <- c('_',':',',',';','-','+','*','&','%','#')
  chars.in.name <- unlist(strsplit(fname, split=vector()))
  splt.id <- which(!(symbs %in% chars.in.name))[[1]]
  splt <- symbs[splt.id]
  
  # Transform the name into something easy to split and the separate it into the elements
  replacement <- paste(paste("\\", 1:length(names), sep=""), collapse=splt)
  params <- strsplit(gsub(fname.pattern, replacement, fname), splt)[[1]]
  names(params) <- names
  if(is.null(col.names)){
    header <- TRUE
  } else { 
    header <- FALSE
  }
  
  data <- read.csv(file, header=header, ...)
  if (!is.null(col.names)) {
    if (ncol(data) != length(col.names)) {
      stop ("The size of the table and the number of column names do not match")
    }
    names(data) <- col.names
  }
  if(is.character(alg.cols)) {
    aux <- which(names(data) %in% alg.cols)
  }else{
    aux <- alg.cols
  }
  
  id.alg <- subset(aux, subset=((aux > 0) & (aux <= ncol(data))))
  if (length(id.alg) != length(alg.cols)) {
    stop ("Not all the algorithm names provided have been found in the ",
          "file header")
  }
  aux.matrix <- matrix(rep(params, nrow(data)), ncol=length(params), byrow=T)
  res <- cbind(aux.matrix, data[, -id.alg], data[, id.alg])
  names(res) <- c(names(params), names(data)[-id.alg], names(data)[id.alg])
  return(res)
}


# EXPORTED FUNCTIONS -----------------------------------------------------------


#' @title Read data from an experiment-like file
#'
#' @export
#' @description This function reads the data from a file where each row is an experiment characterized by some variables (one of which should be the algorithm used) and with one and only one numeric result. for files where there is more than one result per line see \code{\link{readComparisonFile}}.
#' @param file Path to the file to read.
#' @param alg.col Name or index of the column corresponding to the algorithm used in the experiment.
#' @param value.col Name or index of the column corresponding to the numerical result of the experiment.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{readExperimentDir}}, \code{\link{readComparisonFile}}, \code{\link{readComparisonDir}}
#' @examples
#' dir <- system.file("loading_tests",package="scmamp")
#' file <- paste(dir , "beta_complete_experiment.out" , sep="/")
#' data <- readExperimentFile (file = file , alg.col = 'algorithm' , value.col = 'error')
#' dim(data)
#' head(data)

readExperimentFile <- function (file, alg.col, value.col, col.names=NULL, ...) {
  rcsv.args <- list(...)
  if (!is.null(rcsv.args$header)) {
    stop("The argument header cannot be set by hand. It depends on whether the ",
         "col.names argument is passed or not")
  }
  if (!is.null(col.names)) {
    if (length(col.names) != ncol(data)) {
      stop("The number of columns (", ncol(data), 
           ") in the file does not match the length of 'col.names' (",
           length(col.names), ")")
    }
    data <- read.csv(file, header=FALSE, ...)
    colnames(data) <- col.names
  } else {
    data <- read.csv(file, header=TRUE, ...)
  }
  data <- processExperimentMatrix(data, alg.col, value.col)
  return(data)
}


#' @title Read data from an experiment-like files in a directory
#'
#' @export
#' @description This function reads the data from all the files in a directory. Only one column is expected in each file, correspondig to the results obtained by a single algorithm. If the files contain the results of two or more algorithms, see function \code{\link{readComparisonFile}}. The function can extract information from the file name.
#' @param directory Directory with the files to load. It should only contain files to load, no other kind of file.
#' @param names List of names for the variables to be extracted from the file name
#' @param alg.var.name Name of the variable that defines the algorithm used in the experiment. It can be either one of the variables extracted from the file name or the header of the data in the file.
#' @param value.col Name or index (referred to the column in the file) of the column indicating the value for the algorithm
#' @param fname.pattern Regular expression to extract information from the file names. It has to be a regular expression that matches the name of the files and where the information to be extrcted has to be between brakets. As an example, if the whole file name wants to be used, the expression \code{'([.]*)'} can be used. For more example see the examples below or the vignette covering the data loading.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @details Note that all the files should have the same format (only one column with the same header)
#' @seealso \code{\link{readExperimentFile}}, \code{\link{readComparisonFile}}, \code{\link{readComparisonDir}}
#' @examples
#' dir <- paste(system.file("loading_tests",package="scmamp") , "experiment_files" , sep="/")
#' ## The format of the files is beta_ALPHA,BETA_size_SIZE_ESTIMATOR.out, where variables to extract are in capital letters. 
#' list.files(dir)[1]
#' ## The regular expresion can ba as simple as substituting each variable name in the expression above by ([XXX]*), where XXX is the list of symbols that appear in the name.
#' pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*)_([a-z]*).out"
#' var.names <- c('alpha' , 'beta' , 'size', 'estimator')
#' data <- readExperimentDir (dir , var.names , 'estimator' , pattern)
#' dim(data)
#' head(data)

readExperimentDir <- function(directory, names, fname.pattern, alg.var.name, 
                              value.col, col.names=NULL, ...){
  rcsv.args <- list(...)
  if (!is.null(rcsv.args$header)) {
    stop("The argument header cannot be set by hand. It depends on whether a ",
         "col.names argument is passed or not")
  }
  
  if(!is.character(alg.var.name)) {
    stop("This function only accepts a name as the column indicating the ",
         "algorithm ('alg.var.name' argument)")
  }
  
  if (length(alg.var.name) != 1 | length(value.col) != 1) {
    stop ("The 'alg.var.name' and 'value.col' have to be of dimension 1")
  }
  
  if (substring(fname.pattern, 1, 1) == "(") {
    first <- 1
  } else {
    first <- 2
  }

  # Load the first file to check the header name
  f <- list.files(directory)[1]
  if (is.null(col.names)) {
    data <- read.csv(paste(directory, f, sep="/"), header=TRUE, ...)
  } else {
    data <- read.csv(paste(directory, f, sep="/"), header=FALSE, ...)
    if (length(col.names) != ncol(data)) {
      stop("The number of columns (", ncol(data), ") in the file does not ", 
           "match the length of 'col.names' (", length(col.names), ")")
    }
    colnames(data) <- col.names
  }
  
  if ((!alg.var.name %in% names) & (!alg.var.name %in% colnames(data))) {
      stop("The name ", alg.var.name, " not found neither in the file name ","
           nor in the header.", sep="")
  }
  
  if (is.character(value.col)) {
    if (!value.col %in% colnames(data)) {
      stop("Column named ", value.col, " not found in the files")
    }
  } else {
    if (value.col < 1 | value.col > ncol(data)) {
      stop("The column index ", value.col, " is out of the range of the file")
    }else{
      value.col <- colnames(data)[value.col]
    }
  }
  
  data <- data.frame()
  for (file in list.files(directory)) {
    data.new <- processExpFile(file=paste(directory, file, sep="/"), 
                               fname.pattern=fname.pattern, names=names, 
                               alg.var.name=alg.var.name, value.col=value.col,
                               col.names=col.names, ...)
    data <- rbind(data, data.new)
  } 
  
  d <- processExperimentMatrix(data=data, alg.col=alg.var.name, value.col=value.col) 
  return(d)
}


#' @title Read data from a comparison file
#'
#' @export
#' @description This function reads the data from a files where two or more algorithms are compared in different problems. The file can have some columns that characterize the problem and one column per algorithm. If each row contain only the result obtained by one algorithm, use the \code{\link{readExperimentFile}} function.
#' @param file Path of the file to load
#' @param alg.cols A vector column names or indices inicating which columns contain the results. The rest are assumed as descriptors of the problems
#' @param col.names Vector of names of the columns. If not NULL, the files are assumed not to have a header and the columns are named using this vector
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{readExperimentFile}}, \code{\link{readExperimentDir}}, \code{\link{readComparisonDir}}
#' @examples
#' dir <- system.file("loading_tests",package="scmamp")
#' file <- paste(dir , "beta_complete_comparison.out" , sep="/")
#' data <- readComparisonFile (file = file , alg.cols = c('kakizawa','vitale','boundarykernel','betakernel'))
#' dim(data)
#' head(data)
readComparisonFile <- function(file, alg.cols, col.names=NULL, ...) {
  rcsv.args <- list(...)
  if (!is.null(rcsv.args$header)) {
    stop("The argument header cannot be set by hand. It depends on whether a ",
         "col.names argument is passed or not")
  }
  
  if (is.null(col.names)) {
    header <- TRUE
  } else {
    header <- FALSE
  }
  
  data <- read.csv (file, header=header, ...)
  if (!is.null(col.names)) {
    if (ncol(data) != col.names) {
      stop ("The size of the table and the number of column names do not match")
    }
    names(data) <- col.names
  }
  
  if(is.character(alg.cols)) {
    aux <- which(names(data) %in% alg.cols)
  } else {
    aux <- alg.cols
  }
  id.alg <- subset(aux, subset=((aux > 0) & (aux <= ncol(data))))
  if (length(id.alg) != length(alg.cols)) {
    stop ("Not all the algorithm names provided have been found in the file header")
  }
  
  res <- cbind(data[, -id.alg], data[, id.alg])
  return(res)
}


#' @title Read data from a directory of comparison-like files
#'
#' @export
#' @description This function reads the data from all files in a directory. Each file is expected to to be formated as a comparison file, i.e., the file can have some columns that characterize the problem and one column per algorithm. If each row contain only the result obtained by one algorithm, use the \code{\link{readExperimentDir}} function.
#' @param directory Directory where the files are located.
#' @param alg.cols A vector column names or indices inicating which columns contain the results. The rest are assumed as descriptors of the problems
#' @param col.names Vector of names of the columns. If not NULL, the files are assumed not to have a header and the columns are named using this vector.
#' @param names List of names for the variables to be extracted from the file name.
#' @param fname.pattern Regular expression to extract information from the file names. It has to be a regular expression that matches the name of the files and where the information to be extrcted has to be between brakets. As an example, if the whole file name wants to be used, the expression \code{'([.]*)'} can be used. For more example see the examples below or the vignette covering the data loading.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{readExperimentFile}}, \code{\link{readExperimentDir}}, \code{\link{readComparisonDir}}
#' @examples
#' dir <- paste(system.file("loading_tests",package="scmamp") , "comparison_files" , sep="/")
#' ## The format of the files is beta_ALPHA,BETA_size_SIZE.out, where variables to extract are in capital letters. 
#' list.files(dir)[1]
#' ## The regular expresion can ba as simple as substituting each variable name in the expression above by ([XXX]*), where XXX is the list of symbols that appear in the name.
#' pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*).out"
#' var.names <- c('alpha' , 'beta' , 'size')
#' alg.names <- c('kakizawa','vitale','boundarykernel','betakernel')
#' data <- readComparisonDir (directory = dir , alg.cols = alg.names , col.names = NULL , names = var.names , fname.pattern = pattern)
#' dim(data)
#' head(data)

readComparisonDir <- function (directory , alg.cols , col.names , names , fname.pattern , ...){
  rcsv.args <- list(...)
  if (!is.null(rcsv.args$header)) stop("The argument header cannot be set by hand. It depends on whether a col.names argument is passed or not")
  
  first <- ifelse(substring(fname.pattern , 1 , 1)=="(" , 1 , 2)
  ## Load the first file to check the header name
  f <- list.files(directory)[1]
  
  data <- data.frame()
  for (file in list.files(directory)){
    data <- rbind(data, 
                  processCompFile(file=paste(directory,file,sep="/"),
                                  fname.pattern=fname.pattern, names=names, 
                                  alg.cols=alg.cols, col.names=col.names, ...))
  } 
  data
}
