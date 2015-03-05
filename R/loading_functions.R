## Function to process an experiment-type matrix

process.experiment.matrix <- function (data , alg.col , value.col){
  ## Verify the column names/id's  
  if (is.character(alg.col)){
    if (!alg.col %in% colnames(data)){
      stop (paste(alg.col , 
                  " column not found. Column names in the file are (", 
                  paste(colnames(data) , collapse = ",") , 
                  ")",sep=""))
    }else{
      alg.col <- which(colnames(data) %in% alg.col)
    }
  }else if(is.numeric(alg.col)){
    if (alg.col>dim(data)[2] | alg.col<=0){
      stop (paste("The alg.col parameter has to be a valid value (between 1 and " , dim(data)[2] , ")",sep=""))
    }
  }else{
    stop (paste("The alg.col parameter has to be either a number between 1 and ", 
                dim(data)[2] , " or the name of a column (", 
                paste(colnames(data) , collapse = ",") , 
                ")",sep=""))
  }
  
  if (is.character(value.col)){
    if (!value.col %in% colnames(data)){
      stop (paste(value.col , 
                  " column not found. Column names in the file are (", 
                  paste(colnames(data) , collapse = ",") , 
                  ")",sep=""))
    }else{
      value.col <- which(colnames(data) %in% value.col)
    }
  }else if(is.numeric(value.col)){ 
    if (value.col>dim(data)[2] | value.col<=0){
      stop (paste("The value.col parameter has to be a valid value (between 1 and " , dim(data)[2] , ")",sep=""))
    }
  }else{
    stop (paste("The value.col parameter has to be either a number between 1 and ", 
                dim(data)[2] , " or the name of a column (", 
                paste(colnames(data) , collapse = ",") , 
                ")",sep=""))
  }
  
  ## Process the file to build the final matrix
  descriptor.combinations <- unique(data[,-c(alg.col , value.col)])
  algorithms <- as.character(unique(data[,alg.col]))
  
  process.combination <- function (i){
    ## Subset of the whole dataset
    rows <- apply (data[,-c(alg.col , value.col)] , MARGIN = 1 , FUN = function(x)(all(x==descriptor.combinations[i,])))
    sub <- subset (data , rows)
    
    ## Process all the algorithms
    process.alg <- function(alg.name) subset(sub , sub[,alg.col]==alg.name)[,value.col]
    aux<-lapply(algorithms , process.alg)
    ## Check that all the vectors have the same length
    l.aux <- sapply(aux,length)
    if (length(unique(l.aux))>1){
      comb <- paste("(",paste(colnames(descriptor.combinations),collapse=",") , ") = (" ,paste(descriptor.combinations[i,],collapse=",") , ")",sep="")
      lengths <- paste("(",paste(algorithms,collapse=",") , ") = (" , paste(l.aux , collapse=",") , ")",sep="")
      mssg <- paste("Problems while parsing the file. For every combination of the parameters the algorithms should have the same number of values. In the combination ", comb , " the lengths associated to each algoithm are " , lengths  , sep="")
      stop(mssg)
    }
    res <- do.call(cbind,aux)
    colnames(res) <- algorithms
    suppressWarnings (cbind(descriptor.combinations[i,] , res))
  }
  
  aux <- lapply(1:nrow(descriptor.combinations) , FUN = process.combination)
  do.call(rbind , aux)  
}



#' @title Read data from an experiment-like file
#'
#' @export
#' @description This function reads the data from a file where each row is an experiment characterized by some variables (one of which should be the algorithm used) and with one and only one numeric result. for files where there is more than one result per line see \code{\link{read.comparison.file}}.
#' @param file Path to the file to read.
#' @param alg.col Name or index of the column corresponding to the algorithm used in the experiment.
#' @param value.col Name or index of the column corresponding to the numerical result of the experiment.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{read.experiment.dir}}, \code{\link{read.comparison.file}}, \code{\link{read.comparison.dir}}
#' @examples
#' dir <- system.file("loading_tests",package="scma")
#' file <- paste(dir , "beta_complete_experiment.out" , sep="/")
#' data <- read.experiment.file (file = file , alg.col = 'algorithm' , value.col = 'error')
#' dim(data)
#' head(data)

read.experiment.file <- function (file , alg.col , value.col , ...){
  data <- read.csv (file , ...)
  process.experiment.matrix(data , alg.col , value.col)
}

process.exp.file.in.dir <- function(file , names , alg.var.name , fname.pattern , first=2 , ...){
  
  symbs <- c('_',':',',',';','-','+','*','&','%','#')
  splt.id <- which(sapply(1:length(symbs) , function(i) length(grep(symbs[i] , fname))==0))[1]
  splt <- symbs[splt.id]
  
  fname <- basename(file)
  ## Transform the name into something easy to split and the separate it into the elements
  replacement <- paste(paste("\\",1:length(names),sep=""),collapse=splt)
  params <- strsplit(gsub(fname.pattern , replacement , fname),splt)[[1]]
  names(params)<-names
  data <- read.csv(file)
  if (ncol(data)!=1) stop("Only one column expected. For more than one column, trye 'read.comparison.dir' function") 
  if (!alg.var.name %in% names){
    if (alg.var.name != colnames(data)){
      stop(paste("The name " , alg.var.name , " not found neither in the file name nor in the header." , sep=""))
    } else{
      value.col <- 'value'
      alg.col <- 'Algorithm'
      alg <- colnames(data)
      colnames(data) <- value.col
      output <- cbind(matrix(rep(params,nrow(data)),ncol=length(params),byrow=T) , rep(alg , nrow(data)) , data)
      colnames(output) <- c(names , 'Algorithm' , value.col)
    }
  }else{
    value.col <- colnames(data)
    alg.col <- alg.var.name
    output <- cbind(matrix(rep(params,nrow(data)),ncol=length(params),byrow=T),data)
    colnames(output) <- c(names , value.col)
  }
  output
}

#' @title Read data from an experiment-like files in a directory
#'
#' @export
#' @description This function reads the data from all the files in a directory. Only one column is expected in each file, correspondig to the results obtained by a single algorithm. If the files contain the results of two or more algorithms, see function \code{\link{read.comparison.file}}. The function can extract information from the file name.
#' @param directory Directory with the files to load. It should only contain files to load, no other kind of file.
#' @param names List of names for the variables to be extracted from the file name
#' @param alg.var.name Name of the variable that defines the algorithm used in the experiment. It can be either one of the variables extracted from the file name or the header of the data in the file.
#' @param fname.pattern Regular expression to extract information from the file names. It has to be a regular expression that matches the name of the files and where the information to be extrcted has to be between brakets. As an example, if the whole file name wants to be used, the expression \code{'([.]*)'} can be used. For more example see the examples below or the vignette covering the data loading.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @details Note that all the files should have the same format (only one column with the same header)
#' @seealso \code{\link{read.experiment.file}}, \code{\link{read.comparison.file}}, \code{\link{read.comparison.dir}}
#' @examples
#' dir <- paste(system.file("loading_tests",package="scma") , "experiment_files" , sep="/")
#' ## The format of the files is beta_ALPHA,BETA_size_SIZE_ESTIMATOR.out, where variables to extract are in capital letters. 
#' list.files(dir)[1]
#' ## The regular expresion can ba as simple as substituting each variable name in the expression above by ([XXX]*), where XXX is the list of symbols that appear in the name.
#' pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*)_([a-z]*).out"
#' var.names <- c('alpha' , 'beta' , 'size', 'estimator')
#' data <- read.experiment.dir (dir , var.names , 'estimator' , pattern)
#' dim(data)
#' head(data)

read.experiment.dir <- function(directory , names , alg.var.name , fname.pattern , ...){
  
  first <- ifelse(substring(fname.pattern , 1 , 1)=="(" , 1 , 2)
  ## Load the first file to check the header name
  f <- list.files(directory)[1]
  data <- read.csv(paste(directory,f,sep="/") , ...)
  if (!alg.var.name %in% names){
    if (alg.var.name != colnames(data)){
      stop(paste("The name " , alg.var.name , " not found neither in the file name nor in the header." , sep=""))
    } else{
      value.col <- 'value'
      alg.col <- 'Algorithm'
    }
  }else{
    value.col <- colnames(data)
    alg.col <- alg.var.name
  }
  
  data <- data.frame()
  for (file in list.files(directory)){
    data <- rbind(data , process.exp.file.in.dir(paste(directory,file,sep="/") , names, alg.var.name , fname.pattern , first , ...))
  } 
  
  process.experiment.matrix(data , alg.col , value.col) 
}


#' @title Read data from a comparison file
#'
#' @export
#' @description This function reads the data from a files where two or more algorithms are compared in different problems. The file can have some columns that characterize the problem and one column per algorithm. If each row contain only the result obtained by one algorithm, use the \code{\link{read.experiment.file}} function.
#' @param file Path of the file to load
#' @param alg.names List of the column names that include the results. The rest are assumed as descriptors of the problems
#' @param col.names Vector of names of the columns. If not NULL, the files are assumed not to have a header and the columns are named using this vector
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{read.experiment.file}}, \code{\link{read.experiment.dir}}, \code{\link{read.comparison.dir}}
#' @examples
#' dir <- system.file("loading_tests",package="scma")
#' file <- paste(dir , "beta_complete_comparison.out" , sep="/")
#' data <- read.comparison.file (file = file , alg.names = c('kakizawa','vitale','boundarykernel','betakernel'))
#' dim(data)
#' head(data)
read.comparison.file <- function(file , alg.names , col.names=NULL , ...){
  header <- ifelse (is.null(col.names) , TRUE , FALSE)
  data <- read.csv (file , header = header , ...)
  if (!is.null(col.names)){
    if (ncol(data)!=col.names) stop ("The size of the table and the number of column names do not match")
    names(data)<-col.names
  }
  if (!all(alg.names %in% names(data))) stop ("Not all the algorithm names provided have been found in the file header")
  id.alg <- which(names(data) %in% alg.names)
  cbind(data[,-id.alg] , data[,id.alg])
}


## auxiliar function for read.comparison.dir
process.comp.file.in.dir <- function(file , col.names , alg.names , names , fname.pattern , first=2 , ...){
  
  symbs <- c('_',':',',',';','-','+','*','&','%','#')
  splt.id <- which(sapply(1:length(symbs) , function(i) length(grep(symbs[i] , fname))==0))[1]
  splt <- symbs[splt.id]
  
  fname <- basename(file)
  ## Transform the name into something easy to split and the separate it into the elements
  replacement <- paste(paste("\\",1:length(names),sep=""),collapse=splt)
  params <- strsplit(gsub(fname.pattern , replacement , fname),splt)[[1]]
  names(params)<-names
  data <- read.csv(file)
  if (!is.null(col.names)){
    if (ncol(data)!=col.names) stop ("The size of the table and the number of column names do not match")
    names(data)<-col.names
  }
  if (!all(alg.names %in% names(data))) stop ("Not all the algorithm names provided have been found in the file header")
  id.alg <- which(names(data) %in% alg.names)
  res <- cbind(matrix(rep(params,nrow(data)),ncol=length(params),byrow=T), data[,-id.alg] , data[,id.alg])
  names(res) <- c(names(params) , names(data)[-id.alg] , names(data)[id.alg])
  res
}


#' @title Read data from a directory of comparison-like files
#'
#' @export
#' @description This function reads the data from all files in a directory. Each file is expected to to be formated as a comparison file, i.e., the file can have some columns that characterize the problem and one column per algorithm. If each row contain only the result obtained by one algorithm, use the \code{\link{read.experiment.dir}} function.
#' @param directory Directory where the files are located.
#' @param alg.names List of the column names that include the results. The rest are assumed as descriptors of the problems.
#' @param col.names Vector of names of the columns. If not NULL, the files are assumed not to have a header and the columns are named using this vector.
#' @param names List of names for the variables to be extracted from the file name.
#' @param fname.pattern Regular expression to extract information from the file names. It has to be a regular expression that matches the name of the files and where the information to be extrcted has to be between brakets. As an example, if the whole file name wants to be used, the expression \code{'([.]*)'} can be used. For more example see the examples below or the vignette covering the data loading.
#' @param ... Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}).
#' @return A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
#' @seealso \code{\link{read.experiment.file}}, \code{\link{read.experiment.dir}}, \code{\link{read.comparison.dir}}
#' @examples
#' dir <- paste(system.file("loading_tests",package="scma") , "comparison_files" , sep="/")
#' ## The format of the files is beta_ALPHA,BETA_size_SIZE.out, where variables to extract are in capital letters. 
#' list.files(dir)[1]
#' ## The regular expresion can ba as simple as substituting each variable name in the expression above by ([XXX]*), where XXX is the list of symbols that appear in the name.
#' pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*).out"
#' var.names <- c('alpha' , 'beta' , 'size')
#' alg.names <- c('kakizawa','vitale','boundarykernel','betakernel')
#' data <- read.comparison.dir (directory = dir , alg.names = alg.names , col.names = NULL , names = names , fname.pattern = pattern)
#' dim(data)
#' head(data)

read.comparison.dir <- function (directory , alg.names , col.names , names , fname.pattern , ...){
  first <- ifelse(substring(fname.pattern , 1 , 1)=="(" , 1 , 2)
  ## Load the first file to check the header name
  f <- list.files(directory)[1]
  
  data <- data.frame()
  for (file in list.files(directory)){
    data <- rbind(data , process.comp.file.in.dir(paste(directory,file,sep="/") , col.names , alg.names , names,  fname.pattern , first , ...))
  } 
  data
}
