## Function to write a single cell
print.cell <- function (value , bold , italic , format , digits , mark , mark.char , ...){
  value <- unlist(value)
  ## Get the string to print
  if (is.factor(value)) value <- as.character(value)
  suppressWarnings(value.num <- as.numeric(value))
  if (!is.na(value.num)){
    value <- as.numeric(value)
    strv <- formatC(value.num , digits=digits , format=format)
  }else{
    strv <- ifelse (is.na(value) , "" , value)
  }
  
  ## Check mark
  if (mark){
    strv <- paste(strv , "$^" , mark.char , "$",sep="")
  }
  
  ## Check bold and/or italic

  if (bold | italic){
    left <- "{"
    if (bold) left <- paste(left , "\\bf " , sep="")
    if (italic) left <- paste(left , "\\it " , sep="")
    
    strv <- paste(left , strv , "}" , sep="")    
  }
  strv
}

## Function to process a line
process.table.row <- function (row , bold , italic , mark , mark.char , format , digits){
  processed <- sapply(1:length(row) , FUN = function(i) print.cell(row[i] , bold[i] , italic[i] , format , digits[i] , mark[i] , mark.char))
  paste(paste(processed, collapse = " & ") , "\\\\")
}

## hrule and vrule, number of cols/rows after which add the rule
## collapse.XX a vector with the rows/cols to collapse. 0==row/file name, if printed.
## scientific puede ser NA
## para formats ver funcion formatC


#' @title Write a table in LaTeX format
#'
#' @description This is a simple function to create tabular environment in LaTeX
#' @param table A data frame with the information to plot
#' @param file Path of a file. If provided, the tabular is wirten in the given file. Otherwise, it is writen to the standard output
#' @param format Format for the numeric values. The accepted formats are those in the function \code{\link{formatC}}. The typical values are \code{'g'} to automatically set theformat, \code{'f'} for a fixed sized floating point format and \code{'e'} or \code{'E'} for scientific notation
#' @param bold A matrix that matches \code{'table'} in size indicating with true those cells that have to be printed in bold font
#' @param italic A matrix that matches \code{'table'} in size indicating with true those cells that have to be printed in bold font
#' @param mark A matrix that matches \code{'table'} in size indicating with true those cells that have to be marked with a superscipt symbol
#' @param mark.char Character to be used to mark cells. Note that the superscript is included in a math environment, so this has to be either a character or a valid math command in LaTeX
#' @param align Character indicating the alignment of the colums (\code{'l'},\code{'r'} or \code{'c'})
#' @param hrule A vector of positions for the horizontal lines in the tabular. All the lines are drawn after the indicated line. When the column names are plotted, 0 means drawing a line after the column names. The maximum value is the number of rows - 1 (for a line after the last line see parametr \code{bty})
#' @param vrule Similar to \code{'hrule'} but for vertical lines. . The maximum value is the number of columns - 1 (for a line after the last columns see parametr \code{bty})
#' @param bty Vector indicating which borders should be printed. The vector can contain any of subset of \code{c('l','r','t','b')}, which represent, respectively, left, right, top and bottom border. If the parameter is set to \code{NULL} no border is printed.
#' @param print.col.names Local value indicating whether the column names have to be printed or not
#' @param print.row.names Local value indicating whether the row names have to be printed or not
#' @param digits A vector with the number of digits in each column. Its size has to match the number of the final table, i.e., the colums in \code{'table'} if the row names are not included or the number of columns + 1 if the row names are printed in the final table
#' @return LaTeX code to print the table
#' @examples
#' data(data.garcia.herrera)
#' ## Fixed width, 4 digits
#' ncol <- dim(data.garcia.herrera)[2]
#' digits <- rep(4,ncol+1)
#' write.tabular(data.garcia.herrera , format = 'f' , hrule = 0 , vrule = 0 , digits = digits)
#' ## Scientific notation, maximum value in each row in bold
#' max.matrix <- t(apply(data.garcia.herrera , MARGIN = 1 , FUN = function(x) x==max(x)))
#' write.tabular(data.garcia.herrera , format = 'E' , bold = max.matrix)

write.tabular <- function (table , file=NULL , format = 'g' , bold=NULL , italic=NULL , mark=NULL, mark.char = '*' , align = 'l' , hrule = NULL , vrule = NULL , bty = c('t','b','l','r') , print.col.names = TRUE , print.row.names = TRUE  , digits = rep(3,ncol(table) + print.row.names)){
  
  rows <- nrow(table)
  cols <- ncol(table)
  
  print.col.names <- print.col.names & !is.null(colnames(table))
  print.row.names <- print.row.names & !is.null(rownames(table))
  
  hmin <- ifelse(print.col.names , 0 , 1)
  vmin <- ifelse(print.row.names , 0 , 1) 
  
  ## Remove any vrule and hrule beyond the limits
  if (!is.null(hrule))
    hrule <- subset(hrule , hrule >= hmin & hrule < rows)
  if (!is.null(vrule))
    vrule <- subset(vrule , vrule >= vmin & vrule < cols)
  
  ## Control of the digits
  if (length(digits) - print.row.names != cols) stop("The number of elements in the digits vector is incorrect. The vector should have length equal to the number of columns in 'table' if 'print.row.names' is false and the number of columns + 1 if 'print.row.names' is true.")
  
  ## Matrices for bold, italic and mark
  if (is.null(bold)) bold <- matrix(rep(FALSE , rows*cols) , ncol = cols)
  if (is.null(italic)) italic <- matrix(rep(FALSE , rows*cols) , ncol = cols)
  if (is.null(mark)) mark <- matrix(rep(FALSE , rows*cols) , ncol = cols)
  
  ## Pass all to character to avoid problems with factors
  table <- apply(table , MARGIN = 2 , FUN = as.character)
  
  ## Include row and/or col names as additional info in the table
  if (print.row.names){
    suppressWarnings(table <- cbind(rownames(table) , table))
    bold <- cbind(rep(FALSE , rows) , bold)
    italic <- cbind(rep(FALSE , rows) , italic)
    mark <- cbind(rep(FALSE , rows) , mark)
    if (!is.null(vrule)) vrule <- vrule + 1
    cols <- cols+1
  }
  
  if (print.col.names){
    suppressWarnings(table <- rbind(colnames(table) , table))
    bold <- rbind(rep(FALSE , cols) , bold)
    italic <- rbind(rep(FALSE , cols) , italic)
    mark <- rbind(rep(FALSE , cols) , mark)
    if (!is.null(hrule)) hrule <- hrule + 1
    rows <- rows + 1 
  }
  
  ## Open the file
  output <- ifelse(is.null(file) , stdout() , file(file , "w"))
  
  ## Begin the tabular
  algn <- ifelse('l' %in% bty, "|" , "")
  if (is.null(vrule)){
    algn <- paste(algn , paste(rep(align , cols) , collapse="") , sep="")
  }else{
    aux <- c(vrule[1] , diff(c(vrule,cols)))
    algn <- paste(algn , paste(c(rep(align,aux[1]) , 
                    unlist(sapply (aux[-1] , FUN = function(x) c("|" , rep(align , x)))))
                  ,collapse = "") , sep = "")
  }
  algn <- paste(algn ,  ifelse('r' %in% bty, "|" , "") , sep = "")
  l <- paste("\\begin{tabular}{" , algn , "}",sep="")
  cat(l, file = output , sep="\n")
  
  if ('t' %in% bty)  cat("\\hline" , file = output , sep="\n")
  
  ## Rows of the table
  current.row<-1
  for (r in 1:rows){
    l <- process.table.row(row = table[r,] , bold = bold[r,] , italic = italic[r,] , mark = mark[r,] , mark.char = mark.char , format = format , digits = digits)
    cat(l, file = output , sep="\n")
    if (current.row %in% hrule) cat("\\hline" , file = output , sep="\n")
    current.row <- current.row+1
  }
    
  if ('b' %in% bty) cat("\\hline", file = output , sep="\n")
  
  ## End of tabular
  cat("\\end{tabular}", file = output , sep="\n")
  ## Bye bye
  if(!is.null(file)) close(output)
}