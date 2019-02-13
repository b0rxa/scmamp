# NON-EXPORTED, AUXILIAR FUNCTIONS --------------------------------------------

printCell <- function (value, bold, italic, format, digits, mark, mark.char, na.as, ...){
  # Auxiliar function to write a single cell of the table
  # Args:
  #   value:     Value to be printed
  #   bold:      Logic value indicating whether the value is in bold font
  #   italic:    Logic value indicating whether the value is in italic font
  #   format:    String indicating the format of the value, in case it is numeric
  #   digits:    Number of digits for numeric values
  #   mark:      Logic value indicating whether the value has to be marked
  #   mark.char: Mark to be used (any LaTeX Math valid symbol)
  #   na.as:     How NA have to be printed
  #
  # Returns:
  #   The text correspondig to the cell
  #
  value <- unlist(value)
  # Get the string to print. Convert it to either character or numeric
  if (is.factor(value)) {
    value <- as.character(value)
  }
  suppressWarnings(expr={
                     value.num <- as.numeric(value)
                   })
  if (!is.na(value.num)) {
    strv  <- formatC(value.num, digits=digits, format=format)
  }else{
    if (is.na(value)) {
      strv <- na.as
    } else {
      strv <- value
    }
  }
  
  # Check mark
  if (mark){
    strv <- paste(strv, "$^", mark.char,"$", sep="")
  }
  
  # Check bold and/or italic
  if (bold | italic){
    left <- "{"
    if (bold) {
      left <- paste(left, "\\bf ", sep="")
    }
    if (italic) {
      left <- paste(left, "\\it ", sep="")
    }
    
    strv <- paste(left, strv, "}", sep="")    
  }
  
  return(strv)
}



## Function to process a line
processTableRow <- function (row, bold, italic, format, digits, 
                             mark, mark.char, na.as) {
  # Auxiliar function to write a table row
  # Args:
  #   row:       Row to be processed
  #   bold:      Logic value indicating whether the value is in bold font
  #   italic:    Logic value indicating whether the value is in italic font
  #   format:    String indicating the format of the value, in case it is numeric
  #   digits:    Number of digits for numeric values
  #   mark:      Logic value indicating whether the value has to be marked
  #   mark.char: Mark to be used (any LaTeX Math valid symbol)
  #   na.as:     How NA have to be printed
  #
  # Returns:
  #   The text correspondig to the cell
  #
  processed <- sapply(1:length(row), 
                      FUN = function(i) {
                        printCell(value=row[i], bold=bold[i], italic=italic[i],
                                  format=format, digits=digits[i], mark=mark[i],
                                  mark.char=mark.char, na.as=na.as)
                        })
  line <- paste(paste(processed, collapse=" & "), "\\\\")
  return(line)
}




# EXPORTED FUNCTIONS -----------------------------------------------------------

#' @title Write a table in LaTeX format
#'
#' @description This is a simple function to create tabular environment in LaTeX
#' @param table A data frame with the information to write
#' @param file Path of a file. If provided, the tabular is wirten in the given file. Otherwise, it is writen to the standard output
#' @param format Format for the numeric values. The accepted formats are those in the function \code{\link{formatC}}. The typical values are \code{'g'} to automatically set the format, \code{'f'} for a fixed sized floating point format and \code{'e'} or \code{'E'} for scientific notation
#' @param bold A matrix that matches \code{'table'} in size indicating with \code{TRUE} those cells that have to be printed in bold font
#' @param italic A matrix that matches \code{'table'} in size indicating with \code{TRUE} those cells that have to be printed in italic
#' @param mark A matrix that matches \code{'table'} in size indicating with \code{TRUE} those cells that have to be marked with a superscipt symbol
#' @param mark.char Character to be used to mark cells. Note that the superscript is included in a math environment, so this has to be either a character or a valid math command in LaTeX
#' @param na.as Character to be used to write NA values in the table
#' @param align Character indicating the alignment of the colums (\code{'l'},\code{'r'} or \code{'c'})
#' @param hrule A vector of positions for the horizontal lines in the tabular. All the lines are drawn after the indicated line. When the column names are included, 0 means drawing a line after the column names. The maximum value is the number of rows - 1 (for a line after the last line see parametr \code{bty})
#' @param vrule Similar to \code{'hrule'} but for vertical lines. . The maximum value is the number of columns - 1 (for a line after the last columns see parametr \code{bty})
#' @param bty Vector indicating which borders should be printed. The vector can contain any of subset of \code{c('l','r','t','b')}, which represent, respectively, left, right, top and bottom border. If the parameter is set to \code{NULL} no border is printed.
#' @param print.col.names Logical value indicating whether the column names have to be printed or not
#' @param print.row.names Logical value indicating whether the row names have to be printed or not
#' @param digits A single number or a numeric vector with the number of digits in each column. Its size has to match the number of the final table, i.e., the colums in \code{'table'} if the row names are not included or the number of columns + 1 if the row names are printed in the final table
#' @param wrap.as.table Logical value indicating whether the latex object has to be wrapped into a table enviroment
#' @param table.position Character indicating the position of the table (\code{'h'}: here, \code{'t'}: top, or \code{'b'}: botton)
#' @param caption Character string containing the caption of the table. If NULL, no caption is printed
#' @param caption.position Character indicating the possition of the caption (\code{t}: top, the caption is printed over the table; \code{b}: botton, the caption is printed under the table)
#' @param centering Logical value indicating whether the table should be centered in the page
#' @param label Character string containing the label of the table for further references. If NULL, no label is used
#' @return LaTeX code to print the table
#' @seealso \code{\link{summarizeData}}, \code{\link{filterData}} and the vignette \code{vignette(topic="Data_loading_and_manipulation", 
#' package="scmamp")}
#' @examples
#' data(data_blum_2015)
#' args <- list()
#' # Write the summarization of the data
#' args$table <- summarizeData(data.blum.2015, group.by=1:2)
#' 
#' # Set in bold the maximum values per row
#' bold <- apply(args$table[, -(1:2)], MARGIN=1, 
#'              FUN=function(x) {
#'                return(x==max(x))
#'              })
#' args$bold <- cbind(FALSE, FALSE, t(bold))
#' # Fixed width, 2 decimals for the values, 0 for the size and 3 for the radius
#' args$format <- "f"
#' args$digits <- c(0,3,rep(2, 8))
#' 
#' # Print the colum names but not the row names
#' args$print.row.names <- FALSE
#' 
#' # Only top and bottom borders
#' args$bty <- c("t","b")
#' 
#' # Add additional horizontal rules to separate the sizes
#' args$hrule <- c(0,10,20,30)
#' 
#' # An additional vertical rule to separate size and radius from the results
#' args$vrule <- 2
#' 
#' # Print the table
#' do.call(writeTabular, args)

writeTabular <- function (table, file=NULL, format="g", bold=NULL, italic=NULL, 
                          mark=NULL, mark.char="*", na.as="n/a", align="l", 
                          hrule=NULL, vrule=NULL, bty=c("t","b","l","r"), 
                          print.col.names=TRUE, print.row.names=TRUE, digits=3,
                          wrap.as.table=FALSE, table.position="h", caption=NULL,
                          caption.position="b", centering=FALSE, label=NULL) {
  
  rows <- nrow(table)
  cols <- ncol(table)
  
  print.col.names <- print.col.names & !is.null(colnames(table))
  print.row.names <- print.row.names & !is.null(rownames(table))
  
  if (length(digits) == 1) {
    digits <- rep(digits, ncol(table) + print.row.names)
  }
  
  # Control of the digits
  if (length(digits) - print.row.names != cols) {
    stop("The number of elements in the digits vector is incorrect. The vector ",
         "should have length equal to the number of columns in 'table' if ",
         "'print.row.names' is false and the number of columns + 1 if ",
         "'print.row.names' is true. Alternatively, it can have a single value ",
         "to be used in all the numeric columns")
  }
  
  if (print.col.names) {
    hmin <- 0
  } else {
    hmin <- 1
  }
  if (print.row.names) {
    vmin <- 0
  } else {
    vmin <- 1
  }
  
  # Remove any vrule and hrule beyond the limits
  if (!is.null(hrule)) {
    hrule <- subset(hrule, 
                    subset=((hrule >= hmin) & (hrule < rows)))
  }
  if (!is.null(vrule)){
    vrule <- subset(vrule, 
                    subset=((vrule >= vmin) & (vrule < cols)))
  }
  
  # In case being NULL, build the bold, italic and mark matrices
  if (is.null(bold)) {
    bold <- matrix(rep(FALSE, rows * cols), ncol=cols)
  }
  if (is.null(italic)) {
    italic <- matrix(rep(FALSE, rows * cols), ncol=cols)
  }
  if (is.null(mark)) {
    mark <- matrix(rep(FALSE, rows * cols), ncol=cols)
  }
  
  # Pass all to character to avoid problems with factors
  nm <- rownames(table)
  table <- apply(table, MARGIN=2, FUN=as.character)
  ## if table was a matrix with only one row, apply will transform it in a vector and it may
  ## cause problems letter. Therefore, it is converted back to a matrix
  if(is.null(dim(table))){
    table <- t(as.matrix(table))
  }
  rownames(table) <- nm
  
  # Include row and/or col names as additional info in the table
  if (print.row.names) {
    suppressWarnings(expr={
      table <- cbind(rownames(table), table)
    })
    bold <- cbind(rep(FALSE, rows), bold)
    italic <- cbind(rep(FALSE, rows), italic)
    mark <- cbind(rep(FALSE, rows), mark)
    if (!is.null(vrule)) {
      vrule <- vrule + 1
    }
    cols <- cols + 1
  }
  
  if (print.col.names) {
    suppressWarnings(expr={
      table <- rbind(colnames(table) , table)
    })
    bold <- rbind(rep(FALSE, cols), bold)
    italic <- rbind(rep(FALSE, cols), italic)
    mark <- rbind(rep(FALSE, cols), mark)
    if (!is.null(hrule)) {
      hrule <- hrule + 1
    }
    rows <- rows + 1 
  }
  
  # Open the file or the stdout
  if(is.null(file)) {
    out.file <- stdout()
  } else {
    out.file <- file(file, "w")
  }
  
  # begin table enviroment if wrap.as.table
  if(wrap.as.table){
    l <- paste0("\\begin{table}[", table.position, "]")
    cat(l, file=out.file, sep="\n")
    # caption on top
    if(caption.position=="t" & !is.null(caption)){
      l <- paste0("\\caption{", caption, "}")
      cat(l, file=out.file, sep="\n")
    }
    # centering
    if(centering){
      cat("\\centering", file=out.file, sep="\n")
    }
  }
  
  
  # Begin the tabular
  if ("l" %in% bty) {
    algn <- "|"
  } else {
    algn <- ""
  }
  
  if (is.null(vrule) | length(vrule) == 0) {
    algn <- paste(algn, paste(rep(align, cols), collapse=""), sep="")
  } else {
    aux <- c(vrule[1], diff(c(vrule, cols)))
    aux2 <- sapply (aux[-1], 
                    FUN=function(x) {
                      c("|", rep(align, x))
                    })
    algn <- paste(algn, paste(c(rep(align, aux[1]), unlist(aux2)), collapse=""),
                  sep="")
  }
  
  if ("r" %in% bty) {
    algn.end <- "|"  
  } else {
    algn.end <- ""
  }
  algn <- paste(algn, algn.end, sep = "")
  
  l <- paste("\\begin{tabular}{", algn, "}", sep="")
  cat(l, file=out.file, sep="\n")
  
  if ('t' %in% bty)  {
    cat("\\hline", file=out.file, sep="\n")
  }
  
  # Rows of the table
  current.row <- 1
  for (r in 1:rows) {
    l <- processTableRow(row=table[r, ], bold=bold[r, ], italic=italic[r, ], 
                         mark=mark[r,], mark.char=mark.char, format=format, 
                         digits=digits, na.as=na.as)
    cat(l, file=out.file, sep="\n")
    if (current.row %in% hrule) {
      cat("\\hline", file=out.file, sep="\n")
    }
    current.row <- current.row+1
  }
  
  if ('b' %in% bty) {
    cat("\\hline", file=out.file, sep="\n")
  }
  
  # End of tabular
  cat("\\end{tabular}", file=out.file, sep="\n")
  
  # end table enviroment if wrap.as.table
  if(wrap.as.table){
    # caption on botton
    if(caption.position=="b" & !is.null(caption)){
      l <- paste0("\\caption{", caption, "}")
      cat(l, file=out.file, sep="\n")
    }
    # label
    if(!is.null(label)){
      l <- paste0("\\label{", label, "}")
      cat(l, file=out.file, sep="\n")
    }
    # end of table
    cat("\\end{table}", file=out.file, sep="\n")
  }
  
  # Bye bye
  if(!is.null(file)) {
    close(out.file)
  }
}
