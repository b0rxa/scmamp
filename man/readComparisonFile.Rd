% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loading_functions.R
\name{readComparisonFile}
\alias{readComparisonFile}
\title{Read data from a comparison file}
\usage{
readComparisonFile(file, alg.cols, col.names = NULL, ...)
}
\arguments{
\item{file}{Path of the file to load}

\item{alg.cols}{A vector column names or indices inicating which columns contain the results. The rest are assumed as descriptors of the problems}

\item{col.names}{Vector of names of the columns. If not NULL, the files are assumed not to have a header and the columns are named using this vector}

\item{...}{Additional parameters for the read.csv function used to load the data. It can be used, for example, to set the separator (e.g., \code{sep="\t"}). Note that the \code{header} argument is automatically set according to the \code{col.names} argument.}
}
\value{
A data.frame where each column represents either a feature of the experiment or the result of running an algorithm. Algorithm columns are placed always at the end of the table.
}
\description{
This function reads the data from a files where two or more algorithms are compared in different problems. The file can have some columns that characterize the problem and one column per algorithm. If each row contain only the result obtained by one algorithm, use the \code{\link{readExperimentFile}} function.
}
\examples{
dir <- system.file("loading_tests",package="scmamp")
file <- paste(dir , "rgg_complete_comparison.out" , sep="/")
data <- readComparisonFile(file=file, alg.cols=3:10)
dim(data)
head(data)

}
\seealso{
\code{\link{readExperimentFile}}, \code{\link{readExperimentDir}}, \code{\link{readComparisonDir}} and the vignette \code{vignette(topic="Data_loading_and_manipulation", package="scmamp")}
}
