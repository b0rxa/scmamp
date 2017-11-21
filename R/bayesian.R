#' @title Bayesian equivalent to the correlated t-test 
#'
#' @description Implementation of the Bayesian version of the correlated t-test as presented in Benavoli \emph{et al.} 2017
#' @param x First sample
#' @param y Second sample (if not provided, x is assumed to be the difference)
#' @param rho Correlation factor (see details)
#' @param rope Interval for the difference considered as "irrelevant"
#' @return A list with the following elements: 
#' \item{\code{method}}{a string with the name of the method used}
#' \item{\code{posterior.probabilities}}{a vector with the left, rope and right probabilities}
#' \item{\code{approximated}}{a logical value, \code{TRUE} if the posterior distribution is approximated (sampled) and \code{FALSE} if it is exact}
#' \item{\code{parameters}}{parameters used by the method}
#' \item{\code{posterior}}{posterior density function if the method is exact and the obtained sample if the method is approximated}
#' \item{\code{additional}}{a list that contains the posterior cumulative function (\code{pposterior}), the posterior quantile function (\code{qposterior}) and the parameters of the posterior density function (\code{posterior.df}, \code{posterior.mean}, \code{posterior.sd})}
#' @details Note that the default value for rho is 0, wich converts the test in the equivalent
#' of the standard t-test. To correct due to correlation you need to set the rho
#' parameter. In the case of classifiers compared using any validation scheme the
#' heuristic typically used is to set rho to num. test instances / total num. of instance
#' The function has been implemented to be used in the comparison of classifiers
#' and, in such situation, the heuristic used to fix the correlation factor
#' is the size of the training set divided by the total size of the data.
#' For the same reason, the reference value (the midpoint for the rope) is 0.
#' However, this may be changed (though it matches the prior for the paremter,
#' which follows a Gaussian distribution of 0 mean).
#' The results includes the typicall information relative to the three areas
#' of the posterior density (left, right and rope probabilities), but also
#' the basic functions (density, cummulative and quantile), as well as the
#' parameters of the posterior density function, which is a Student's t.
#' @references A. Benavoli, G. Corani, J. Demsar, M. Zaffalon (2017) Time for a Change: a Tutorial
#' for Comparing Multiple Classifiers Through Bayesian Analysis. \emph{Journal of Machine Learning Research}, 18, 1-36.
#' @examples
#' sample1 <- rnorm(25, 1, 1)
#' sample2 <- rnorm(25, 1.2, 1)
#' correlatedTtest (x=sample1, y=sample2, rho=0.1, alternative="less")
#' 

bCorrelatedTtest <- function (x, y=NULL, rho=0, rope=c(-0.01, 0.01)) {

  # Some checks
  if (rope[2] < rope[1]) {
    warning("The rope paremeter has to contain the ordered limits of the rope
            (min, max), but the values are not orderd. They will be swapped to
            follow with the procedure")
    
    rope <- sort(rope)
  }
  
  if (rho>1) {
    stop("The correlation factor has to be strictly smaller than 1!!")
  }
  
  # Convert the data to differences
  if (!is.null(y)) {
    sample <- x-y
  } else {
    sample <- x
  }
  
  # Get the 
  sample.mean <- mean(sample)
  sample.sd <- sd(sample)
  n <- length(sample)
  
  tdist.df <- n-1
  tdist.mean <- sample.mean
  tdist.sd <- sample.sd*sqrt(1/n + rho/(1-rho))
  
  dpos <- function(mu) {
    #Standarize the value and get the density
    x <- (mu-tdist.mean)/tdist.sd
    return(dt(x,tdist.df))
  }
  
  ppos <- function(mu) {
    #Standarize the value and get the density
    x <- (mu-tdist.mean)/tdist.sd
    return(pt(x,tdist.df))
  }
  
  qpos <- function(q) {
    return(qt(q,tdist.df) * tdist.sd + tdist.mean)
  }
  
  parameters <- list(rope=rope, rho=rho)  
  left.prob <- ppos(rope[1])
  rope.prob <- ppos(rope[2]) - left.prob
  right.prob <- 1 - ppos(rope[2])
  posterior.probs <- c(left.prob, rope.prob, right.prob)
  names(posterior.probs) <- c("Left", "Rope", "Right")
  
  additional <- list(pposterior=ppos, qposterior=qpos, posterior.df=tdist.df, 
                     posterior.mean=tdist.mean, posterior.sd=tdist.sd)
  
  results <- list()
  results$method                  <- "Bayesian correlated t-test"
  results$parameters              <- parameters
  results$posterior.probabilities <- posterior.probs
  results$approximate             <- FALSE
  results$posterior               <- dpos
  results$additional              <- additional
  
  return (results)
}


#' @title Bayesian equivalent to Wilcoxon's signed-rank test 
#'
#' @description Implementation of the Bayesian version of the signed-rank test presented in Benavoli \emph{et al.} 2017
#' @param x First sample
#' @param y Second sample (if not provided, x is assumed to be the difference)
#' @param s Scale parameter of the prior Dirichlet Process. The default value is set to 0.5
#' @param z0 Position of the pseudo-observation associated to the prior Dirichlet Process. The default value is set to 0 (inside the rope)
#' @param rope Interval for the difference considered as "irrelevant"
#' @param nsim Number of samples used to estimate the posterior distribution
#' @param seed Optional parameter used to fix the random seed
#' @return A list with the following elements: 
#' \item{\code{method}}{a string with the name of the method used}
#' \item{\code{posterior.probabilities}}{a vector with the left, rope and right probabilities}
#' \item{\code{approximated}}{a logical value, \code{TRUE} if the posterior distribution is approximated (sampled) and \code{FALSE} if it is exact}
#' \item{\code{parameters}}{parameters used by the method}
#' \item{\code{posterior}}{Sampled probabilities (see details)}
#' @details The results includes the typical information relative to the three 
#' areas of the posterior density (left, right and rope probabilities), but also  
#' the result of the simulation used to estimate the probabilities.
#' 
#' The posterior contains the sampled probabilities in a matrix where each colums corresponds
#' to the sampled (posterior) probability of the measure (z), falling in that particular area (left to the rope, inside the rope and right to the rope). Conversely, the posterior probabilities refer to the probability of each region having the highest probability
#' 
#' As for the prior parameters, they are set to the default values indicated in Benavoli \emph{et al.} 2017 and you should not modify the unless you know what you are doing.
#' @references A. Benavoli, G. Corani, J. Demsar, M. Zaffalon (2017) Time for a Change: a Tutorial
#' for Comparing Multiple Classifiers Through Bayesian Analysis. \emph{Journal of Machine Learning Research}, 18, 1-36.
#' @examples
#' sample1 <- rnorm(25, 1, 1)
#' sample2 <- rnorm(25, 1.2, 1)
#' results <- bSignedRankTest (x=sample1, y=sample2, s=0.5, z0=0)
#' res$posterior.probabilities
#' 

bSignedRankTest <- function(x, y=NULL, s=0.5, z0=0, rope=c(-0.01, 0.01), nsim=1000,
                            seed=as.numeric(Sys.time())) {
  
  if (!require(MCMCpack)) {
    stop("This function requires the MCMCpack package. Please install it and try again.")
  }  
  
  if (rope[2] < rope[1]) {
    warning("The rope paremeter has to contain the ordered limits of the rope
            (min, max), but the values are not orderd. They will be swapped to
            follow with the procedure")
    
    rope <- sort(rope)
  }
  
  # Convert the data to differences
  if (!is.null(y)) {
    sample <- x-y
  } else {
    sample <- x
  }
  
  # Create the parameter vector for the sampling of the weights
  weights.dir.params <- c(s ,rep(1,length(sample)))
  
  # Add the pseudo-observation due to the prior to the sample vector.
  sample <- c (z0, sample)
  
  set.seed(seed)
  weights <- rdirichlet(nsim, weights.dir.params)
  
  
  # Prepare to do get the terms for all the pairs i,j
  aux <- matrix(rep(sample, length(sample)), ncol=length(sample))
  sample.matrix <- aux + t(aux)
  
  # Get the elements to be summed for each event
  left.matrix  <- sample.matrix < 2*rope[1]
  right.matrix <- sample.matrix > 2*rope[2]
  rope.matrix  <- sample.matrix >= 2*rope[1] & sample.matrix <= 2*rope[2]
  
  # To the get the i-th sample create the corresponding weight product matrix
  # and apply the selection defined in the previous matrices
  f <- function(i) {
    aux <- matrix(rep(weights[i, ], length(sample)), ncol=length(sample))
    weights.matrix <- aux * t(aux)
    return(data.frame(Left=sum(left.matrix * weights.matrix),
                      Rope=sum(rope.matrix * weights.matrix),
                      Right=sum(right.matrix * weights.matrix)))
  }
  
  posterior.distribution <- do.call(rbind, lapply(1:nsim, f))
  
  # For the posterior probabilities, we will check, for each region (left, rope, right) 
  # the probability of being the highest
  f2 <- function(s, col){
    max.values <- s==max(s)
    if (max.values[col]) {
      r <- 1/sum(max.values)
    } else { 
      r <- 0
    }
    return(r)
  }
  
  left.prob  <- mean(apply(posterior.distribution, MARGIN=1, FUN=f2, col=1))
  rope.prob  <- mean(apply(posterior.distribution, MARGIN=1, FUN=f2, col=2))
  right.prob <- mean(apply(posterior.distribution, MARGIN=1, FUN=f2, col=3))
  
  posterior.probs <- c(left.prob, rope.prob, right.prob)
  names(posterior.probs) <- c("Left", "Rope", "Right")
  
  
  parameters <- list(rope=rope, s=s, z0=z0)  
  
  results <- list()
  results$method                  <- "Bayesian signed-rank test"
  results$parameters              <- parameters
  results$posterior.probabilities <- posterior.probs
  results$approximate             <- TRUE
  results$posterior               <- posterior.distribution
  
  return (results)
}




#' @title Bayesian hierarchical model for the analysis of two algorithms in multiple datasets 
#'
#' @description Bayesian hierarchical model for the simulatenous analysis of two algorithms in multiple datasets as presented in Benavoli \emph{et al.} 2017
#' @param x.matrix First sample, a matrix with the results obtained by the first algorithm (each dataset in a row)
#' @param y.matrix Second sample, a matrix with the results obtained by the second algorithm (each dataset in a row) (if not provided, x is assumed to be the difference)
#' @param std.upper Factor to set the upper bound for both sigma_i and sigma_0 (see Benavoli \emph{et al.} 2017 for more details)
#' @param d0.lower Lower bound for the prior for mu_0. If not provided, the smallest observed difference is used
#' @param d0.upper Upper bound for the prior for mu_0. If not provided, the biggest observed difference is used
#' @param alpha.lower Lower bound for the (uniform) prior for the alpha hyperparameter (see Benavoli \emph{et al.} 2017 for more details). Default value set at 0.5, as in the original paper
#' @param alpha.upper Upper bound for the (uniform) prior for the alpha hyperparameter (see Benavoli \emph{et al.} 2017 for more details). Default value set at 5, as in the original paper
#' @param beta.lower Lower bound for the (uniform) prior for the beta hyperparameter (see Benavoli \emph{et al.} 2017 for more details). Default value set at 0.05, as in the original paper
#' @param beta.lower Upper bound for the (uniform) prior for the beta hyperparameter (see Benavoli \emph{et al.} 2017 for more details). Default value set at 0.15, as in the original paper
#' @param z0 Position of the pseudo-observation associated to the prior Dirichlet Process. The default value is set to 0 (inside the rope)
#' @param rope Interval for the difference considered as "irrelevant"
#' @param nsim Number of samples (per chain) used to estimate the posterior distribution. Note that, by default, half the simulations are used for the burn-in
#' @param nchain Number of MC chains to be simulated. As half the simulations are used for the warm-up, the total number of simulations will be \code{nchain}*\code{nsim}/2
#' @param parallel Logical value. If \code{true}, Stan code is executed in parallel
#' @param stan.output.file String containing the base name for the output files produced by Stan. If \code{NULL}, no files are stored.
#' @param seed Optional parameter used to fix the random seed
#' @param ... Additional arguments for the rstan::stan function that runs the analysis 
#' @return A list with the following elements: 
#' \item{\code{method}}{a string with the name of the method used}
#' \item{\code{posterior.probabilities}}{a vector with the left, rope and right probabilities}
#' \item{\code{approximated}}{a logical value, \code{TRUE} if the posterior distribution is approximated (sampled) and \code{FALSE} if it is exact}
#' \item{\code{parameters}}{parameters used by the method}
#' \item{\code{posterior}}{Sampled probabilities (see details)}
#' \item{\code{additional}}{Additional information provided by the model. This includes:\code{per.dataset}, the results per dataset (left, rope and right probabilities together with the expected mean value); \code{global.sin} sampled probabilities of mu_0 being positive or negative and \code{stan.results}, the complete set of results produced by Stan program}
#' @details The results includes the typical information relative to the three 
#' areas of the posterior density (left, right and rope probabilities), both global and per dataset (in the additional information). Also, the simulation results are included.
#' 
#' As for the prior parameters, they are set to the default values indicated in Benavoli \emph{et al.} 2017, except for the bound for the prior distribution of mu_0, which are set to the maximum and minimum values observed in the sample. You should not modify them unless you know what you are doing.
#' @references A. Benavoli, G. Corani, J. Demsar, M. Zaffalon (2017) Time for a Change: a Tutorial
#' for Comparing Multiple Classifiers Through Bayesian Analysis. \emph{Journal of Machine Learning Research}, 18, 1-36.
#' @examples
#' sample1 <- matrix(rnorm(25*5, 1, 1), nrow=5)
#' sample2 <- matrix(rnorm(25*5, 1.2, 1), nrow=5)
#' results <- bHierarchicalTest (x.matrix=sample1, y.matrix=sample2, rho=0, rope=c(-0.05, 0.05))
#' res$posterior.probabilities
#' 


bHierarchicalTest <- function(x.matrix, y.matrix=NULL, rho, std.upper=1000, d0.lower=NULL, d0.upper=NULL, 
                              alpha.lower=0.5, alpha.upper=5, beta.lower=0.05, beta.upper=0.15, 
                              rope=c(-0.01, 0.01), nsim=2000, nchains=8, parallel=TRUE, stan.output.file=NULL,
                              seed=as.numeric(Sys.time()), ...) {
  

  
  if (!require(rstan)) {
    stop("This function requires the rstan package. Please install it and try again.")
  }  
  
  if (!require(metRology)) {
    stop("This function requires the ggplot2 package. Please install it and try again.")
  }
  

  if (parallel) {
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())
  }
  
  if (!is.null(stan.output.file)) {
    rstan_options(auto_write = TRUE)
    if (!dir.exists("./stan_out")) {
      dir.create("./stan_out")
    }
    stan.output.file <- paste0("./stan_out/",stan.output.file,".StanOut")
  }
  
  if (rope[2] < rope[1]) {
    warning("The rope paremeter has to contain the ordered limits of the rope
            (min, max), but the values are not orderd. They will be swapped to
            follow with the procedure")
    
    rope <- sort(rope)
  }

  if (is.null(y.matrix)) {
    sample.matrix <- x.matrix
  } else {
    sample.matrix <- x.matrix - y.matrix
  }
    
  # Check the input data (we need a matrix or a data.frame)
  if (class(sample.matrix) == "numeric") {
    sample.matrix <- matrix(sample.matrix, nrow=1)
  }
  
  
  # Code inherited from BayesianTestML/tutoria/hierarchical/hierarchical_test
  # Not sure the reason for this check (in the context of the original code, x
  # values should be bounded to the (-1, 1) interval)
  if ((max(sample.matrix))>1 & rope[2] < 0.02) {
    stop('value of rope_max  not compatible with scale of provided x')
  }
  
  num.samples  <- ncol(sample.matrix)
  num.datasets <- nrow(sample.matrix)
  
  # Scale the problem according to the mean standard deviation of all the datasets
  dataset.sds <- apply(sample.matrix, MARGIN=1, FUN=sd)
  mean.dataset.sd <- mean(dataset.sds)
  
  # Save this value for undoing the scaling in the results
  scale.factor <- mean.dataset.sd
  
  sample.matrix <- sample.matrix / mean.dataset.sd
  
  # Update also the rope
  rope <- rope / mean.dataset.sd
  
  # Also update the limits for d0 in case they are provided
  if (!is.null(d0.lower)) {
    d0.lower <- d0.lower / mean.dataset.sd
  }
  
  if (!is.null(d0.upper)) {
    d0.upper <- d0.upper / mean.dataset.sd
  }
  
  
  # In case there is any dataset with 0 variance, add a small value to avoid problems
  # taking care to not alter the mean value
  for (id in which(dataset.sds==0)) {
    noise <- runif(num.samples/2, rope[1], rope[2])
    sample.matrix[id, 1:(num.samples/2)] <- sample.matrix[id, 1:(num.samples/2)] + noise
    sample.matrix[id, (num.samples/2 + 1):num.samples] <- sample.matrix[id, (num.samples/2 + 1):num.samples] - noise
  }
  
  # Just in case, as the sd of those datasets with sd=0 has changed
  # The mean sd of the dasets is related with the upper bound of the sd
  dataset.sds     <- apply(sample.matrix, MARGIN=1, FUN=sd)
  mean.dataset.sd <- mean(dataset.sds)
  
  # For the upper bound of the sd0, we get the sd of the mean values per
  # dataset. In case we only have one, the upper bound is set using its sd
  if (num.samples==1) {
    dataset.mean.sd <- sd(sample.matrix)
  } else {
    dataset.mean.sd <- sd(apply(sample.matrix, MARGIN=1, FUN=mean))
  }
  
  if (is.null(d0.lower)) {
    d0.lower <- -max(abs(sample.matrix))
  }
  
  if (is.null(d0.upper)) {
    d0.upper <- max(abs(sample.matrix))
  }
  
  data <- list()
  
  data$deltaLow   <- d0.lower
  data$deltaHi    <- d0.upper
  data$stdLow     <- 0
  data$stdHi      <- mean.dataset.sd*std.upper
  data$std0Low    <- 0
  data$std0Hi     <- dataset.mean.sd*std.upper
  data$Nsamples   <- num.samples
  data$q          <- num.datasets 
  data$x          <- sample.matrix 
  data$rho        <- rho
  data$upperAlpha <- alpha.upper
  data$lowerAlpha <- alpha.lower
  data$upperBeta  <- beta.upper
  data$lowerBeta  <- beta.lower
  
  stan.program <- system.file("stan/hierarchical_t_test.stan", package="scmamp")
  
  stan.fit <-  stan(file=stan.program, data=data, iter=nsim, chains=nchains,
                    seed=seed, sample_file=stan.output.file, ...)
  
  stan.results<- extract(stan.fit, permuted = TRUE)
  
  # Remove irrelevant variables
  stan.results$diff<-NULL
  stan.results$diagQuad<-NULL
  stan.results$oneOverSigma2<-NULL
  stan.results$nuMinusOne<-NULL
  stan.results$log_lik<-NULL
  
  
  # Once simulated we now get the relevant probabilities, starting with each dataset
  aux <- apply(stan.results$delta, MARGIN=2, 
               FUN=function(i){
                 res <- data.frame(Left=mean(i<rope[1]),
                                   Right=mean(i>rope[2]),
                                   Rope=mean(i<=rope[2] & i>=rope[1]))
               })
  
  probs.per.dataset <- do.call(rbind, aux)
  
  # Analyse the posterior distribution of delta parameter to check the most probable 
  # outcome for the next delta (future datasets)
  aux <- pt.scaled(rope[2], df=stan.results$nu, mean=stan.results$delta0, sd=stan.results$std0)
  cum.left <- pt.scaled(rope[1], df=stan.results$nu, mean=stan.results$delta0, sd=stan.results$std0)
  cum.rope <- aux - cum.left
  cum.right <- 1 - aux
  
  posterior.distribution <- data.frame("Left"=cum.left, "Rope"=cum.rope, "Right"=cum.right)
  
  left.wins <- cum.left > cum.right & cum.left > cum.rope
  right.wins <- cum.right > cum.left & cum.right > cum.rope
  # Difference with the original code in BayesianTestsML/tutorial. In case of ties
  # the point goes to the rope, in order to be more conservative
  rope.wins <- !(left.wins | right.wins)
  
  positive.d0 <- stan.results$delta0 > 0 
  
  # Get the probabilities according to the counts
  prob.left.win  <- mean(left.wins)
  prob.right.win <- mean(right.wins)
  prob.rope.win  <- mean(rope.wins)
  
  prob.positive  <- mean(positive.d0)
  prob.negative  <- 1 - prob.positive
  
  # Get the results ready
  
  # Remember that the differences had been scaled, so we need to revert that scaling
  per.dataset <- cbind("mean.delta" = colMeans(stan.results$delta)*scale.factor,
                       probs.per.dataset)
  
  global.sign <- c(prob.negative, prob.positive)
  names(global.sign) <- c("Negative", "Positive")
  
  global.wins <- c(prob.left.win, prob.rope.win, prob.right.win)
  names(global.wins) <- c("Left",  "Rope", "Right")
  
  parameters <- list(rho=rho, std.upper=std.upper, d0.lower=d0.lower, d0.upper=d0.upper,
                     alpha.lower=alpha.lower, alpha.upper=alpha.upper, 
                     beta.lower=beta.lower, beta.upper=beta.upper,
                     rope=rope, nsim=nsim, nchains=nchains)
 
  additional <- list(per.dataset=per.dataset, global.sign=global.sign, stan.results=stan.results)
  
  
  results <- list()
  results$method                  <- "Hierarchical Bayesian correlated model"
  results$parameters              <- parameters
  results$posterior.probabilities <- global.wins
  results$approximate             <- TRUE
  results$posterior               <- posterior.distribution
  results$additional              <- additional
  
  
  return(results)
}
