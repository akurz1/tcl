######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# change_test
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothesis that a constant shift parameter equals zero,
# where one assumes all items presented twice to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 09/09/2019
######################################################################
#' Tests in context of measurement of change using LLTM.
#'
#'Computes gradient (GR), likelihood ratio (LR), Rao score (RS) and Wald (W) test statistics
#'  for hypotheses on parameters expressing change between two time points.
#'
#' Assume all items be presented twice (2 time points) to the same persons.
#'   The data matrix X has n rows (number of persons) and 2k columns considered as virtual items.
#'   Assume a constant shift of item difficulties of each item between the 2 time points represented
#'   by one parameter. The shift parameter is the only parameter of interest.
#'   Test of hypothesis that shift parameter equals zero against the two-sided alternative that shift parameter
#'   is not equal to zero.
#'
#' @param X data matrix containing the responses of n persons to 2k binary items.
#'   Columns 1 to k contain the responses to k items at time point 1,
#'   and columns (k+1) to 2k the responses to the same k items at time point 2.
#' @return A list of test statistics, degrees of freedom, and p-values.
#'  \item{test}{a numeric vector of gradient (GR), likelihood ratio (LR), Rao score (RS), and Wald test statistics.}
#'  \item{df}{degrees of freedom.}
#'  \item{pvalue}{a numeric vector of corresponding p-values.}
#'  \item{call}{the matched call.}
#' @references{
#'  Fischer, G. H. (1995). The Linear Logistic Test Model. In G. H. Fischer & I. W. Molenaar (Eds.),
#'  Rasch models: Foundations, Recent Developments, and Applications (pp. 131-155). New York: Springer.
#'
#'  Fischer, G. H. (1983). Logistic Latent Trait Models with Linear Constraints. Psychometrika, 48(1), 3-26.
#'  }
#' @keywords htest
#' @export
#' @seealso \code{\link{invar_test}}, and \code{\link{LLTM_test}}.
#' @examples
#' # Numerical example with 400 persons and 4 items
#' # presented twice, thus 8 virtual items
#'
#' # Data y generated under the assumption that shift parameter equals 0
#' # (no change from time point 1 to 2)
#'
#' # design matrix W used only for exmaple data generation
#' #     (not used for estimating in change_test function)
#' W <- rbind(c(1,0,0,0,0),
#'   c(0,1,0,0,0),
#'   c(0,0,1,0,0),
#'   c(0,0,0,1,0),
#'   c(1,0,0,0,1),
#'   c(0,1,0,0,1),
#'   c(0,0,1,0,1),
#'   c(0,0,0,1,1))
#'
#' # eta Parameter, first 4 are nuisance
#' # (easiness parameters of the 4 items at time point 1),
#' # last one is the shift parameter
#' eta <- c(-2,-1,1,2,0)
#'
#' y <- eRm::sim.rasch(persons = rnorm(400), items = colSums(eta * t(W)))
#'
#' res <- change_test(X = y)
#'
#' res$test # test statistics
#' res$df # degrees of freedoms
#' res$pvalue # p-values
#'
# #' #' @importFrom eRm LLTM sim.rasch

change_test <- function(X) {

  # X = observed data matrix comprised of 2k columns (twice the number of times k)
  # representing 2k so-called virtual items

  call <- match.call()

  model <- "LLTM"  # model used in computing scorefunction and Hessian matrix

  ###################################################
  #               main programm
  ###################################################
  # design matrix with (k+1) columns: k eta parameter + shift parameter
  # with assumption shift par. equals 0 (no change from time point 1 to 2)
  y <- X
  k <- dim(y)[2] / 2  # number of items

  r  <- eRm::LLTM(y, mpoints = 2,se = TRUE, sum0 = FALSE)

  W1 <- r$W[,1:(k - 1)]                       # design matrix W1 restricting shift parameter to 0

  r1 <- eRm::LLTM( y , W = W1, sum0 = FALSE)


  ###################################################
  ##        compute  GR, LR, RS and Wald test      ##
  ###################################################

  eta.unrest <- r$etapar         # unrestricted CML estimates of eta Parameters
  eta.rest <- c(r1$etapar, 0)    # restricted CML estimates of eta parameters

  ## restricted scorefun and hessian
  rest.1 <- eRm_cml(X = y, eta = eta.rest, W = r$W, model = model)

  LR <- (r1$loglik -  r$loglik) * (-2) # LR test

  RS <- sum ( colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun ) # score test

  GR <- sum(rest.1$scorefun * -eta.unrest)    # gradient test statistic

  Wald <- (eta.unrest[k] / r$se.eta[k])^2 # Wald test

  test.stats <- c( GR, LR, RS, Wald)
  names(test.stats) <- c("GR", "LR", "RS", "W")

  df <- 1
  pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))

  res.list <- list("test" = round(test.stats, digits = 3),
                   "df" = df,
                   "pvalue" = round(pvalue, digits = 3),
                   "call" = call)

  return(res.list)
}
