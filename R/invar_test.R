######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# invar_test
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes Wald (W), Score (RS),
# likelihood ratio (LR) and gradient (GR) test
# for hypothesis that the difference of two group parameters equals zero,
# where one assumes all items presented to all persons.
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 10/09/2019
######################################################################
#' Test of invariance of item parameters between two groups.
#'
#' Computes gradient (GR), likelihood ratio (LR), Rao score (RS) and Wald (W) test statistics
#'   for hypothesis of equality of item parameters between two groups of persons against a two-sided
#'  alternative that at least one item parameter differs between the two respected groups.
#'
#' @param X data matrix.
#' @param splitcr split criterion which is either "mean", "median" or a numeric vector x.
#' \describe{
#'  \item{"mean"}{corresponds to division of the sample according to the mean of the person score.}
#'  \item{"median"}{corresponds to division of the sample according to the median of the person score.}
#'   \item{x}{has length equal to number of persons and contains zeros and ones indicating group membership of the persons.}
#'   }
#'
#' @param model RM, PCM, RSM
#' @return A list of test statistics, degrees of freedom, and p-values.
#'  \item{test}{a numeric vector of gradient (GR), likelihood ratio (LR), Rao score (RS), and Wald test statistics.}
#'  \item{df}{degrees of freedom.}
#'  \item{pvalue}{a numeric vector of corresponding p-values.}
#'  \item{call}{the matched call.}
#'
#' @references{
#' Draxler, C. (2010). Sample Size Determination for Rasch Model Tests. Psychometrika, 75(4), 708–724.
#'
#' Draxler, C., & Alexandrowicz, R. W. (2015). Sample Size Determination Within the Scope of Conditional Maximum Likelihood Estimation
#' with Special Focus on Testing the Rasch Model. Psychometrika, 80(4), 897–919.
#'
#' Draxler, C., Kurz, A., & Lemonte, A. J. (2019). The Gradient Test and its Finite Sample Size Properties in a Conditional Maximum Likelihood
#' and Psychometric Modeling Context. Submitted for publication.
#'
#' Glas, C. A. W., & Verhelst, N. D. (1995a). Testing the Rasch Model. In G. H. Fischer & I. W. Molenaar (Eds.),
#' Rasch Models: Foundations, Recent Developments, and Applications (pp. 69–95). New York: Springer.
#'
#' Glas, C. A. W., & Verhelst, N. D. (1995b). Tests of Fit for Polytomous Rasch Models. In G. H. Fischer & I. W. Molenaar (Eds.),
#' Rasch Models: Foundations, Recent Developments, and Applications (pp. 325-352). New York: Springer.
#'
#' Lemonte, A. J. (2016). The Gradient Test. Another Likelihood-Based Test. London:Academic Press.
#'
#' Terrell, G. R. (2002). The Gradient Statistic. Computing Science and Statistics, 34(34), 206–215.
#'  }
#' @keywords htest
#' @export
#' @seealso \code{\link{change_test}}, and \code{\link{LLTM_test}}.
#' @examples
#'##### Rasch Model #####
#'y <- eRm::sim.rasch(persons = rnorm(400), c(0,-3,-2,-1,0,1,2,3))
#'x <- c(rep(1,200),rep(0,200))
#'
#'res <- invar_test(y, splitcr = x, model = "RM")
#'
#'res$test # test statistics
#'res$df # degrees of freedoms
#'res$pvalue # p-values
#'
# #' @importFrom eRm RM PCM RSM sim.rasch

invar_test <- function(X, splitcr = "median", model = "RM"){
  # X = observed data matrix comprised of k columns
  # splitcr... splitting criterion for 2 covariate groups.
  #   "median"  corresponds to a median person score split,
  #   "mean" corresponds to the mean person score split.
  #    vector of length n containing zeros or ones only for sample split
  #        (group 1 = '1', group 2 = '0')
  # model = RM, PCM, RSM

  call<-match.call()

  ###################################################
  #               main programm
  ###################################################

  # restricted CML estimates
  if (model == "RM")  {
    r  <- eRm::RM( X, se = FALSE, sum0 = FALSE)   # total sample
  } else if (model == "PCM")  {
    r  <- eRm::PCM( X, se = FALSE, sum0 = FALSE)   # total sample
  } else if (model == "RSM")  {
    r  <- eRm::RSM( X, se = FALSE, sum0 = FALSE)   # total sample
  } # end if

  ###################################################
  ######         compute eRm LR test               ##
  ###################################################
  # object =  Object of class "Rm"

  e <- tcl_LRtest.Rm(object = r, splitcr = splitcr, se = FALSE, sum0 = FALSE)
  #### LRtest.Rm.tcl

  LR <- e$LR

  df <- e$df
  pvalue.LR <- e$pvalue

  ###################################################
  #####         computes GR, RS and Wald test         ##
  ###################################################
  if ( is.numeric(splitcr) ) {
    y1 <- e$X.list$`1`   # data group 1
    y2 <- e$X.list$`0`   # data group 2

    ## CML estimates of eta parameters
    eta.rest     <- r$etapar      # CML eta estimates from total sample
    eta.unrest.1 <- e$etalist$`1` # CML eta estimates group 1
    eta.unrest.2 <- e$etalist$`0` # CML eta estimates group 2

  } else if ( !is.numeric(splitcr) && (splitcr=="median" || splitcr=="mean")  ) {
    y1 <- e$X.list$low    # data group 1 = low
    y2 <- e$X.list$high   # data group 2 = high

    ## CML estimates of eta parameters
    eta.rest     <- r$etapar         # CML eta estimates from total sample
    eta.unrest.1 <- e$etalist$low    # CML eta estimates group 1 low person score
    eta.unrest.2 <- e$etalist$high   # CML eta estimates group 2 high person score
  }

  ####################################
  ####      test statistics       ####
  ####################################

  # score function and Hessian matrix evaluated at unrestricted eta estimates
  unrest.1 <- eRm_cml(X = y1, eta = eta.unrest.1, model = model)
  unrest.2 <- eRm_cml(X = y2, eta = eta.unrest.2, model = model)

  # score function and Hessian matrix evaluated at restricted eta estimates
  rest.1 <- eRm_cml(X = y1, eta = eta.rest, model = model)
  rest.2 <- eRm_cml(X = y2, eta = eta.rest, model = model)

  # Fisher informaion matrix
  i <- solve( solve(unrest.1$hessian ) + solve(unrest.2$hessian ))

  # difference of eta parameters in groups
  delta <- eta.unrest.1 - eta.unrest.2

  Wald <- sum( colSums(delta * i) * delta) # classical Wald test

  GR <- sum( (rest.1$scorefun * delta) )  # gradient test

  RS <- sum (colSums( (rest.1$scorefun * solve(rest.1$hessian) )) * rest.1$scorefun) +
    sum (colSums( (rest.2$scorefun * solve(rest.2$hessian) )) * rest.2$scorefun)

  #  } # end if

  test.stats <- c( GR, LR, RS, Wald)
  names(test.stats) <- c("GR", "LR", "RS", "W")

  pvalue <- 1 - (sapply(test.stats, stats::pchisq, df = df))


  res.list <- list("test" = round(test.stats, digits = 3),
                   "df" = df,
                   "pvalue" = round(pvalue, digits = 3),
                   "call" = call)

  return(res.list)
}
