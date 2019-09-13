######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# tcl_hessian
#
# Part of R/tlc - Testing in Conditional Likelihood Context package
#
# This file contains a routine that computes the Hessian matrix
# using function eRm_cml
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2019, Last Modified 06/09/2019
######################################################################
#' Computation of Hessian matrix.
#'
#' Uses function \code{hessian()} from numDeriv package to compute (approximate numerically) Hessian matrix
#'   evaluated at arbitrary values of item easiness parameters.
#'
#' @param X data matrix.
#' @param eta numeric vector of item easiness parameters.
#' @param W design matrix.
#' @param model RM, PCM, RSM, LLTM.
#' @return Hessian matrix evaluated at eta
#'
#' @references{
#' Gilbert, P., Gilbert, M. P., & Varadhan, R. (2016). numDeriv: Accurate Numerical Derivatives. R package
#' version 2016.8-1.1. url: https://CRAN.R-project.org/package=numDeriv
#'  }
#' @keywords htest
#' @export
#' @examples
#' # Rasch model with beta_1 restricted to 0
#' y <- eRm::raschdat1
#' res <- eRm::RM(X = y, sum0 = FALSE)
#' mat <- tcl_hessian(X = y, eta = res$etapar, model = "RM")

tcl_hessian <-function( X, eta, W, model = "RM") {
# X = observed data matrix
# eta - numeric vector of item easiness parameters
# W = design matrix
# model = RM, PCM, RSM, LLTM

  if (!is.numeric(eta)) stop("eta needs to be a numeric vector!")

  e <- eRm_cml(X = X, eta = eta, W = W, model = model)$hessian

  return(e)

}
