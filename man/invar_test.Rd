% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/invar_test.R
\name{invar_test}
\alias{invar_test}
\title{Test of invariance of item parameters between two groups.}
\usage{
invar_test(X, splitcr = "median", model = "RM")
}
\arguments{
\item{X}{data matrix.}

\item{splitcr}{split criterion which is either "mean", "median" or a numeric vector x.
\describe{
 \item{"mean"}{corresponds to division of the sample according to the mean of the person score.}
 \item{"median"}{corresponds to division of the sample according to the median of the person score.}
  \item{x}{has length equal to number of persons and contains zeros and ones indicating group membership of the persons.}
  }}

\item{model}{RM, PCM, RSM}
}
\value{
A list of test statistics, degrees of freedom, and p-values.
 \item{test}{a numeric vector of gradient (GR), likelihood ratio (LR), Rao score (RS), and Wald test statistics.}
 \item{df}{degrees of freedom.}
 \item{pvalue}{a numeric vector of corresponding p-values.}
 \item{call}{the matched call.}
}
\description{
Computes gradient (GR), likelihood ratio (LR), Rao score (RS) and Wald (W) test statistics
  for hypothesis of equality of item parameters between two groups of persons against a two-sided
 alternative that at least one item parameter differs between the two respected groups.
}
\examples{
##### Rasch Model #####
y <- eRm::sim.rasch(persons = rnorm(400), c(0,-3,-2,-1,0,1,2,3))
x <- c(rep(1,200),rep(0,200))

res <- invar_test(y, splitcr = x, model = "RM")

res$test # test statistics
res$df # degrees of freedoms
res$pvalue # p-values

}
\references{
{
Draxler, C. (2010). Sample Size Determination for Rasch Model Tests. Psychometrika, 75(4), 708–724.

Draxler, C., & Alexandrowicz, R. W. (2015). Sample Size Determination Within the Scope of Conditional Maximum Likelihood Estimation
with Special Focus on Testing the Rasch Model. Psychometrika, 80(4), 897–919.

Draxler, C., Kurz, A., & Lemonte, A. J. (2019). The Gradient Test and its Finite Sample Size Properties in a Conditional Maximum Likelihood
and Psychometric Modeling Context. Submitted for publication.

Glas, C. A. W., & Verhelst, N. D. (1995a). Testing the Rasch Model. In G. H. Fischer & I. W. Molenaar (Eds.),
Rasch Models: Foundations, Recent Developments, and Applications (pp. 69–95). New York: Springer.

Glas, C. A. W., & Verhelst, N. D. (1995b). Tests of Fit for Polytomous Rasch Models. In G. H. Fischer & I. W. Molenaar (Eds.),
Rasch Models: Foundations, Recent Developments, and Applications (pp. 325-352). New York: Springer.

Lemonte, A. J. (2016). The Gradient Test. Another Likelihood-Based Test. London:Academic Press.

Terrell, G. R. (2002). The Gradient Statistic. Computing Science and Statistics, 34(34), 206–215.
 }
}
\seealso{
\code{\link{change_test}}, and \code{\link{LLTM_test}}.
}
\keyword{htest}
