\name{confidence_set}
\alias{confidence_set}
\alias{print.cpt.confset}
\title{
Confidence set of segmentations for a PELT changepoint fit
}
\description{
Computes a confidence set for a changepoint fit produced with \code{method="PELT"}, stores it in the \code{conf.set} slot of the object keyed by the confidence level, and returns the updated object.

Unlike a classical confidence interval, this does not give a range around each individual changepoint.  It computes a set of whole \emph{segmentations} that are plausible at the requested level, so the uncertainty covers both the number of changepoints and their locations at the same time.

\code{\link{confint}} is deliberately identical to this function, it exists so that users can reach the same feature through the name they already know from \code{lm} and \code{glm}.
}
\usage{
confidence_set(object, level = 0.95)
}
\arguments{
  \item{object}{
	An object of class \code{cpt} fitted with \code{method="PELT"}.  Objects of class \code{cpt.reg} and \code{cpt.range} (for example CROPS fits) are not supported.
}
  \item{level}{
	Confidence level, a single number strictly between 0 and 1.  Default is 0.95.  A higher level gives a larger set.
}
}
\details{
	The function re-runs PELT on the data stored in \code{object}, this time asking the C code to also report the near-miss candidate changepoints it would normally prune away, together with their costs.  Every candidate whose cost is within a budget of the optimal cost is kept, and each one is expanded back into a full segmentation.  The budget grows as \code{level} increases, so a higher confidence level accepts more segmentations.

	The optimal segmentation is always a member of the set.

	The cost values come straight from the cost function PELT itself used, they are never recomputed here.  That is what makes the calculation distribution agnostic: \code{Normal} (mean, variance, mean and variance), \code{Exponential}, \code{Gamma} and \code{Poisson} all flow through the same code path, with or without the MBIC penalty.

	WARNING: the budget formula was derived for the Normal case.  The non-Normal distributions will run and return a set, but the coverage of that set has not yet been verified, so treat the level as provisional for those.

	The expensive part of the computation (the PELT re-run and the capture of the raw candidates) does not depend on the level, so it is only ever paid once per fit: the first call parks these raw ingredients in the \code{conf.set} slot under the key \code{"ingredients"}, and every later call at any level (through either \code{confidence_set} or \code{confint}) reuses them, making additional levels essentially free.  If the object is refitted so that the cached ingredients no longer match its changepoints or penalty, they are recomputed automatically.

	Because R functions cannot modify an object in place, the storing only takes effect when the result is assigned back, \code{fit = confidence_set(fit, level=0.95)}.
}
\value{
	The input object with its \code{conf.set} slot updated.  The slot holds one entry per requested level, keyed by the level as a string, plus the cached \code{"ingredients"}.  Levels accumulate: requesting a new level keeps the earlier ones.

	Each level entry is a list of class \code{"cpt.confset"} with elements:
	\item{level}{The confidence level used.}
	\item{cpts}{The optimal changepoints, ending with n.}
	\item{segmentations}{A list of the plausible segmentations, named by how many changepoints each contains (n counts as one).  Each element is a matrix with one segmentation per row.  The optimal segmentation is included.}

	Retrieve an entry with \code{conf.set(fit)[["0.95"]]}.  A \code{print} method shows the optimal changepoints first and then the alternatives.
}
\references{
PELT Algorithm: Killick R, Fearnhead P, Eckley IA (2012) Optimal detection of changepoints with a linear computational cost, \emph{JASA} \bold{107(500)}, 1590--1598
}
\author{
Arslan Batyrov
}
\seealso{
\code{\link{confint}},\code{\link{conf.set}},\code{\link{cpt.mean}},\code{\linkS4class{cpt}}
}
\examples{
set.seed(1)
x = c(rnorm(50,0,1),rnorm(50,5,1),rnorm(50,0,1))
fit = cpt.mean(x, method="PELT")

# compute and store the set, then read it back out of the slot
fit = confidence_set(fit, level=0.95)
conf.set(fit)[["0.95"]]

# a second level reuses the cached ingredients, no PELT re-run.
# confint is an identical twin of confidence_set, both levels now live on the fit
fit = confint(fit, level=0.99)
names(conf.set(fit))
}
\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
