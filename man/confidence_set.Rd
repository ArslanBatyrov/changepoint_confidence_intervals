\name{confidence_set}
\alias{confidence_set}
\alias{print.cpt.confset}
\alias{print.cpt.confsets}
\title{
Confidence set of segmentations for a PELT changepoint fit
}
\description{
Finds which segmentations, other than the best one, are still plausible at a given confidence level.

Rather than a range around each changepoint, this gives a set of whole \emph{segmentations}, so it covers uncertainty in how many changepoints there are as well as where they sit.

The set is stored on the fit and the fit is returned, so remember to assign it back:
\code{fit = confidence_set(fit, level=0.95)}.
}
\usage{
confidence_set(object, level = 0.95)

\method{print}{cpt.confset}(x, n = getOption("changepoint.confset.n", 20), ...)
\method{print}{cpt.confsets}(x, n = getOption("changepoint.confset.n", 20), ...)
}
\arguments{
  \item{object}{
	A \code{cpt} object fitted with \code{method="PELT"}.  BinSeg, SegNeigh and CROPS fits are not supported yet.
}
  \item{level}{
	Confidence level, a single number between 0 and 1.  Default 0.95.  A higher level gives a bigger set.
}
  \item{x}{
	A stored confidence set, or the whole \code{conf.set} slot.
}
  \item{n}{
	How many alternative segmentations to print.  Default 20, use \code{n=Inf} to see them all, or set \code{options(changepoint.confset.n=)} to change the default.  This only affects printing, the full set is always kept.
}
  \item{\dots}{
	Ignored.
}
}
\details{
	PELT throws away the near-miss changepoints it prunes.  This function re-runs it asking for those to be kept, then keeps every candidate whose cost is within a budget of the best one and rebuilds each into a full segmentation.  A higher \code{level} means a bigger budget, so more segmentations qualify.  The best segmentation is always in the set.

	Costs come straight from the cost function PELT used and are never recomputed, so Normal (mean, variance, mean and variance), Exponential, Gamma and Poisson all work, with or without the MBIC penalty.

	WARNING: the budget formula was derived for the Normal case.  The other distributions run and return a set, but their coverage has not been checked yet, so treat the level as provisional there.

	The slow part (the PELT re-run) does not depend on the level, so it is only done once per fit and cached in the slot.  Extra levels after the first are nearly free, whether you ask through \code{confidence_set} or \code{\link{confint}}.  The cache is refreshed automatically if the fit no longer matches it.

	\code{\link{confint}} does exactly the same thing under the familiar name from \code{lm} and \code{glm}.
}
\value{
	The fit, with its \code{conf.set} slot updated.  Levels build up, so asking for a new one keeps the old ones.  Read one back with \code{conf.set(fit)[["0.95"]]}, which gives a list with:
	\item{level}{The confidence level.}
	\item{cpts}{The best changepoints, ending with n.}
	\item{segmentations}{The plausible segmentations, grouped by how many changepoints they have (n counts as one), one per row.  The best one is included.}
	\item{gaps}{Same shape as \code{segmentations}: for each segmentation, how much more it costs than the best one.  0 for the best segmentation, larger for less plausible ones, never more than \code{budget}.  This is a cost difference on the same scale as the budget, NOT a log-likelihood (the cost functions drop a constant term), so read it only as "how far from optimal".}
	\item{budget}{The cost budget for this level.  Every gap lies in [0, budget]; a bigger level gives a bigger budget and so a bigger set.}
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

# store a confidence set, then look at it
fit = confidence_set(fit, level=0.95)
conf.set(fit)[["0.95"]]

# add a second level, no PELT re-run this time
fit = confint(fit, level=0.99)
conf.set(fit)                              # every level at once

print(conf.set(fit)[["0.99"]], n=Inf)      # show all the alternatives
}
\keyword{methods}
\keyword{univar}
\keyword{models}
\keyword{ts}
