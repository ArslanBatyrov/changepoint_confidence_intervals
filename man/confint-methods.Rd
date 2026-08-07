\name{confint-methods}
\docType{methods}
\alias{confint-methods}
\alias{confint,cpt-method}
\title{ ~~ Methods for Function confint  ~~}
\description{
 ~~ Methods for function \code{confint}  ~~
}
\section{Methods}{
\describe{

\item{\code{signature(object = "cpt")}}{
	Identical to \code{\link{confidence_set}} by design, provided so the feature is reachable through the name users know from lm/glm. Computes the confidence set of segmentations for a cpt object fitted with PELT and stores it in the conf.set slot, keyed by the confidence level. Returns the updated object, not an interval matrix.
}
}}
\keyword{methods}
\keyword{cpt}
\keyword{internal}
