\name{conf.set<-}
\alias{conf.set<-}
\title{
Generic Function - conf.set<-
}
\description{
Generic function
}
\usage{
conf.set(object)<-value
}
\arguments{
  \item{object}{
	Depending on the class of \code{object} depends on the method used (and if one exists)
}
  \item{value}{
	Replacement value
}
}
\details{
	Generic Function
}
\value{
	Depends on the class of \code{object}, see individual methods
}
\author{
Rebecca Killick
}


\seealso{
\code{\link{conf.set<--methods}}
}
\examples{
x=new("cpt") # new cpt object
conf.set(x)<-list("0.95"=c(46,54)) # replaces the current conf.set list in x
}

\keyword{methods}
\keyword{cpt}
\keyword{internal}
