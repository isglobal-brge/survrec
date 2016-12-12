\name{q.search}
\alias{q.search}

\title{ Calculate the survival time of a selected quantile }
\description{
  Auxiliary function called from survdiffr function. Given a survfitr object we obtain the 
  quantile from a survival function
}
\usage{
q.search(f, q = 0.5)
}

\arguments{
  \item{f}{ survdifr object }
  \item{q}{ quantile. Default is 0.5 }
}

\value{
 Returns the time in a selected quantile
}
\examples{

data(MMC)
fit<-survfitr(Survr(id,time,event)~1,data=MMC)

# 75th percentile from the survival function
q.search(fit,q=0.75)

}
\keyword{survival}
