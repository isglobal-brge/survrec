\name{summary.survfitr}
\alias{summary.survfitr}
\alias{print.summary.survfitr}

\title{ Summary of a Survival of Recurrences Curve }
\description{
  Returns a matrix containing the survival curve and other information. If 
there are multiple curves, returns a list that contains the previous matrix for
each curve.  
}

\usage{
\method{summary}{survfitr}(object,...)
}

\arguments{
  \item{object}{output from a call to survfitr, psh.fit, wc.fit or mlefrailty.fit.}
  \item{...}{other unused arguments.}
}

\value{
  For one survival curve returns a matrix, and for multiple curves a list with the
same matrix for each curve. This matrix contains the distinct failure times, and the number
of events, at risk subjects, survival and standard error for each distinct time
}

\seealso{ \code{\link{survfitr}} }

\examples{
data(MMC)
summary(survfitr(Survr(id,time,event)~group,data=MMC))
}
\keyword{ survival}
