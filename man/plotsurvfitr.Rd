\name{plot.survfitr}
\alias{plot.survfitr}
\alias{lines.survfitr}

\title{Plots estimated survival function from an object of class `survrec'.}

\description{Additional plots can be added to the same of axes using `lines.survrec'.}

\usage{
\method{plot}{survfitr}(x, conf.int=TRUE, prob = FALSE, ...)
}

\arguments{
 \item{x}{
    Object of class \code{survrec} (output from calling \code{survrec}
    function).}
 \item{conf.int}{Print the pointwise confidence intervals of the probability or 
                 survival function if its value is TRUE or FALSE.}
 \item{prob}{Print of the probability or survival function if its value is 
              TRUE or FALSE respectively.} 
 \item{...}{additional arguments passed to the plot function.}
}

\value{
  Print a plot of class \code{survrec}
}

\seealso{\code{\link{psh.fit}} \code{\link{mlefrailty.fit}}  \code{\link{wc.fit}}}

\keyword{survival}
