\name{print.survfitr}
\alias{print.survfitr}
\title{ Print a Short Summary of a Survival Recurrent Curve}
\description{
  Print number of observations, number of events, the restricted mean survival
and its standard error, the median survial and the minimum, maximum and median 
number of recurrences for each subject.
}
\usage{
\method{print}{survfitr}(x, scale=1, digits=max(options()$digits - 4, 3), ...)
}
\arguments{
  \item{x}{the result of a call to the survfit, psh.fit, wc.fil or mlefrailty.fit functions}
  \item{scale}{a numeric value to rescale the survival time, e.g., if the input data to
survfit were in days, scale=365 would scale the printout to years}
  \item{digits}{number of digits to print}
  \item{\dots}{ other unusued arguments }
}
\details{
   The restricted mean and its standard error are based on a truncated estimator. If the last
observation(s) is not a death, then the survival curve estimate does not fot to zero and the 
mean survival time cannot be estimated. Instead, the quantity reported is the mean of survival
restricted to the time before the last censoring. When the las censoring time is not random 
this quantity is occasionally of interest.

The median are defined by drawing a horizontal line at 0.5 on the plot of the survival curve.
  
}
\value{
x, with the invisible flag set.
 
The number of observations, the number of events, the restricted mean survival
and its standard error, the median survival and the minimum, maximum and median number of
recurrences are printed. If there are multiple curves, there is one line of output for each.
}


\seealso{ \code{\link{summary.survfitr}},  \code{\link{survfitr}} }

\examples{
data(MMC)
fit<-survfitr(Survr(id,time,event)~group,data=MMC)
print(fit)
}
\keyword{survival}
