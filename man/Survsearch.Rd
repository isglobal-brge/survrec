
\name{surv.search}
\alias{surv.search}

\title{Calculate the survival in selected times}

\description{

Auxiliary function called from pshPLE, wcPLE
and MLEFrailty.


The estimation using PLE (e.g. Kaplan-Meier) is a decreasing
constant piecewise function with jumps in the times with events. Thus,
to estimate the survival at any time we take the time of the precious event.


}

\usage{
  surv.search(tvals,time,surv) 
 }

\arguments{

\item{tvals}{vector of times where the survival function has to be estimated}
\item{time}{vector of failures times (distinct)}
\item{surv}{vector of survival of each time}
}

\value{
 Returns the survival in each selected time (tvals)
 from a vector of survival values
 }



\examples{
# we have the times 4,7,9,15,21,67
time<-c(4,7,9,15,21,67)

# and its survival (note: in this example there may be more
#                         than one event in some times)
surv<-c(0.8,0.7,0.65,0.55,0.43,0.22)

# We want to calculated the survival at times 1, 10, 32,64
surv.search(c(1,10,32,74),time,surv)

}


\keyword{survival}
