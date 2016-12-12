\name{psh.fit}

\alias{psh.fit}

\title{Survival function estimator for recurrence time data using the estimator 
developed by Peña, Strawderman and Hollander}

\description{
   Estimation of survival function for recurrence time data by means 
   the generalized product limit estimator (PLE) method developed by 
   Pe a, Strawderman and Hollander. The resulting object of class "survfitr" 
   is plotted by `plot.survfitr', before it is returned.
}

\usage{
   psh.fit(x,tvals) 
}

\arguments{
\item{x}{a survival recurrent event object}
\item{tvals}{vector of times where the survival function can be estimated. }
}

\details{ 
The estimator computed by this object is the nonparametric
estimator of the inter-event time survivor function under the
assumption of a renewal or IID model. This generalizes the
product-limit estimator to the situation where the event is
recurrent. For details and the theory behind this estimator,
please refer to Pe\~na, Strawderman and Hollander (2001, JASA).
}

\value{
\item{n}{number of unit or subjects observed.}
\item{m}{vector of number of recurrences in each subject (length n).}
\item{failed}{vector of number of recurrences in each subject (length n*m). Vector ordered (e.g. times of first unit, times of second unit, ..., times of n-unit)}
\item{censored}{vector of times of censorship for each subject (length n).}
\item{numdistinct}{number of distinct failures times}
\item{distinct}{vector of distinct failures times}
\item{AtRisk}{matrix of number of persons-at-risk at each distinct time and for each subject}
\item{survfunc}{vector of survival estimated in distinct times}
\item{tvals}{copy of argument}
\item{PSHpleAttvals}{vector of survival estimated in tvals times}
}

\references{Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric Estimation with Recurrent Event Data.  \emph{J. Amer. Statist. Assoc} \bold{96}, 1299-1315.}

\seealso{ \code{\link{survfitr}} \code{\link{Survr}}}

\examples{

data(MMC)
fit<-psh.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
plot(fit,conf.int=FALSE)

# compare with MLE Frailty

fit<-mlefrailty.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=2)

# and with wang-chang

fit<-wc.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=3)

}

\keyword{survival}
