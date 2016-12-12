\name{wc.fit}

\alias{wc.fit}

\title{Survival function estimator for recurrence time data using the 
estimator developed by Wang and Chang.}

\description{Estimation of survival function for correlated or i.i.d. recurrence 
time data by means of the product limit estimator (PLE) method developed by Wang and Chang.
The resulting object of class "survfitr" is plotted by `plot.survfitr', before it is returned.
}

\usage{
wc.fit(x,tvals) 
}

\arguments{
\item{x}{a survival recurrent event object.}
\item{tvals}{vector of times where the survival function can be estimated. }
}

\value{
\item{n}{number of unit or subjects observed.}
\item{m}{vector of number of recurrences in each subject (length n).}
\item{failed}{vector of number of recurrences in each subject (length n*m). Vector ordered (e.g. times of first unit, times of second unit, ..., times of n-unit)}
\item{censored}{vector of times of censorship for each subject (length n).}
\item{numdistinct}{number of distinct failures times.}
\item{distinct}{vector of distinct failures times.}
\item{AtRisk}{matrix of number of persons-at-risk at each distinct time and for each subject.}
\item{survfunc}{vector of survival estimated in distinct times.}
\item{tvals}{copy of argument.}
\item{PSHpleAttvals}{vector of survival estimated in tvals times.}
}

\details{
Wang and Chang (1999) proposed an estimator of the common marginal survivor function in the 
case where within-unit interrocurence times are correlated. The correlation structure 
considered by Wang and Chang (1999) is quite general and contains, in particular, 
both the i.i.d. and multiplicative (hence gamma) frailty model as special cases. 

This estimator removes the bias noted for the produc-limit estimator developed by 
Pea, Strawderman and Hollander (PSH, 2001) when interrocurence times are correlated within 
units. However, when applied to i.i.d. interocurrence times, this estimator is not expected 
to perform as well as the PSH estimator, especially with regard to efficiency.
}

\references{

Wang, M.-C. and Chang, S.-H. (1999). Nonparametric Estimation of a Recurrent Survival Function. \emph{J. Amer. Statist. Assoc} \bold{94}, 146-153.

}

\note{ The mantainer wishes to thank Professors Chiung-Yu Huang and Shu-Hui Chang
for their help for providing us with the Fortran code which computes
standard errors of Wang and Chang's estimator.

}


\seealso{ \code{\link{survfitr}} \code{\link{Survr}}}

\examples{

data(MMC)

fit<-wc.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
plot(fit,conf.int=FALSE)

# compare with pena-straderman-hollander

fit<-psh.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=2)

# and with MLE frailty

fit<-mlefrailty.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=3)

}



\keyword{survival}

