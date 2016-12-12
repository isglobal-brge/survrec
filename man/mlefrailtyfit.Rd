\name{mlefrailty.fit}
\alias{mlefrailty.fit}
\title{Survival function estimator for correlated recurrence time data 
under a Gamma Frailty Model}

\description{
  Estimation of survival function for correlated recurrence time data under a Gamma Frailty model using the maximum likelihood criterion. 
  The resulting object of class "survfitr" is plotted by `plot.survfitr', before it is returned.
}

\usage{
mlefrailty.fit(x,tvals, lambda=NULL, alpha=NULL, alpha.min, alpha.max, 
tol=1e-07, maxiter=500,alpha.console=TRUE)
}


\arguments{
  \item{x}{a survival recurrent event object.}
  \item{tvals}{vector of times where the survival function can be estimated. }
  \item{lambda}{optional vector  of baseline hazard probabilities at t (see details).
   Default is numdeaths/apply(AtRisk,2,sum).}
  \item{alpha}{optional parameter of shape and scale for the frailty distribution. If this 
               parameter is unknown is estimate via EM algorithm. In order to obtain the 
               convergence of this algorithm a seed is calculated (see details).}
  \item{alpha.min}{optional left bound of the alpha parameter in order to obtain a seed to 
                   estimate alpha parameter. Default value is 0.5.}
  \item{alpha.max}{optional rigth bound of the alpha parameter in order to obtain a seed to 
                    estimate alpha  parameter. Default value is the maximum of distinct times 
                    of events.}
  \item{tol}{optional tolerance of the EM algorithm used to estimate the alpha parameter. Default is 10e-7}
  \item{maxiter}{optional maximum number of iterations of the EM algorithm used to estimate 
                 the alpha parameter. Default is 500. }
  \item{alpha.console}{if TRUE prints in the console the estimates initial value for alpha and
                       the alpha estimate via the EM algorithm, if FALSE not.}
}

\details{
The product limit estimator developed by Peña, Strawderman and Hollander (2001)
are valid when the interoccurrence times are assumed to represent an IID sample from 
some underlying distribution F. This assumption is clearly restrictive in biomedical
applications, and one obvious generalization that allows association between interocurrence
times is a frailty model. 

A common and convenient choice of frailty distribution is a gamma distribution with shape and scale parameters
set equal to an unknown parameter \eqn{\alpha}. The common marginal survival function can be written
as following

\deqn{\bar {F}(t) =  {\left[ {{\frac{{\alpha }}{{\alpha + \Lambda
_{0} \left( {t} \right)}}}} \right]}^{\alpha}}{%
         1-F(t)=(alpha/(alpha+Lambda_0(t))^alpha}

The parameter \eqn{\alpha} controls the degree of association between interoccurrence 
times within a unit. Peña, Strawderman and Hollander (2001) showed that the estimation of 
\eqn{\alpha} and \eqn{\Lambda_0} can be obtained via the maximisation of the marginal likelihood function 
and the expectation-maximisation (EM) algorithm. For details and the theory behind this
estimator, please refer to Peña, Strawderman and Hollander (2001, JASA). 

In order to obtain a good convergence, \eqn{\alpha} is estimated previously. This 
estimation is used as a initial value in the EM procedure and it's carried out by the 
maximisation of the profile likelihood for \eqn{\alpha}. In this case the arguments 
of \bold{mlefrailty.fit} function called alpha.min and alpha.max are the boundaries 
of this maximisation. The maximum is obtained using the golden section search method.
}

\value{

If the convergence of EM algorithm is not obtained, the initial value of alpha can
be used as a alpha.min argument and recalculate.

\item{n}{number of unit or subjects observed.}
\item{m}{vector of number of recurrences in each subject (length n).}
\item{failed}{vector of number of recurrences in each subject (length n*m). 
              Vector ordered (e.g. times of first unit, times of second unit
              , ..., times of n-unit).}
\item{censored}{vector of times of censorship for each subject (length n).}
\item{numdistinct}{number of distinct failures times.}
\item{distinct}{vector of distinct failures times.}
\item{status}{0 if the estimation is can be provided and 1 if not depending 
              if alpha could be estimate or not. }
\item{alpha}{parameter of Gamma Frailty Model.}
\item{lambda}{Estimates of the hazard probabilities at distinct failures times.} 
\item{survfunc}{vector of survival estimated in distinct times.}
\item{tvals}{copy of argument.}
\item{MLEAttvals}{vector of survival estimated in tvals times.}
}

\references{ Peña, E.A., Strawderman, R. and Hollander, M. (2001). Nonparametric 
             Estimation with Recurrent Event Data. \emph{J. Amer. Statist. Assoc} 
             \bold{96}, 1299-1315.}

\seealso{ \code{\link{survfitr}} \code{\link{Survr}}}

\examples{

data(MMC)
fit<-mlefrailty.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
plot(fit)

# compare with pena-straderman-hollander

fit<-psh.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=2)

# and with wang-chang

fit<-wc.fit(Survr(MMC$id,MMC$time,MMC$event))
fit
lines(fit,lty=3)

}


\keyword{survival}
