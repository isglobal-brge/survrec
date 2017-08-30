\name{survdiffr}
\alias{survdiffr}

\title{ Test median survival differences (or other quantile)}
\description{
  Obtain bootstrap replicates of the median survival time for different groups of
  subjects. We can compute confidence intervals using boot package. 
}
\usage{
survdiffr(formula, data, q, B = 500, boot.F = "WC", boot.G = "none", ...)
}

\arguments{
  \item{formula}{ A formula object. If a formula object is supplied it must have
a Survr object as the response on the left of the ~ operator and a term on the
right. For a single bootstrap median survival the "~1" part of the formula is required. }
  \item{data}{A data frame in wich to interpret the variables named in the formula.}
  \item{q}{ Quantile that we are interested in to obtain a bootstrap sample 
            from survival function }
  \item{B}{ Number of boostrap samples }
  \item{boot.F}{ a character string specifying the boostrap procedure. 
  Possible value are either "PSH" or "WC" for nonparametric boostrap 
  or "semiparametric" for semiparametric boostrap. The default
  is "WC". Only the first words are required, e.g "P","W","se" }
  \item{boot.G}{a character string specifying if we also resample form censored empirical 
  distribution. Possible value are either "none" or "empirical". The default is "none".
  Only the first words are required, e.g "n","e"}
  \item{\dots}{additional arguments passed to the type of estimator.}
}
\details{
  See reference. Some procedures can be slow
}
\value{
  A boot object. Bootstrap confidence intervals can be computed using boot.ci function from
  boot package
}
\references{
Gonzalez JR, Peña EA. Bootstraping median survival with recurrent event data. 
IX Conferencia Española de Biometría; 2003 May 28-30; A Coruña, España.

Paper available upon request to the mantainer

}


\seealso{\code{\link{survfitr}},\code{\link{boot.ci}}}
 
 
\examples{

data(colon)

# We will compare the median survival time for three dukes stages
fit<-survdiffr(Survr(hc,time,event)~as.factor(dukes),
               data=colon, q=0.5)
boot.ci(fit$"1", type=c("norm","basic", "perc"))
boot.ci(fit$"2", type=c("norm","basic", "perc"))
boot.ci(fit$"3", type=c("norm","basic", "perc"))

# 75th quantile of survival function
fit<-survdiffr(Survr(hc,time,event)~as.factor(dukes),data=colon,q=0.75)
# bootstrap percentile confidence interval
quantile(fit$"1"$t,c(0.025,0.975))
quantile(fit$"2"$t,c(0.025,0.975))
quantile(fit$"3"$t,c(0.025,0.975))

# We could execute this if there is none Inf value
# boot.ci(fit$"1")
# boot.ci(fit$"2")
# boot.ci(fit$"3")


# We can modify the bootstrap procedure modifiying boot.F parameter
fit<-survdiffr(Survr(hc,time,event)~as.factor(dukes),data=colon,q=0.5,boot.F="PSH")
# bootstrap percentile confidence interval
quantile(fit$"1"$t,c(0.025,0.975))
quantile(fit$"2"$t,c(0.025,0.975))
quantile(fit$"3"$t,c(0.025,0.975))


}

\keyword{survival}
