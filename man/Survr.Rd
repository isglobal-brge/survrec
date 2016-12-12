\name{Survr}
\alias{Survr}
\alias{is.Survr}
\title{ Create a Survival recurrent object}
\description{
 Create a survival recurrent object, usually used as a response variable in a model formula 
}
\usage{
Survr(id, time, event)
is.Survr(x) 
}
\arguments{
  \item{id}{Identificator of each subject. This value is the same for all recurrent times of each subject.}
  \item{time}{time of reccurence. For each subject the last time are censored.}
  \item{event}{The status indicator, 0=no recurrence 1=recurrence. Only this values are accepted.}
  \item{x}{any R object.}
}
\value{
An object of class Survr. Survr objects are implemented as a matrix of 3 colummns. No method for print.

In the case of is.Survr, a logical value T if x inherits from class "Survr", otherwise an F.

}

\seealso{  \code{\link{survfitr}},\code{\link{psh.fit}},\code{\link{wc.fit}},\code{\link{mlefrailty.fit}}  }

\examples{
 data(MMC)
 Survr(MMC$id,MMC$time,MMC$event)
}

\keyword{survival}
