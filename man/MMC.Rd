\name{MMC}
\alias{MMC}
\non_function{}
\title{ Migratory Motor Complex}
\usage{data(MMC)}
\description{
 This contains the Migratoty Motor Complex data
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{ID of each subject. Repeated for each recurrence}
    \item{time}{recurrence o censoring time}
    \item{event}{censoring status. All event are 1 for each subject 
                 excepting last one that it is 0}
    \item{group}{a factor with levels}
  }
 Note: The group have been created (at random) to illustrate a group comparison

}
\source{
 Husebye E, Skar V, Aalen O and Osnes M (1990), Digestive Diseases and Sciences, p1057
}

\keyword{datasets}
