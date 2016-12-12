\name{colon}
\alias{colon}
\non_function{}
\title{Rehospitalization colorectal cancer}
\usage{data(colon)}
\description{
 This contains rehospitalization times after surgery in patients with colorectal cancer
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{hc}{identificator of each subject. Repeated for each recurrence}
    \item{time}{rehospitalization o censoring time}
    \item{event}{censoring status. All event are 1 for each subject 
                 excepting last one that it is 0}
    \item{chemoter}{Did patient receive chemotherapy? 1: No 2:Yes}
    \item{dukes}{Dukes' tumoral stage: 1:A-B 2:C 3:D}
    \item{distance}{distance from living place to hospital 1:<=30 Km. 2:>30 Km.}
  }
 
}
\source{
 Gonz\'alez, JR., Fernandez, E., Moreno, V. et al. Gender differences in hospital 
 readmission among colorectal cancer patients. Currently submited to JCO
}

\keyword{datasets}
