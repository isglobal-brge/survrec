#' Migratory Motor Complex
#'
#' Times of the Migratory Motor Complex (MMC) of 19 healthy individuals,
#' measured through small bowel manometry.
#'
#' @format A data frame with 99 rows and 4 columns:
#' \describe{
#'   \item{id}{subject identifier, repeated for each recurrence}
#'   \item{time}{recurrence or censoring gap time}
#'   \item{event}{censoring status: 1 for every recurrence, 0 for the last
#'     (censored) time of each subject}
#'   \item{group}{a factor with levels `Males` and `Females`. Note: the
#'     groups were created at random to illustrate a group comparison}
#' }
#'
#' @source Husebye, E., Skar, V., Aalen, O.O. and Osnes, M. (1990).
#'   Digital ambulatory manometry of the small intestine in healthy adults.
#'   *Digestive Diseases and Sciences* **35**, 1057--1065.
#' @keywords datasets
"MMC"

#' Rehospitalization in colorectal cancer
#'
#' Rehospitalization times after surgery in patients diagnosed with
#' colorectal cancer.
#'
#' @format A data frame with 861 rows and 6 columns:
#' \describe{
#'   \item{hc}{subject identifier, repeated for each recurrence}
#'   \item{time}{rehospitalization or censoring gap time}
#'   \item{event}{censoring status: 1 for every rehospitalization, 0 for
#'     the last (censored) time of each subject}
#'   \item{chemoter}{did the patient receive chemotherapy? 1: no, 2: yes}
#'   \item{dukes}{Dukes' tumoral stage: 1: A-B, 2: C, 3: D}
#'   \item{distance}{distance from residence to hospital: 1: <= 30 km,
#'     2: > 30 km}
#' }
#'
#' @source González, J.R., Fernandez, E., Moreno, V. et al. (2005). Sex
#'   differences in hospital readmission among colorectal cancer patients.
#'   *Journal of Epidemiology and Community Health* **59**, 506--511.
#' @keywords datasets
"colon"
