.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ExperimentHub::createHubAccessors(pkgname, titles)
}

#' Get annnotation data for all 303 human serine/threonine kinase PWMs
#'
#' The annotation data records for each of the 303 human serine/threonine kinase
#' PWMs originally published in Johnson et al. the PWM matrix name, gene symbol
#' and description, Uniprot ID and Entrez ID as well as the kinase family.
#'
#' @return A data frame with columns MatrixName, GeneName, UniprotID, EntrezID,
#'   Description and KinaseFamily
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom ExperimentHub createHubAccessors
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' anno <- getKinaseAnnotation()
getKinaseAnnotation <- function() {
  read.csv( JohnsonKinaseAnnotation() )
}

#' Get precomputed PWM scores for a large set of curated human phosphosites
#'
#' The background scores are derived from matching each PWM to the 85'603 unique
#' phosphosites published in Johnson et al. 2023. The data frame contains the
#' log2-odds score per phosphosite and PWM.
#'
#' @return A data frame with log2-odds scores per phosphosite (rows) and PWMs
#'   (columns)
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom ExperimentHub createHubAccessors
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' bg <- getBackgroundScores()
getBackgroundScores <- function() {
  read.csv( JohnsonKinaseBackgroundScores() )
}


