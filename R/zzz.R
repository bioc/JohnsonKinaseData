.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ExperimentHub::createHubAccessors(pkgname, titles)
}

#' get_kinase_annotation
#'
#' Annotation data for all 303 human kinase PWMs
#'
#' @return A data frame with columns MatrixName, GeneName, UniprotID, EntrezID, Description and KinaseFamily
#' @export
#'
#' @examples
#' anno <- get_kinase_annotation()
get_kinase_annotation <- function() {
  read.csv( JohnsonKinaseAnnotation() )
}

#' get_background_scores
#'
#' A data frame with precomputed PWM scores for a large set of curated
#' phosphosites
#'
#' @return Data frame with match scores per phosphosite (rows) and PWMs (columns)
#' @export
#'
#' @examples
#' bg <- get_background_scores()
get_background_scores <- function() {
  read.csv( JohnsonKinaseBackgroundScores() )
}


