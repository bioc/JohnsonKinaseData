.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ExperimentHub::createHubAccessors(pkgname, titles)
}

#' Get annnotation data for all 303 human serine/threonine kinase PWMs
#'
#' @return A data frame with columns MatrixName, GeneName, UniprotID, EntrezID,
#'   Description and KinaseFamily
#' @export
#'
#' @examples
#' anno <- getKinaseAnnotation()
getKinaseAnnotation <- function() {
  read.csv( JohnsonKinaseAnnotation() )
}

#' Get precomputed PWM scores for a large set of curated human phosphosites
#'
#' @return A data frame with log2-odds scores per phosphosite (rows) and PWMs
#'   (columns)
#' @export
#'
#' @examples
#' bg <- getBackgroundScores()
getBackgroundScores <- function() {
  read.csv( JohnsonKinaseBackgroundScores() )
}


