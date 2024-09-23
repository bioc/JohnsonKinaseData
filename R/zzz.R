.onLoad <- function(libname, pkgname) {
    fl <- system.file("extdata", "metadata.csv", package=pkgname)
    titles <- read.csv(fl, stringsAsFactors=FALSE)$Title
    ## create hub accessors from titles as internal functions    
    ns <- asNamespace(pkgname)
    for (title in titles) {
        interal <- paste0('.',title)
        assign(interal, ExperimentHub:::.hubAccessorFactory(pkgname, title), 
               envir = ns)
        namespaceExport(ns, interal)
    }
}

#' Get annnotation data for all kinase PWMs
#'
#' The annotation data records for each kinase PWM, the PWM matrix name, gene
#' symbol and description, Uniprot ID, Entrez ID, acceptor specificity, kinase
#' sub-type, as well as the kinase family.
#'
#' Kinase PWMs are either serine/threonine or tyrosine specific in their central
#' phospho-acceptor. For dual-specific kinases, the tyrosine-specific PWM is
#' indicated by the '_TYR' suffix in the PWM name.
#'
#' Tyrosine kinases are further distinguished by sub-type and include receptor
#' tyrosine kinases (RTK), non-receptor tyrosine kinases (nRTK) and
#' non-canonical tyrosine kinases (ncTK) with dual-specificity.
#'
#' @return A data frame with columns MatrixName, GeneName, UniprotID, EntrezID,
#'   Description, AcceptorSpecificity, KinaseSubType and KinaseFamily
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom ExperimentHub createHubAccessors
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' anno <- getKinaseAnnotation()
getKinaseAnnotation <- function() {
    st <- read.csv( .JohnsonKinaseAnnotation() )
    st <- dplyr::mutate(AcceptorSpecificity = 'Ser/Thr', 
                        KinaseSubType = NA,
                        .after='EntrezID')
    
    ty <- read.csv( .TyrosineKinaseAnnotation() )
    ty <- dplyr::mutate(AcceptorSpecificity = 'Tyr', 
                        .after='EntrezID') |>

    dplyr::bind_rows(st, ty)
 }

#' Get precomputed PWM scores for two sets of curated human phosphosites
#'
#' Two sets of background scores are provided:
#'
#' 1. Ser/Thr phosphosites published in Johnson et al. 2023 2. Tyr phosphosites
#' published in Yaron-Barir et al. 2024
#'
#' The background scores are derived from matching the corresponding PWMs to
#' each of the sets. The resulting data frames contain the log2-odds score per
#' phosphosite and PWM.
#'
#' @param phosphoAcceptor Return background scores for either Ser/Thr or Tyr
#'   phosphosites
#'
#' @return A data frame with log2-odds scores per phosphosite (rows) and PWMs
#'   (columns)
#'
#' @references Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al. An atlas of
#'   substrate specificities for the human serine/threonine kinome. Nature 613,
#'   759–766 (2023). https://doi.org/10.1038/s41586-022-05575-3
#'
#'   Yaron-Barir, T.M., Joughin, B.A., Huntsman, E.M. et al. The intrinsic
#'   substrate specificity of the human tyrosine kinome. Nature 629, 1174–1181
#'   (2024). https://doi.org/10.1038/s41586-024-07407-y
#'
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom ExperimentHub createHubAccessors
#' @importFrom checkmate assert_choice
#' @importFrom utils read.csv
#'
#' @export
#'
#' @examples
#' bg <- getBackgroundScores()
getBackgroundScores <- function(phosphoAcceptor = c('Ser/Thr','Tyr')) {
    
    checkmate::assert_choice(phosphoAcceptor = c('Ser/Thr','Tyr'))
    
    if (phosphoAcceptor = 'Ser/Thr') {
        read.csv( .JohnsonKinaseBackgroundScores() )        
    } else {
        read.csv( .TyrosineKinaseBackgroundScores() )   
    }
}


