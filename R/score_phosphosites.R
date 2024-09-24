#' Get a list of position specific weight matrices (PWMs) for human kinases
#'
#' The function returns a named list of 396 kinase PWMs. Among these are 303
#' serine/threonine kinases, 78 canonical tyrosine kinases and 15 non-canonical
#' tyrosine kinases i.e. dual-specific kinases, indicated by the '_TYR' suffix.
#'
#' Each PWM stores the log2-odds score per amino acid (23 rows) and position (10
#' or 11 columns) in matrix format. Beside the 20 standard amino acids also
#' phosphorylated serine, threonine and tyrosine residues are included.
#'
#' The central phospho-acceptor position of each PWM is at position `0` (column
#' 6). For serine/threonine specific kinases this position quantifies the
#' favorability of serine over threonine, but can be omitted when setting
#' 'includeSTfavorability=FALSE'.
#'
#' The specificity of a kinase PWM is controlled by parameter
#' 'matchAcceptorSpecificity'. It is set to 'TRUE' by default. Sites without a
#' matching acceptor are scored with -Inf in this case.
#'
#' @param includeSTfavorability Include serine vs. threonine favorability for
#'   the central phospho-acceptor?
#'
#' @param matchAcceptorSpecificity Only score sites with a matching central
#'   phospho-acceptor?
#'
#' @return A named list of numeric matrices (PWMs).
#'
#' @importFrom checkmate assert_logical
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @references Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al. An atlas of
#'   substrate specificities for the human serine/threonine kinome. Nature 613,
#'   759–766 (2023). https://doi.org/10.1038/s41586-022-05575-3
#'
#'   Yaron-Barir, T.M., Joughin, B.A., Huntsman, E.M. et al. The intrinsic
#'   substrate specificity of the human tyrosine kinome. Nature 629, 1174–1181
#'   (2024). https://doi.org/10.1038/s41586-024-07407-y
#'
#' @examples
#' pwms <- getKinasePWM()
getKinasePWM <- function(includeSTfavorability=TRUE,
                         matchAcceptorSpecificity=TRUE) {
  
    checkmate::assert_logical(includeSTfavorability)
    checkmate::assert_logical(matchAcceptorSpecificity)
    
    .convert_to_matrix <- function(x,
                                   includeSTfavorability=TRUE,
                                   acceptorSpecificity=NULL) {
        x <- x |> 
            dplyr::select(-Matrix) |> 
            tidyr::pivot_wider(values_from=Score, 
                               names_from=Position) 
        y <- as.matrix(x |> dplyr::select(-AA))
        rownames(y) <- x |> dplyr::pull(AA)
        if (!includeSTfavorability)
            y[c('S','T'),'0'] <- NA_real_
        y
        if (!is.null(acceptorSpecificity))
            y[setdiff(rownames(y), acceptorSpecificity),'0'] <- -Inf
        y
    }
    
    pwm <- read.csv( .JohnsonKinasePWM() )
    st <- lapply(split(pwm, pwm$Matrix), .convert_to_matrix, 
                 includeSTfavorability=includeSTfavorability,
                 acceptorSpecificity=if (matchAcceptorSpecificity) c('S','T','s','t')
                 )

    pwm <- read.csv( .TyrosineKinasePWM() )
    ty <- lapply(split(pwm, pwm$Matrix), .convert_to_matrix, 
                 includeSTfavorability=FALSE,
                 acceptorSpecificity=if (matchAcceptorSpecificity) c('Y','y')
    )
    
    c(st, ty)
}

#' Map log2-odds score to percentile rank
#'
#' For each kinase PWM, get a function that maps its log2-odds score to the
#' percentile rank in the background score distribution. The percentile rank of
#' a given score is the percentage of scores in corresponding background score
#' distribution that are less than or equal to that score. The background score
#' distribution per PWM is derived from matching each PWM to either the 85'603
#' unique phosphosites published in Johnson et al. 2023 (serine/threonine PWMs)
#' or the 6659 unique phosphosites published in Yaron-Barir et al. 2024
#' (tyrosine PWMs).
#'
#' Note: since the background sites don't contain non-central phosphorylated
#' residues (phospho-priming), the percentile rank of an input site which
#' includes phospho-priming will be capped to 100, if its PWM score exceeds the
#' largest observed background score for that PWM.
#'
#' Internally, [stats::approxfun] is used to linearly interpolate between the
#' PWM score and its 0.1% - quantile in the distribution over background scores.
#' This approximation allows for a lower memory footprint compared with the full
#' set of background scores.
#'
#' @return A named list of functions, one for each kinase PWM. Each function is
#'   taking a vector of log2-odds scores and maps them to a percentile rank in
#'   the range 0 to 100.
#'
#' @importFrom stats approxfun
#'
#' @export
#'
#' @examples
#' maps <- getScoreMaps()
getScoreMaps <- function() {
    PR <- cbind(
        read.csv( .JohnsonKinaseBackgroundQuantiles() ),
        read.csv( .TyrosineKinaseBackgroundQuantiles() )[,-1]
    )
    lapply(PR[,-1], function(score, quant) {
        stats::approxfun(score, quant, 
                         yleft=0, yright=100, 
                         ties=min)
    }, quant = 100 * PR[,"Quantiles"])
}


#' Low level single PWM scoring function. All unmatched characters (like "_" or
#' ".") evaluate to NA and don't contribute to the score sum. All strings must
#' be left aligned with the PWM. Characters extending the width of the PWM are 
#' excluded.
#' @noRd
.scoreSinglePWM <- function(sites, pwm) {
    vapply(sites, 
           FUN=function(aa) {
               aa <- substr(aa, 1, ncol(pwm))
               aa_score <- pwm[cbind(base::match(aa,rownames(pwm)), 
                                     seq_along(aa))]
               sum(aa_score, na.rm=TRUE) 
           }, 
           NA_real_)
}

#' Low level multiple PWM scoring function
#' @noRd
.scoreMultiplePWM <- function(sites, pwms, BPPARAM) {
    BiocParallel::bplapply(
        pwms,
        FUN=function(pwm, sites) {
            .scoreSinglePWM(sites, pwm)
        },
        sites=sites,
        BPPARAM=BPPARAM
    )
}

#' log2-odds to percentile rank
#' @noRd
.lodToPR <- function(scores, scoreMaps, BPPARAM) {
    BiocParallel::bpmapply(
        FUN=function(score, scoreMap) {
            scoreMap(score)
        },
        scores,
        scoreMaps[names(scores)],
        SIMPLIFY=FALSE,
        USE.NAMES=TRUE,
        BPPARAM=BPPARAM)
}

#' Match kinase PWMs to processed phosphosites
#'
#' `scorePhosphosites` takes a list of kinase PWMs and a vector of processed
#' phosphosites as input and returns a matrix of match scores per PWM and site.
#'
#' The match score is either the log2-odds score (`lod`) or the percentile rank
#' (`percentile`) in the background score distribution.
#'
#' @param pwms List with kinase PWMs as returned by [getKinasePWM].
#' @param sites A character vector with phosphosites. Check
#'   [processPhosphopeptides] for the correct phosphosite format.
#' @param scoreType Log2-odds score or percentile rank.
#' @param BPPARAM A [BiocParallelParam] object specifying how parallelization
#'   should be performed.
#'
#' @return A numeric matrix of size `length(sites)` times `length(pwms)`.
#'
#' @importFrom checkmate assert_list
#' @importFrom checkmate assert_character
#' @importFrom checkmate assert_class
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel SerialParam
#'
#' @export
#'
#' @seealso [getKinasePWM] for getting a list of kinase PWMs,
#'   [processPhosphopeptides] for the correct phosphosite format, and
#'   [getScoreMaps] for mapping PWM scores to percentile ranks
#'
#' @examples
#' score <- scorePhosphosites(getKinasePWM(), c("TGRRHTLAEV", "LISAVSPEIR"))
scorePhosphosites <- function(pwms, sites, 
                              scoreType=c('lod', 'percentile'), 
                              BPPARAM=BiocParallel::SerialParam()) {
    
    checkmate::assert_list(pwms, types=c("numeric", "matrix"), 
                           unique=TRUE, any.missing=FALSE, names='named')
    
    checkmate::assert_character(sites)
    
    scoreType <- match.arg(scoreType)
    
    checkmate::assert_class(BPPARAM, "BiocParallelParam")
    
    splitted <- strsplit(sites, split="")
    
    scores <- .scoreMultiplePWM(splitted, pwms, BPPARAM)
    
    if (scoreType == "percentile") {
        scoreMaps <- getScoreMaps()
        
        if (!all(names(scores) %in% names(scoreMaps))) {
            stop('Not all PWMs have score maps!')
        } 
        
        scores <- .lodToPR(scores, scoreMaps, BPPARAM)
    }
    
    scores <- do.call(cbind, scores)
    dimnames(scores) <- list(sites, names(pwms))
    
    scores
}
