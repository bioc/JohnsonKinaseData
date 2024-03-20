#' Get a list of position specific weight matrices (PWMs)
#'
#' The function returns a list of PWMs for the 303 human serine/threonine
#' kinases originally published in Johnson et al. 2023. Each PWM stores the
#' log2-odds score per amino acid (23 rows) and position (10 columns) in matrix
#' format. Beside the 20 standard amino acids also phosphorylated serine,
#' threonine and tyrosine residues are included.
#'
#' The central phospho-acceptor position of each PWM is either serine or
#' threonine. By default, this position quantifies the favorability of serine
#' over threonine. This favorability can be omitted when setting
#' 'includeSTfavorability=FALSE' in which case the central position doesn't
#' contribute to the PWM score.
#'
#' @param includeSTfavorability Include serine vs. threonine favorability for
#'   the central phospho-acceptor?
#'
#' @return A named list of numeric matrices (PWMs).
#'
#' @importFrom checkmate assert_logical
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom tidyr pivot_wider
#'
#' @export
#'
#' @references Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al. An atlas of
#'   substrate specificities for the human serine/threonine kinome. Nature 613,
#'   759–766 (2023). https://doi.org/10.1038/s41586-022-05575-3
#'
#' @examples
#' pwms <- getKinasePWM()
getKinasePWM <- function(includeSTfavorability=TRUE) {
  
  checkmate::assert_logical(includeSTfavorability)
  
  pwms <- read.csv( JohnsonKinasePWM() )
  
  ## convert to a list of numeric matrices
  lapply(split(pwms, pwms$Matrix), function(x) {
    x <- x |> 
      dplyr::select(-Matrix) |> 
      tidyr::pivot_wider(values_from=Score, 
                         names_from=Position) 
    y <- as.matrix(x |> dplyr::select(-AA))
    rownames(y) <- x |> dplyr::pull(AA)
    if (!includeSTfavorability)
      y[,'0'] <- NA_real_
    y
  })
}

#' Map log2-odds score to percentile rank
#'
#' For each kinase PWM, get a function that maps its log2-odds score to the
#' percentile rank in the background score distribution. The percentile rank of
#' a given score is the percentage of scores in corresponding background score
#' distribution that are less than or equal to that score. The background score
#' distribution per PWM is derived from matching each PWM to the 85'603 unique
#' phosphosites published in Johnson et al. 2023.
#'
#' The background sites used by Johnson et al. don't contain non-central
#' phosphorylated residues (phospho-priming). Therefore any input sites which
#' include phospho-priming will be capped to 100 percentile rank, if their PWM
#' score exceeds the largest observed background score for that PWM.
#'
#' Internally, [stats::approxfun] is used to linearly interpolate between the
#' PWM score and its 0.1% - quantile in the distribution over background scores.
#' This approximation allows for a lower memory footprint compared with the full
#' set of background scores.
#'
#' @return A named list of functions, one for each kinase PWM. Each function is
#'   taking a vector of PWM log2-odds scores and maps them to a percentile rank
#'   in the range 0 to 100.
#'
#' @importFrom stats approxfun
#'
#' @export
#'
#' @examples
#' maps <- getScoreMaps()
getScoreMaps <- function() {
  PR <- read.csv( JohnsonKinaseBackgroundQuantiles() )
  lapply(PR[,-1], function(score, quant) {
    stats::approxfun(score, quant, 
                     yleft=0, yright=100, 
                     ties=min)
  }, quant = 100 * PR[,"Quantiles"])
}


#' Low level single PWM scoring function. All unmatched characters (like "_" or
#' ".") evaluate to NA and don't contribute to the score sum
#' @noRd
.scoreSinglePWM <- function(sites, pwm) {
  vapply(sites, 
         FUN=function(aa) {
           aa_score <- pwm[cbind(base::match(aa,rownames(pwm)), seq_along(aa))]
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
  
  if (!all(vapply(pwms, ncol, NA_integer_) == 10)) {
    stop('PWMs have incorrect dimension!')
  }
  
  checkmate::assert_character(sites)
  
  splitted <- strsplit(sites, split="")
  
  if (!all(vapply(splitted, length, NA_integer_) == 10)) {
    stop("'sites' elements are expected to be of length 10. Consider using 'processPhosphopeptides()'.")
  }
  
  scoreType <- match.arg(scoreType)
  
  checkmate::assert_class(BPPARAM, "BiocParallelParam")
  
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
