#' get_kinase_pwms
#'
#' Get a list of position weight matrices (PWMs) for the 303 human
#' serine/threonine kinases originally published in Johnson et al. 2023.
#'
#' @param include_ST_favorability Include serine vs. threonine favorability for
#'   the central phospho-acceptor?
#'
#' @return A named list of numeric matrices (PWMs).
#' @export
#'
#' @examples
#' pwms <- get_kinase_pwms()
get_kinase_pwms <- function(include_ST_favorability = TRUE) {
  
  checkmate::assert_logical(include_ST_favorability)
  
  pwms <- read.csv( JohnsonKinasePWM() )
  
  ## convert to a list of numeric matrices
  lapply(split(pwms, pwms$Matrix), function(x) {
    x <- x |> 
      dplyr::select(-Matrix) |> 
      tidyr::pivot_wider(values_from = Score, 
                         names_from = Position) 
    y <- as.matrix(x |> dplyr::select(-AA))
    rownames(y) <- x |> dplyr::pull(AA)
    if (!include_ST_favorability)
      y[,'0'] <- NA
    y
  })
}

#' get_score_maps
#'
#' For each kinase PWM, get a function that maps its log2-odds score to the
#' percentile rank in the background score distribution. The percentile rank of
#' a given score is the percentage of scores in corresponding background score
#' distribution that are less than or equal to that score. The background score
#' distribution per PWM is derived from matching each PWM to the 85'603 unique
#' phosphosites published in Johnson et al. 2023.
#'
#' Note that the background sites used by Johnson et al. do not contain any
#' non-central phosphorylated residues (phospho-priming). Therefore any input
#' sites which include phospho-priming will be capped to 100 percentile rank, if
#' their PWM score exceeds the largest observed background score for that PWM.
#'
#' Internally, stats::approxfun is used to linearly interpolate between the PWM
#' score and its 0.1% - quantile in the distribution over background scores.
#' This approximation allows for a lower memory footprint compared with the full
#' set of background scores.

#' @return A named list of functions, one for each kinase PWM. Each function is
#'   taking a vector of PWM log2-odds scores and maps them to a percentile rank
#'   in the range 0 to 100.
#' @export
#'
#' @examples
#' maps <- get_score_maps()
get_score_maps <- function() {
  PR <- read.csv( JohnsonKinaseBackgroundQuantiles() )
  lapply(PR[,-1], function(score, quant) {
    stats::approxfun(score, quant, 
                     yleft = 0, yright = 100, 
                     ties = min)
  }, quant = 100 * PR[,"Quantiles"])
}


#' Low level site scoring function
#' @export
.score_phosphosites <- function(sites, pwm) {
  sapply(sites, function(aa) {
    aa_score <- pwm[cbind(base::match(aa,rownames(pwm)), seq_along(aa))]
    # all unmatched characters (like "_" or ".") evaluate to NA and don't
    # contribute to the score sum
    sum(aa_score, na.rm = TRUE) 
  })
}

#' Match kinase PWMs to processed phosphosites
#'
#' `score_phosphosites` takes a list of kinase PWMs and a vector of processed
#' phosphosites as input and returns a matrix of match scores per PWM and site.
#'
#' The score is either the PWM match score (`log2_odds`) or the percentile rank
#' (`percentile`) in the background score distribution.
#'
#' @param pwms List with kinase PWMs as returned by [get_kinase_pwms()].
#' @param sites A character vector with phosphosites. Check
#'   [process_phosphosites()] for the correct phosphosite format.
#' @param score_type Percentile rank or log2-odds score.
#' @param BPPARAM A BiocParallelParam object specifying how parallelization
#'   should be performed.
#'
#' @return A numeric matrix of size `length(sites)` times `length(kinases)`.
#' @export
#'
#' @seealso [get_kinase_pwms()] for getting a list of kinase PWMs,
#'   [process_phosphosites()] for the correct phosphosite format, and
#'   [get_score_maps()] for mapping PWM scores to percentile ranks
#'
#' @examples
#' score <- score_phosphosites(get_kinase_pwms(), c("TGRRHTLAEV", "LISAVSPEIR"))
score_phosphosites <- function(pwms, sites, score_type = c('percentile', 'log2_odds'), 
                               BPPARAM = BiocParallel::MulticoreParam(1)) {

  checkmate::assert_list(pwms, types = c("numeric", "matrix"), 
                         unique = TRUE, any.missing = FALSE, names = 'named')
  
  if (!all(sapply(pwms, ncol) == 10)) {
    stop('PWMs have incorrect dimension!')
  }
  
  checkmate::assert_character(sites)
  
  splitted <- strsplit(sites, split="")
  
  if (!all(sapply(splitted, length) == 10)) {
    stop("'sites' elements are expected to be of length 10. Consider using 'process_phosphosites()'.")
  }
  
  score_type <- match.arg(score_type)
  
  checkmate::assert_class(BPPARAM, "BiocParallelParam")
  
  scores <- BiocParallel::bplapply(
    pwms, 
    FUN = function(pwm, sites) {
      .score_phosphosites(sites, pwm)
    },
    sites = splitted,
    BPPARAM = BPPARAM)
  
  if (score_type == "percentile") {
    score_maps <- get_score_maps()
    
    if (!all(names(scores) %in% names(score_maps))) {
      stop('Not all PWMs have score maps!')
    } 

    scores <- BiocParallel::bpmapply(
      FUN = function(score, score_map) {
        score_map(score)
      },
      scores,
      score_maps[names(scores)],
      BPPARAM = BPPARAM)
  } else {
    scores <- do.call(cbind, scores)    
  }
  
  dimnames(scores) <- list(sites, names(pwms))
  scores
}
