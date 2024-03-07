#' Process phosphosites to a common format used for PWM matching
#'
#' The central phospho-acceptor of each phosphosite is recognized in two ways.
#' Either the acceptor residue is followed by an asterisk (*) character, e.g.,
#' "SAGLLS*DEDC". Alternatively, it is defined as the central residue i.e. at
#' position `floor(width(site))+1`.
#'
#' Non-central phospho-acceptors (phospho-priming) should be indicated by the
#' lower case letters "s", "t" or "y". If phospho-priming is disabled these
#' residues are converted to upper case letters.
#'
#' The input sites are truncated and/or padded such that the processed sites are
#' of width 10 and have the central phospho-acceptor surrounded by 5 upstream
#' and 4 downstream residues.
#'
#' A warning is raised if the central phospho-acceptor is not serine or
#' threonine, as these sites are not covered by the Johnson PWMs.
#'
#' @param sites Character vector with phosphosites
#' @param allow_phospho_priming Allow phospho-acceptors at non-central
#'   positions? These should be indicated by the lower case letters "s", "t" or
#'   "y".
#'
#' @return A tibble with columns: `sites`, `processed`, `residue`
#' @export
#'
#' @examples
#' proc_sites <- process_phosphosites(c("SAGLLS*DEDC", "GDS*ND", "EKGDSN__", "___LySDEDC", "EKGtS*N"))
process_phosphosites <- function(sites,
                                 allow_phospho_priming = TRUE) {

  checkmate::assert_character(sites)
  checkmate::assert_logical(allow_phospho_priming)
  
  data <- tidyr::tibble(sites) |>
    dplyr::mutate(width = stringr::str_width(sites),
                  count = stringr::str_count(sites, "\\*")) 
  
  if ( any(data$count > 1) )
    stop('Data contains sites with multiple phosphorylations. Considered replacing non-central S*/T*/Y* by the lower case letters s/t/y.')
  
  data <- data |>
    dplyr::mutate(left = ifelse(count > 0,
                                stringr::str_replace(sites, stringr::regex(paste0("\\*", '(\\w*)$')), ''),
                                stringr::str_sub(sites, end = floor(width/2) + 1)),
                  right = ifelse(count > 0,
                                 stringr::str_replace(sites, stringr::regex(paste0('^(\\w*)', "\\*")), ''),
                                 stringr::str_sub(sites, start = floor(width/2) + 2))) |>
    dplyr::mutate(left = ifelse(stringr::str_width(left) <= 5, 
                                stringr::str_pad(left, pad = '_', width = (5+1), side = "left"),
                                stringr::str_trunc(left, width = (5+1), side = "left", ellipsis = "")),
                  right = ifelse(stringr::str_width(right) < 4,
                                 stringr::str_pad(right, pad = '_', width = 4, side = "right"),
                                 stringr::str_trunc(right, width = 4, side = "right", ellipsis = ""))) |>
    dplyr::mutate(residue = stringr::str_extract(left, '\\w$'))
  
  if (any(!data$residue %in% c('S','T')))
    warning('No S/T at central phospho-acceptor position.')
  
  data <- data |> 
    dplyr::mutate(processed = stringr::str_c(left, right)) |>
    dplyr::select(sites, processed, residue)
  
  if (!allow_phospho_priming) {
    data <- data |> 
      dplyr::mutate(processed = stringr::str_to_upper(processed)) 
  }
  
  data
}
