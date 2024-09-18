#' processPhosphopetides
#'
#' Process phospho-peptides to a common format used for PWM matching
#'
#' Phosphorylated residues are recognized either by lower case letters (s, t or
#' y) or the phosphorylated residue is followed by an asterisk (S*, T* or Y*).
#'
#' If a peptide reports several phosphorylated residues, parameter
#' `onlyCentralAcceptor` allows for two processing options: (1) By default, only
#' the central phospho-acceptor of each phospho-peptide is considered. Here
#' central is defined as the left-closest position to `floor(nchar(site)/2)+1`.
#' (2) All phospho-acceptors are considered as central in which case the
#' phospho-peptide is replicated and aligned for each of its phosphorylated
#' residues. In this case the output sites are not in parallel to the input
#' peptides.
#'
#' In both cases, non-central phospho-acceptors are indicated by lower case
#' letters (s, t, or y). These residues enable phospho-priming of the site. If
#' phospho-priming is disabled (parameter `allowPhosphoPriming`) these residues
#' are converted to upper case letters.
#'
#' If a site does not follow the phosphorylation patterns described above, the 
#' central residue defined by position `floor(nchar(site)/2)+1` is considered 
#' the default phospho-acceptor site.
#'
#' The input sites are truncated and/or padded such that the processed sites are
#' of width `upstream`+`downstream`+1. By default the central phospho-acceptor 
#' is surrounded by 5 upstream and 5 downstream residues.
#'
#' A warning is raised if the central phospho-acceptor is not serine, threonine 
#' or tyrosine.
#'
#' @param sites Character vector with phospho-peptides
#' @param onlyCentralAcceptor Process only the central phospho-acceptor residue?
#' @param allowPhosphoPriming Allow phospho-acceptors at non-central positions?
#'   These should be indicated by the lower case letters s, t or y.
#' @param upstream Number of bases upstream of central phospho-acceptor in the 
#'   processed output
#' @param downstream Number of bases downstream of central phospho-acceptor in 
#'   the processed output
#'
#' @return A tibble with columns: `sites`, `processed`, `acceptor`
#'
#' @importFrom checkmate assert_character
#' @importFrom checkmate assert_logical
#' @importFrom tidyr tibble
#' @importFrom tidyr unnest
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_min
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr full_join
#' @importFrom dplyr join_by
#' @importFrom stringr str_count
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_sub
#' @importFrom stringr regex
#' @importFrom stringr str_pad
#' @importFrom stringr str_trunc
#' @importFrom stringr str_extract
#' @importFrom stringr str_c
#' @importFrom stringr str_to_upper
#' @importFrom purrr map
#'
#' @export
#'
#' @examples
#' procSites <- processPhosphopeptides(c("AGLLS*DEDC", "RtEKGS*N", "ETGKDN"))
processPhosphopeptides <- function(sites,
                                   onlyCentralAcceptor=TRUE,
                                   allowPhosphoPriming=TRUE,
                                   upstream=5L,
                                   downstream=5L) {
    
    checkmate::assert_character(sites)
    checkmate::assert_logical(allowPhosphoPriming)
    checkmate::assert_numeric(upstream, lower=0, len=1)
    checkmate::assert_numeric(downstream, lower=0, len=1)
    
    ## unify acceptor pattern
    data <- tidyr::tibble(sites) |>
        dplyr::mutate(
            modified=stringr::str_replace_all(sites, 
                                              c('[S,s]\\*'='s', '[T,t]\\*'='t', 
                                                '[Y,y]\\*'='y')),
            center_pos=floor(nchar(modified)/2) + 1,
            is_lower=stringr::str_count(modified, "[s,t,y]") > 0)
    
    ## find possible acceptor locations
    locs <- data |>
        dplyr::mutate(hits=stringr::str_locate_all(modified, "[s,t,y]"),  
                      acceptor_pos=purrr::map(hits, function(df) df[,1])) |> 
        tidyr::unnest(acceptor_pos, keep_empty=TRUE) |>
        dplyr::mutate(acceptor_pos=dplyr::if_else(is_lower, 
                                                  acceptor_pos, 
                                                  center_pos)) 
    
    if (onlyCentralAcceptor) {
        ## select left-closest acceptor to central position 
        locs <- locs |>
            dplyr::group_by(modified) |>
            dplyr::mutate(diff=abs(acceptor_pos - center_pos)) |>
            dplyr::slice_min(diff, with_ties=TRUE) |>
            dplyr::arrange(acceptor_pos) |>
            dplyr::slice_head() |>
            dplyr::ungroup()
    }
    
    data <- data |>
        dplyr::full_join(locs |> dplyr::select(modified, acceptor_pos),
                         by = dplyr::join_by(modified))
    
    ## split sites in left and right parts
    data <- data |>
        dplyr::mutate(left=stringr::str_sub(modified, 
                                            end=acceptor_pos),
                      right=stringr::str_sub(modified, 
                                             start=acceptor_pos+1))
                      
    ## pad left/right parts to match upstream/downstream conditions
    data <- data |>
        dplyr::mutate(left=stringr::str_pad(left, pad='_', 
                                            width=upstream+1, side="left"),
                      left=stringr::str_trunc(left, width=upstream+1, 
                                              side="left", ellipsis=""),
                      right=stringr::str_pad(right, pad='_', 
                                             width=downstream, side="right"),
                      right=stringr::str_trunc(right, width=downstream, 
                                               side="right", ellipsis=""),
                      processed=stringr::str_c(left,right), 
                      acceptor=stringr::str_to_upper(
                          stringr::str_sub(processed, 
                                           start=upstream+1, 
                                           end=upstream+1)))
    
    stringr::str_sub(data$processed, start=upstream+1, end=upstream+1) <- 
        (data |> dplyr::pull(acceptor))
    
    if (any(!data$acceptor %in% c('S','T','K')))
        warning('No S/T/K at central phospho-acceptor position.')
    
    if (!allowPhosphoPriming) {
        data <- data |> 
            dplyr::mutate(processed=stringr::str_to_upper(processed)) 
    }
    
    data |> dplyr::select(sites, processed, acceptor)
}
