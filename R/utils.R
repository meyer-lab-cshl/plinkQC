#' Check PLINK software access
#'
#' checkPlink checks that the PLINK software
#' (\url{https://www.cog-genomics.org/plink/1.9/}) required for perMarkerQC and
#' check_sex can be found from system call
#'
#' @param path2plink [character, optional] Absolute path to plink software. If
#' not provided, assumed that PATH set-up works and plink will be found by
#' system("plink").

checkPlink <- function(path2plink) {
    if (is.null(path2plink)) {
        findPlink <- tryCatch(system("plink -h", ignore.stdout=TRUE,
                                     ignore.stderr=TRUE),
                              warning=function(w) w)
        if("simpleWarning" %in% is(findPlink)) {
            stop("PLINK software required for running this function cannot be ",
                 "found in current PATH setting. Try to set path2plink.")
        }
    } else {
        findPlink <- tryCatch(system(paste(path2plink, "/plink -h", sep=""),
                                     ignore.stdout=TRUE, ignore.stderr=TRUE),
                              warning=function(w) w)
        if("simpleWarning" %in% is(findPlink)) {
            stop("PLINK software required for running this function cannot be ",
                 "found in path2plink.")
        }
    }
}
