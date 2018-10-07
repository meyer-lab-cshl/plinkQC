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

#' Test lists for different properties of numerics
#'
#' Test all elements of a list if they are numeric, positive numbers, integers
#' or proportions (range 0-1)
#'
#' @param numbers List whose elements are tested for being numeric
#' @param positives List whose elements are tested for being positive numbers
#' @param integers List whose elements are tested for being integers
#' @param proportions List whose elements are tested for being proportions
#' between 0 and 1
#' @return NULL
testNumerics <- function(numbers, positives=NULL, integers=NULL,
                         proportions=NULL) {
    notNullnotNumeric <- function(x) !is.null(x) && !is.numeric(x)
    notNullnotPositive <- function(x) !is.null(x)  && x <= 0
    notNullnotInteger <- function(x) !is.null(x)  && x%%1!=0
    notNullNotInRange <- function(x) !is.null(x)  && (x < 0 || x > 1)
    if (any(sapply(numbers, notNullnotNumeric))) {
        notNumbers <- which(sapply(numbers, notNullnotNumeric))
        stop(paste(names(numbers)[notNumbers], collapse=","),
             " is/are not numeric (",
             paste(numbers[notNumbers], collapse=","), ")", sep="")
    }
    if (!is.null(positives)) {
        if (any(sapply(positives, notNullnotPositive))) {
            notPositives <- which(sapply(positives, notNullnotPositive))
            stop(paste(names(positives)[notPositives], collapse=","),
                 " has/have to be greater than zero (",
                 paste(positives[notPositives], collapse=","), ")",
                 sep="")
        }
    }
    if (!is.null(integers)) {
        if (any(sapply(integers, notNullnotInteger))) {
            notIntegers <- which(sapply(integers, notNullnotInteger))
            stop(paste(names(integers)[notIntegers], collapse=","),
                 " has/have to be integers, given ",
                 paste(integers[notIntegers], collapse=","),
                 sep="")
        }
    }
    if (!is.null(proportions)) {
        if (any(sapply(proportions, notNullNotInRange))) {
            outOfRange <- which(sapply(proportions, notNullNotInRange))
            stop("Proportions have to be specified between 0 and 1: ",
                 paste(names(proportions)[outOfRange], collapse=","),
                 " are outside of this range (",
                 paste(proportions[outOfRange], collapse=","), ")", sep="")
        }
    }
}
