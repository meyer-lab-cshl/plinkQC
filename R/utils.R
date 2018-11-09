#' Check PLINK software access
#'
#' checkPlink checks that the PLINK software
#' (\url{https://www.cog-genomics.org/plink/1.9/})
#' can be found from system call.
#'
#' @param path2plink [character] Absolute path to PLINK executable
#' (\url{https://www.cog-genomics.org/plink/1.9/}) i.e.
#' plink should be accesible as path2plink -h. The full name of the executable
#' should be specified: for windows OS, this means path/plink.exe, for unix
#' platforms this is path/plink. If not provided, assumed that PATH set-up works
#' and PLINK will be found by \code{\link[sys]{exec_wait}}('plink').
#' @return Path to PLINK executable.
#' @export

checkPlink <- function(path2plink=NULL) {
    if (grepl("~", path2plink)) {
        stop("Path to plink (", path2plink,
             ") contains ~: please supply full path, not relying",
             " on tilde extension.")
    }
    if (is.null(path2plink)) {
        path2plink <- 'plink'
        preset <- FALSE
    } else {
        preset <- TRUE
    }
    if (.Platform$OS.type == 'windows') {
        path2plink <- paste(gsub('\\.exe', '', path2plink), '.exe', sep="")
    }
    findPlink <- tryCatch(sys::exec_wait(path2plink, args="-h", std_out=FALSE,
                                        std_err=TRUE),
                          warning=function(w) w,  error = function(e) e)
    if("simpleError" %in% is(findPlink)) {
        if (!preset) {
            stop("PLINK software required for running this function cannot be ",
                 "found in current PATH setting. Error message:", findPlink,
                 ". Try to set path2plink.")
        } else {
            stop("PLINK software required for running this function cannot be ",
                 "found in path2plink (", path2plink, ") . Error message: ",
                 findPlink, ".")
        }
    }
    return(path2plink)
}

#' Test lists for different properties of numerics
#'
#' Test all elements of a list if they are numeric, positive numbers, integers
#' or proportions (range 0-1).
#'
#' @param numbers [list] whose elements are tested for being numeric.
#' @param positives [list] whose elements are tested for being positive numbers.
#' @param integers [list] whose elements are tested for being integers.
#' @param proportions [list] whose elements are tested for being proportions.
#' between 0 and 1.
#' @return NULL
#' @export
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

#' Remove related individuals while keeping maximum number of individuals
#'
#' \code{filterRelatedness} takes a data.frame with pair-wise relatedness
#' measures of samples and returns pairs of individual IDs that are related as
#' well as a list of suggested individual IDs to remove.
#' \code{filterRelatedness} finds pairs of samples whose relatedness estimate is
#' larger than the specified relatednessTh. Subsequently, for pairs of
#' individual that do not have additional relatives in the dataset, the
#' individual with the worse otherCriterionMeasure (if provided) or arbitrarily
#' individual 1 of that pair is selected and returned as the individual failing
#' the relatedness check. For more complex family structures, the unrelated
#' individuals per family are selected (e.g. in a parents-offspring trio, the
#' offspring will be marked as fail, while the parents will be kept in the
#' analysis).
#' @param relatedness [data.frame] containing pair-wise relatedness estimates
#' (in column [relatednessRelatedness]) for individual 1 (in column
#' [relatednessIID1] and individual 2 (in column [relatednessIID1]). Columns
#' relatednessIID1, relatednessIID2 and relatednessRelatedness have to present,
#' while additional columns such as family IDs can be present. Default column
#' names correspond to column names in output of plink --genome
#' (\url{https://www.cog-genomics.org/plink/1.9/ibd}). All original
#' columns for pair-wise highIBDTh fails will be returned in fail_IBD.
#' @param relatednessTh [double] Threshold for filtering related individuals.
#' Individuals, whose pair-wise relatedness estimates are greater than this
#' threshold are considered related.
#' @param otherCriterion [data.frame] containing a QC measure (in column
#' [otherCriterionMeasure]) per individual (in column [otherCriterionIID]).
#' otherCriterionMeasure and otherCriterionIID have to present, while additional
#' columns such as family IDs can be present. IIDs in relatednessIID1 have to
#' be present in otherCriterionIID.
#' @param otherCriterionTh [double] Threshold for filtering individuals based on
#' otherCriterionMeasure. If related individuals fail this threshold they will
#' automatically be excluded.
#' @param otherCriterionThDirection [character] Used to determine the direction
#' for failing the otherCriterionTh. If 'gt', individuals whose
#' otherCriterionMeasure > otherCriterionTh will automatically be excluded.
#' For pairs of individuals that have no other related samples in the cohort:
#' if both otherCriterionMeasure < otherCriterionTh, the individual with the
#' larger otherCriterionMeasure will be excluded.
#' @param relatednessIID1 [character] Column name of column containing the IDs
#' of the first individual.
#' @param relatednessIID2 [character] Column name of column containing the IDs
#' of the second individual.
#' @param relatednessFID1 [character, optional] Column name of column containing
#' the family IDs of the first individual; if only relatednessFID1 but not
#' relatednessFID2 provided, or none provided even though present in
#' relatedness, FIDs will not be returned.
#' @param relatednessFID2 [character, optional] Column name of column containing
#' the family IDs of the second individual; if only relatednessFID2 but not
#' relatednessFID1 provided, or none provided even though present in
#' relatedness, FIDs will not be returned.
#' @param relatednessRelatedness [character] Column name of column containing
#' the relatedness estimate.
#' @param otherCriterionIID [character] Column name of column containing the
#' individual IDs.
#' @param otherCriterionMeasure [character] Column name of the column containing
#' the measure of the otherCriterion (for instance SNP missingness rate).
#' @return named [list] with i) relatednessFails, a [data.frame] containing the
#' reordered input data.frame relatedness after filtering for pairs of
#' individuals in relatednessIID1 and relatednessIID2, that fail the relatedness
#' QC and ii) failIDs, a [data.frame] with the [IID]s (and [FID]s if provided)
#' of the individuals that fail the relatednessTh.
#' @export

relatednessFilter <- function(relatedness, otherCriterion=NULL,
                              relatednessTh, otherCriterionTh=NULL,
                              otherCriterionThDirection=
                                  c("gt", "ge", "lt", "le", "eq"),
                              relatednessIID1="IID1", relatednessIID2="IID2",
                              relatednessFID1=NULL, relatednessFID2=NULL,
                              relatednessRelatedness="PI_HAT",
                              otherCriterionIID="IID",
                              otherCriterionMeasure=NULL ) {

    if (!(relatednessIID1 %in% names(relatedness))) {
        stop(paste("Column", relatednessIID1, "not found!"))
    }
    if (!(relatednessIID2 %in% names(relatedness))) {
        stop(paste("Column", relatednessIID1, "not found!"))
    }
    if (!(relatednessRelatedness %in% names(relatedness))) {
        stop(paste("Column", relatednessRelatedness, "not found!"))
    }
    relatedness_names <- names(relatedness)
    names(relatedness)[names(relatedness) == relatednessIID1] <- "IID1"
    names(relatedness)[names(relatedness) == relatednessIID2] <- "IID2"
    names(relatedness)[names(relatedness) == relatednessRelatedness] <- "M"

    if (!is.null(relatednessFID1) && is.null(relatednessFID2) ||
        is.null(relatednessFID1) && !is.null(relatednessFID2) ) {
        stop("Either none or both, relatednessFID1 and relatednessFID2 have to
             be provided")
    }
    if (!is.null(relatednessFID1) && !is.null(relatednessFID2)) {
        if (!(relatednessFID1 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID1, "not found!"))
        }
        if (!(relatednessFID2 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID2, "not found!"))
        }
        names(relatedness)[names(relatedness) == relatednessFID1] <- "FID1"
        names(relatedness)[names(relatedness) == relatednessFID2] <- "FID2"
    }
    relatedness_original <- relatedness
    relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~M)


    # individuals with at least one pair-wise comparison > relatednessTh
    highRelated <- dplyr::filter_(relatedness, ~M > relatednessTh)
    uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))

    fail_other <- NULL
    if (!is.null(otherCriterion)) {
        otherCriterionThDirection <- match.arg(otherCriterionThDirection)
        if (!(otherCriterionMeasure %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionMeasure, "not found!"))
        }
        if (!(otherCriterionIID %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionIID, "not found!"))
        }
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionMeasure] <- "M"
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionIID] <- "IID"

        fail_other <- apply(highRelated, 1, function(pair) {
            failID1 <- NULL
            failID2 <- NULL
            one <- evaluateDirection(
                otherCriterion$M[otherCriterion$IID == pair[1]],
                otherCriterionTh,
                direction=otherCriterionThDirection)
            two <- evaluateDirection(
                otherCriterion$M[otherCriterion$IID == pair[2]],
                otherCriterionTh,
                direction=otherCriterionThDirection)
            if (one) failID1 <- pair[1]
            if (two) failID2 <- pair[2]
            return(c(failID1, failID2))
        })
        fail_other <- unique(unlist(fail_other))

        # remove samples that fail other criterion
        highRelated <- highRelated[!(highRelated$IID1 %in% fail_other |
                                         highRelated$IID2 %in% fail_other), ]
        uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
    }
    # all samples with related individuals
    allRelated <- c(highRelated$IID1, highRelated$IID2)

    # Further selection of samples with relatives in cohort
    duplicateIDs <- unique(allRelated[duplicated(allRelated)])

    if (length(duplicateIDs) == 0) {
        # Only one related samples per individual
        if (!is.null(otherCriterion)) {
            # choose one with best otherCriterion, if equal arbitrarily return
            # the first one as fail
            failIDs <- apply(highRelated, 1, function(pair) {
                one_two <- evaluateDirection(
                    otherCriterion$M[otherCriterion$IID == pair[1]],
                    otherCriterion$M[otherCriterion$IID == pair[2]],
                    direction=otherCriterionThDirection)
                two_one <- evaluateDirection(
                    otherCriterion$M[otherCriterion$IID == pair[2]],
                    otherCriterion$M[otherCriterion$IID == pair[1]],
                    direction=otherCriterionThDirection)
                if (one_two) failID <- pair[1]
                else if (two_one) failID <- pair[2]
                else failID <- pair[1]
                return(failID)
            })
        } else {
            # randomly choose one to exclude
            failIDs <- highRelated[1,]
        }
    } else {
        # Get all related samples per individual
        relatedPerID <- lapply(duplicateIDs, function(x) {
            tmp <- highRelated[rowSums(cbind(highRelated$IID1 %in% x,
                                             highRelated$IID2 %in% x)) != 0,1:2]
            rel <- unique(unlist(tmp))
            return(rel)
        })
        names(relatedPerID) <- duplicateIDs

        # Find all non-related samples from family structures
        nonRelated <- lapply(relatedPerID, function(rel) {
            pairwise <- t(combn(rel, 2))
            combination <- (pairwise[,1] %in% highRelated$IID1 &
                            pairwise[,2] %in% highRelated$IID2) |
                            (pairwise[,2] %in% highRelated$IID1 &
                            pairwise[,1] %in% highRelated$IID2)
            return(pairwise[!combination,])
        })
        nonRelated <- unique(unlist(nonRelated))
        failIDs <- duplicateIDs[!duplicateIDs %in% nonRelated]
    }
    allFailIIDs <- c(failIDs, fail_other)
    relatednessFails <- lapply(allFailIIDs, function(id) {
        fail_inorder <- relatedness_original$IID1 == id &
            relatedness_original$M > relatednessTh
        fail_inreverse <- relatedness_original$IID2 == id &
            relatedness_original$M > relatednessTh
        if (any(fail_inreverse)) {
            inreverse <- relatedness_original[fail_inreverse, ]
            if (relatednessFID2 %in% names(relatedness_original)) {
                id1 <-
                    c(which(colnames(relatedness_original) == relatednessIID1),
                      which(colnames(relatedness_original) == relatednessFID1))
                id2 <-
                    c(which(colnames(relatedness_original) == relatednessIID2),
                      which(colnames(relatedness_original) == relatednessFID2))
            } else {
                id1 <- which(colnames(relatedness_original) == relatednessIID1)
                id2 <- which(colnames(relatedness_original) == relatednessIID2)
            }
            inreverse[,c(id1, id2)] <- inreverse[,c(id2, id1)]
            names(inreverse) <- relatedness_names
        } else {
            inreverse <- NULL
        }
        inorder <- relatedness_original[fail_inorder, ]
        names(inorder) <- relatedness_names
        return(rbind(inorder, inreverse))
    })
    relatednessFails <- do.call(rbind, relatednessFails)
    if (nrow(relatednessFails) == 0) {
        relatednessFails <- NULL
        failIDs <- NULL
    } else {
        names(relatednessFails) <- relatedness_names
        rownames(relatednessFails) <- 1:nrow(relatednessFails)
        uniqueFails <- relatednessFails[!duplicated(relatednessFails$IID1),]
        if (relatednessFID2 %in% names(relatedness_original)) {
            failIDs <- data.frame(FID=uniqueFails$FID1,
                              IID=uniqueFails$IID1)
        } else {
            failIDs <- data.frame(IID=uniqueFails$IID1)
        }
    }
    return(list(relatednessFails=relatednessFails, failIDs=failIDs))
}

evaluateDirection <- function(x, y, direction) {
    if (direction == 'ge') x >= y
    else if (direction == 'le') x <= y
    else if (direction == 'gt') x > y
    else if (direction == 'lt') x < y
    else if (direction == 'eq') x == y
    else stop(direction, " as direction in evaluateDirection not known.")
}

makepath <- function(directory, name) {
    path <- file.path(directory, name)
    path <- gsub('\\\\', '/', path)
    return(path)
}

checkFormat <- function(prefix) {
    if (!file.exists(paste(prefix, ".fam", sep=""))){
        stop("plink family file: ", prefix, ".fam does not exist.")
    }
    if (!file.exists(paste(prefix, ".bim", sep=""))){
        stop("plink snp file: ", prefix, ".bim does not exist.")
    }
    if (!file.exists(paste(prefix, ".bed", sep=""))){
        stop("plink binary file: ", prefix, ".bed does not exist.")
    }
}
