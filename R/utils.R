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

#' Remove related samples while keeping maximum number of samples in cohort
#'
#' @param ID [N] vector with N samples IDs of study cohort
#' @param relatedness relatedness data.frame as obtained from ukbgene -rel
#' containing pairwise relatedness estimates for individuals in the 500 UKB
#' genotypes. Has to have original header i.e. ID1 ID2 HetHet IBS0 Kinship
#' @return named list with related2filter, related2decide and related2keep.
#' related2filter contains IDs of related samples in the study cohort that
#' should be removed. related2decide contains a dataframe with pairs of related
#' that are only related to one another. One of each pair should be removed -
#' if no other reason is given, simply choose on column as samples to keep and
#' the other as samples to filter. related2keep, contains the IDs
#' of samples in the study cohort that have related individuals in UKB but those
#' are not in the current study - these IDs are just provided for information
#' and nothing needs to be done.

relatednessFilter <- function(ID, relatedness,  otherCriterion=NULL,
                              relatednessTh, otherCriterionTh=NULL,
                              relatednessIID1="IID1", relatednessIID2="IID2",
                              relatednessFID1="FID1", relatednessFID2="FID2",
                              columnRelatedness="PI_HAT"
                              otherCriterionIID="IID",
                              otherCriterionFID="FID",
                              otherCriterionMeasure=NULL ) {

    if (!(columnIID1 %in% names(relatedness))) {
        stop(paste("Column", columnIID1, "not found!"))
    }
    if (!(columnIID2 %in% names(relatedness))) {
        stop(paste("Column", columnIID2, "not found!"))
    }
    if (!(columnRelatedness %in% names(relatedness))) {
        stop(paste("Column", columnRelatedness, "not found!"))
    }
    names(relatedness)[names(relatedness) == columnIID1] <- "IID1"
    names(relatedness)[names(relatedness) == columnIID2] <- "IID2"
    names(relatedness)[names(relatedness) == columnrelatedness] <- "M"

    if (!(is.null(columnFID1) && is.null(columnFID2)) && !is.null(familyInfo)) {
        message("Both, ", columnFID1, "/", columnFID2, "and familyInfo provided",
                "will replace",  columnFID1, "/", columnFID2,
                "with information in familyInfo.")
        if (!(relatednessFID1 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID1, "in relatedness not found!"))
        }
        if (!(relatednessFID2 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID2, "in relatedness not found!"))
        }
        if (!(familyIID %in% names(familyInfo))) {
            stop(paste("Column", familyIID, "in familyInfo not found!"))
        }
        if (!(familyFID %in% names(familyInfo))) {
            stop(paste("Column", familyFID, "in familyInfo not found!"))
        }

        names(relatedness)[names(relatedness) == relatednessFID1] <- "FID1"
        names(relatedness)[names(relatedness) == relatednessFID2] <- "FID2"
        names(familyInfo)[names(familyInfo) == familyFID] <- "FID"
        names(familyInfo)[names(familyInfo) == familyIID] <- "IID"

        familyInfo <- dplyr::select_(familyInfo, ~IID1, ~FID2)
        relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~M)

        familiyIDs <- apply(relatedness, 1, function(pair) {
            FID1 <- ifelse(pair[1] %in% familyInfo$IID, familyInfo$FID, pair[1])
            FID2 <- ifelse(pair[2] %in% familyInfo$IID, familyInfo$FID, pair[2])
            return(c(FID1, FID2))
        })
        relatedness <- data.frame(relatedness, FID1=familiyIDs[,1],
                                  FID2=familiyIDs[,2])
        relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~FID1, ~FID2,
                                      ~M)
    } else if (!(is.null(columnFID1) && is.null(columnFID2))) {
        if (!(columnFID1 %in% names(relatedness))) {
            stop(paste("Column", columnFID1, "not found!"))
        }
        if (!(columnFID2 %in% names(relatedness))) {
            stop(paste("Column", columnFID2, "not found!"))
        }
        names(relatedness)[names(relatedness) == columnFID1] <- "FID1"
        names(relatedness)[names(relatedness) == columnFID2] <- "FID2"
        relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~FID1, ~FID2,
                                      ~M)
    } else {
        relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~relatedness)
    }

    # individuals with at least one pair-wise comparison > relatednessTh
    highRelated <- dplyr::filter_(relatedness, ~M > relatednessTh)
    uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
    uniqueFIDs <- unique(c(highRelated$FID1, highRelated$FID2))

    if (is.null(otherCriterion)) {
        if (!(otherCriterionMeasure %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionMeasure, "not found!"))
        }
        if (!(otherCriterionID %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionID, "not found!"))
        }
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionMeasure] <- "M"
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionID] <- "ID"

        fail_other <- apply(highRelated, 1, function(pair) {
            failID1 <- NULL
            failID2 <- NULL
            if (otherCriterion$M[otherCriterion$ID == pair[1]] >
                otherCriterionTh) {
                failID1 <- pair[1]
            }
            if (otherCriterion$M[otherCriterion$ID == pair[2]] >
                otherCriterionTh) {
                failID2 <- pair[2]
            }
            return(c(failID1, failID2))
        })
        fail_other <- unique(unlist(fail_other))

        # remove samples that fail other criterion
        highRelated <- highRelated[!(highRelated$IID1 %in% fail_other |
                                         highRelated$IID2 %in% fail_other), ]
        uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
        uniqueFIDs <- unique(c(highRelated$FID1, highRelated$FID2))

    }
    # unique IID to FID mapping
    allRelated <- rbind(data.frame(IID=highRelated$IID1, FID=highRelated$FID1,
                               stringsAsFactors=FALSE),
                    data.frame(IID=highRelated$IID2, FID=highRelated$FID2,
                               stringsAsFactors=FALSE))

    # Further selection of samples with relatives in cohort
    duplicateIDs <- unique(allRelated$IID[duplicated(allRelated$IID)])

    if (length(duplicateIDs) == 0) {
        # Only one related samples per individual
        if (!is.null(otherCriterion)) {
            # choose one with best otherCriterion
            failIDs <- apply(highRelated, 1, function(pair) {
                if (otherCriterion$M[otherCriterion$ID == pair[1]] >
                    otherCriterion$M[otherCriterion$ID == pair[2]]) {
                    failID <- pair[1]
                }
                if (otherCriterion$M[otherCriterion$ID == pair[1]] <
                    otherCriterion$M[otherCriterion$ID == pair[2]]) {
                    failID <- pair[2]
                }
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
            return(pairwise[!combination])
        })
        nonRelated <- unique(unlist(nonRelated))
        failIDs <- duplicateIDs[!duplicateIDs %in% nonRelated]
    }
    return(failIDs)
}


