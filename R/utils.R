#' Check PLINK software access
#'
#' checkPlink checks that the PLINK software
#' (\url{https://www.cog-genomics.org/plink/1.9/})
#' can be found from system call.
#'
#' @param path2plink [character] Absolute path to PLINK executable
#' (\url{https://www.cog-genomics.org/plink/1.9/}) i.e.
#' plink should be accessible as path2plink -h. The full name of the executable
#' should be specified: for windows OS, this means path/plink.exe, for unix
#' platforms this is path/plink. If not provided, assumed that PATH set-up works
#' and PLINK will be found by \code{\link[sys]{exec}}('plink').
#' @return Path to PLINK executable.
#' @export

checkPlink <- function(path2plink=NULL) {
    if (is.null(path2plink)) {
        path2plink <- 'plink'
        preset <- FALSE
    } else {
        preset <- TRUE
    }
    if (grepl("~", path2plink)) {
        stop("Path to plink (", path2plink,
             ") contains ~: please supply full path, not relying",
             " on tilde extension.")
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
#' \code{relatednessFilter} takes a data.frame with pair-wise relatedness
#' measures of samples and returns pairs of individual IDs that are related as
#' well as a list of suggested individual IDs to remove.
#' \code{relatednessFilter} finds pairs of samples whose relatedness estimate
#' is larger than the specified relatednessTh. Subsequently, for pairs of
#' individual that do not have additional relatives in the dataset, the
#' individual with the worse otherCriterionMeasure (if provided) or arbitrarily
#' individual 1 of that pair is selected and returned as the individual failing
#' the relatedness check. For more complex family structures, the unrelated
#' individuals per family are selected (e.g. in a simple case of a
#' parents-offspring trio, the offspring will be marked as fail, while the
#' parents will be kept in the analysis). Selection is achieved by constructing
#' subgraphs of clusters of individuals that are related.
#' \code{relatednessFilter} then finds the maximum independent set of vertices
#' in the subgraphs of related individuals. If all individuals are related (i.e.
#' all maximum independent sets are 0), one individual of that cluster will be
#' kept and all others listed as failIDs.
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
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return named [list] with i) relatednessFails, a [data.frame] containing the
#' data.frame relatedness after filtering for pairs of individuals in
#' relatednessIID1 and relatednessIID2, that fail the relatedness
#' QC; the data.frame is reordered with the fail individuals in column 1 and
#' their related individuals in column 2 and ii) failIDs, a [data.frame] with
#' the [IID]s (and [FID]s if provided) of the individuals that fail the
#' relatednessTh.
#' @export

relatednessFilter <- function(relatedness, otherCriterion=NULL,
                              relatednessTh, otherCriterionTh=NULL,
                              otherCriterionThDirection=
                                  c("gt", "ge", "lt", "le", "eq"),
                              relatednessIID1="IID1", relatednessIID2="IID2",
                              relatednessFID1=NULL, relatednessFID2=NULL,
                              relatednessRelatedness="PI_HAT",
                              otherCriterionIID="IID",
                              otherCriterionMeasure=NULL,
                              verbose=FALSE) {

    if (!(relatednessIID1 %in% names(relatedness))) {
        stop(paste("Column", relatednessIID1, "for relatedness not found!"))
    }
    if (!(relatednessIID2 %in% names(relatedness))) {
        stop(paste("Column", relatednessIID1, "for relatedness not found!"))
    }
    if (!(relatednessRelatedness %in% names(relatedness))) {
        stop(paste("Column", relatednessRelatedness,
                   "for relatedness not found!"))
    }

    iid1_index <- which(colnames(relatedness) == relatednessIID1)
    iid2_index <- which(colnames(relatedness) == relatednessIID2)

    relatedness[,iid1_index] <- as.character(relatedness[,iid1_index])
    relatedness[,iid2_index] <- as.character(relatedness[,iid2_index])

    relatedness_names <- names(relatedness)
    names(relatedness)[iid1_index] <- "IID1"
    names(relatedness)[iid2_index] <- "IID2"
    names(relatedness)[names(relatedness) == relatednessRelatedness] <- "M"

    relatedness_original <- relatedness

    if (!is.null(relatednessFID1) && is.null(relatednessFID2) ||
        is.null(relatednessFID1) && !is.null(relatednessFID2) ) {
        stop("Either none or both, relatednessFID1 and relatednessFID2 have to
             be provided")
    }
    if (!is.null(relatednessFID1) && !is.null(relatednessFID2)) {
        if (!(relatednessFID1 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID1, "for relatedness not found!"))
        }
        if (!(relatednessFID2 %in% names(relatedness))) {
            stop(paste("Column", relatednessFID2, "for relatedness not found!"))
        }
        fid1_index <- which(colnames(relatedness) == relatednessFID1)
        fid2_index <- which(colnames(relatedness) == relatednessFID2)
        names(relatedness)[fid1_index] <- "FID1"
        names(relatedness)[fid2_index] <- "FID2"
        relatedness$FID1 <- as.character(relatedness$FID1)
        relatedness$FID2 <- as.character(relatedness$FID2)
    }

    # Remove symmetric IID rows
    relatedness <- dplyr::select_(relatedness, ~IID1, ~IID2, ~M)

    sortedIDs <- data.frame(t(apply(relatedness, 1, function(pair) {
        c(sort(c(pair[1], pair[2])))
        })), stringsAsFactors=FALSE)
    keepIndex <- which(!duplicated(sortedIDs))

    relatedness_original <- relatedness_original[keepIndex,]
    relatedness <- relatedness[keepIndex,]

    # individuals with at least one pair-wise comparison > relatednessTh
    highRelated <- dplyr::filter_(relatedness, ~M > relatednessTh)
    if (nrow(highRelated) == 0) {
        return(list(relatednessFails=NULL, failIDs=NULL))
    }

    uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))

    failIDs_other <- NULL
    if (!is.null(otherCriterion)) {
        otherCriterionThDirection <- match.arg(otherCriterionThDirection)
        if (!(otherCriterionMeasure %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionMeasure,
                       "for otherCriterion not found!"))
        }
        if (!(otherCriterionIID %in% names(otherCriterion))) {
            stop(paste("Column", otherCriterionIID,
                       "for otherCriterion not found!"))
        }
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionMeasure] <- "M"
        names(otherCriterion)[names(otherCriterion) ==
                                  otherCriterionIID] <- "IID"

        if (any(!uniqueIIDs %in% otherCriterion$IID)) {
            stop("Not all IIDs provided in relatedness are contained in" ,
                 "otherCriterion")
        }
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
        failIDs_other <- unique(unlist(fail_other))

        # remove samples that fail other criterion
        highRelated <- highRelated[!(highRelated$IID1 %in% failIDs_other |
                                         highRelated$IID2 %in% failIDs_other), ]
        if (nrow(highRelated) == 0) {
            if (verbose) {
                message("Relatedness cannot be evaluated as all individuals ",
                        "involved fail due to otherCriterion")
            }
            return(list(relatednessFails=NULL, failIDs=NULL))
        }
        uniqueIIDs <- unique(c(highRelated$IID1, highRelated$IID2))
    }
    # all samples with related individuals
    allRelated <- c(highRelated$IID1, highRelated$IID2)

    # Further selection of samples with relatives in cohort
    multipleRelative <- unique(allRelated[duplicated(allRelated)])
    singleRelative <- uniqueIIDs[!uniqueIIDs %in% multipleRelative]

    highRelatedMultiple <- highRelated[highRelated$IID1 %in% multipleRelative |
                                        highRelated$IID2 %in% multipleRelative,]
    highRelatedSingle <- highRelated[highRelated$IID1 %in% singleRelative &
                                       highRelated$IID2 %in% singleRelative,]

    # Only one related samples per individual
    if(length(singleRelative) != 0) {
        if (!is.null(otherCriterion)) {
            # choose one with best otherCriterion, if equal arbitrarily return
            # the first one as fail
            failIDs_single <- apply(highRelatedSingle, 1, function(pair) {
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
            failIDs_single <- highRelatedSingle[,1]
        }
    } else {
        failIDs_single <- NULL
    }

    if(length(multipleRelative) != 0) {
        relatedPerID <- lapply(multipleRelative, function(x) {
            tmp <- highRelatedMultiple[rowSums(
                cbind(highRelatedMultiple$IID1 %in% x,
                      highRelatedMultiple$IID2 %in% x)) != 0,1:2]
            rel <- unique(unlist(tmp))
            return(rel)
        })
        names(relatedPerID) <- multipleRelative

        keepIDs_multiple <- lapply(relatedPerID, function(x) {
            pairwise <- t(combn(x, 2))
            index <- (highRelatedMultiple$IID1 %in% pairwise[,1] &
                          highRelatedMultiple$IID2 %in% pairwise[,2]) |
                (highRelatedMultiple$IID1 %in% pairwise[,2] &
                     highRelatedMultiple$IID2 %in% pairwise[,1])
            combination <- highRelatedMultiple[index,]
            combination_graph <- igraph::graph_from_data_frame(combination,
                                                               directed=FALSE)
            all_iv_set <- igraph::ivs(combination_graph)
            length_iv_set <- sapply(all_iv_set, function(x) length(x))

            if (all(length_iv_set == 1)) {
                # check how often they occurr elsewhere
                occurrence <- sapply(x, function(id) {
                    sum(sapply(relatedPerID, function(idlist) id %in% idlist))
                })
                # if occurrence the same everywhere, pick the first, else keep
                # the one with minimum occurrence elsewhere
                if (length(unique(occurrence)) == 1) {
                    nonRelated <- sort(x)[1]
                } else {
                    nonRelated <- names(occurrence)[which.min(occurrence)]
                }
            } else {
                nonRelated <- all_iv_set[which.max(length_iv_set)]
            }
            return(nonRelated)
        })
        keepIDs_multiple <- unique(unlist(keepIDs_multiple))
        failIDs_multiple <- c(multipleRelative[!multipleRelative %in%
                                                   keepIDs_multiple])
    } else {
        failIDs_multiple <- NULL
    }
    allFailIIDs <- c(failIDs_single, failIDs_multiple, failIDs_other)
    relatednessFails <- lapply(allFailIIDs, function(id) {
        fail_inorder <- relatedness_original$IID1 == id &
            relatedness_original$M > relatednessTh
        fail_inreverse <- relatedness_original$IID2 == id &
            relatedness_original$M > relatednessTh
        if (any(fail_inreverse)) {
            inreverse <- relatedness_original[fail_inreverse, ]
            if (!is.null(relatednessFID2)) {
                id1 <- c(iid1_index,fid1_index)
                id2 <- c(iid2_index,fid2_index)
            } else {
                id1 <- iid1_index
                id2 <- iid2_index
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
        uniqueFails <- relatednessFails[!duplicated(relatednessFails[,iid1_index]),]
        if (!is.null(relatednessFID2)) {
            failIDs <- data.frame(FID=uniqueFails[,fid1_index],
                              IID=uniqueFails[,iid1_index],
                              stringsAsFactors=FALSE)
        } else {
            failIDs <- data.frame(IID=uniqueFails[,iid1_index],
                                  stringsAsFactors=FALSE)
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

#' Check and construct PLINK sample and marker filters
#'
#' checkFiltering checks that the file names with the individuals and markers to
#' be filtered can be found. If so, it constructs the command for filtering
#'
#' @param keep_individuals [character] Path to file with individuals to be
#' retained in the analysis. The file has to be a space/tab-delimited text file
#' with family IDs in the first column and within-family IDs in the second
#' column. All samples not listed in this file will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#indiv}.
#' Default: NULL, i.e. no filtering on individuals.
#' @param remove_individuals [character] Path to file with individuals to be
#' removed from the analysis. The file has to be a space/tab-delimited text file
#' with family IDs in the first column and within-family IDs in the second
#' column. All samples listed in this file will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#indiv}.
#' Default: NULL, i.e. no filtering on individuals.
#' @param extract_markers [character] Path to file with makers to be
#' included in the analysis. The file has to be a text file with a list of
#' variant IDs (usually one per line, but it's okay for them to just be
#' separated by spaces). All unlisted variants will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#snp}.
#' Default: NULL, i.e. no filtering on markers.
#' @param exclude_markers [character] Path to file with makers to be
#' removed from the analysis. The file has to be a text file with a list of
#' variant IDs (usually one per line, but it's okay for them to just be
#' separated by spaces). All listed variants will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#snp}.
#' Default: NULL, i.e. no filtering on markers.
#' @return Vector containing args in sys::exec_wait format to enable filtering
#' on individuals and/or markers.
#' @export
checkFiltering <- function(keep_individuals=NULL,
                           remove_individuals=NULL,
                           extract_markers=NULL,
                           exclude_markers=NULL
                           ) {
    args <- NULL
    if (!is.null(keep_individuals)) {
        if (!file.exists(keep_individuals)) {
            stop("File with individuals to keep in analysis does not exist: ",
                 keep_individuals)
        }
        args = c(args, "--keep", keep_individuals)
    }
    if (!is.null(remove_individuals)) {
        if (!file.exists(remove_individuals)) {
            stop("File with individuals to remove from analysis does not exist: ",
                 remove_individuals)
        }
        args = c(args, "--remove", remove_individuals)
    }
    if (!is.null(extract_markers)) {
        if (!file.exists(extract_markers)) {
            stop("File with markers to extract for analysis does not exist: ",
                 extract_markers)
        }
        args = c(args, "--extract", extract_markers)
    }
    if (!is.null(exclude_markers)) {
        if (!file.exists(exclude_markers)) {
            stop("File with markers to exclude from analysis does not exist: ",
                 exclude_markers)
        }
        args = c(args, "--exclude", exclude_markers)
    }
    return(args)
}


#' Check and construct individual IDs to be removed
#'
#' checkRemoveIDs checks that the file names with the individuals to be filtered
#' can be found. It reads the corresponding files, combines the selected
#' individuals into one data.frame and compares these to all individuals in
#' the analysis.
#' @param prefix [character] Prefix of PLINK files, i.e. path/2/name.bed,
#' path/2/name.bim and path/2/name.fam.
#' @param keep_individuals [character] Path to file with individuals to be
#' retained in the analysis. The file has to be a space/tab-delimited text file
#' with family IDs in the first column and within-family IDs in the second
#' column. All samples not listed in this file will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#indiv}.
#' Default: NULL, i.e. no filtering on individuals.
#' @param remove_individuals [character] Path to file with individuals to be
#' removed from the analysis. The file has to be a space/tab-delimited text file
#' with family IDs in the first column and within-family IDs in the second
#' column. All samples listed in this file will be removed from the current
#' analysis. See \url{https://www.cog-genomics.org/plink/1.9/filter#indiv}.
#' Default: NULL, i.e. no filtering on individuals.
#' @return data.frame containing family (FID) and individual (IID) IDs of
#' individuals to be removed from analysis.
checkRemoveIDs <- function(prefix, remove_individuals=NULL, keep_individuals) {
    removeIDs <- NULL
    allIDs <- data.table::fread(paste(prefix, ".fam", sep=""),
                                data.table=FALSE, stringsAsFactors=FALSE,
                                header=FALSE)
    allIDs <- allIDs[,1:2]

    if (!is.null(remove_individuals)) {
        if (!file.exists(remove_individuals)) {
            stop("File with individuals to remove from analysis does not exist: ",
                 remove_individuals)
        }
        removeIDs <- data.table::fread(remove_individuals, data.table=FALSE,
                                       stringsAsFactors=FALSE,
                                       header=FALSE)
    }
    if (!is.null(keep_individuals)) {
        if (!file.exists(keep_individuals)) {
            stop("File with individuals to keep in analysis does not exist: ",
                 keep_individuals)
        }
        pre_keepIDs <- data.table::fread(keep_individuals, data.table=FALSE,
                                         stringsAsFactors=FALSE,
                                         header=FALSE)
        if(ncol(pre_keepIDs) != 2) {
            stop("File keep_individual is not in the right format; should be ",
                 "two columns separated by space/tab.")
        }
        removeIDs <- rbind(removeIDs,
                           allIDs[!allIDs[,2] %in% pre_keepIDs[,2],])
    }
    return(removeIDs)
}
