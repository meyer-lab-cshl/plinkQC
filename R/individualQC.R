#' Quality control for all individuals in plink-dataset
#'
#' perIndividualQC checks the samples in the plink dataset for their total
#' missingness and heterozygosity rates, the concordance of their assigned sex
#' to their SNP sex, their relatedness to other study individuals and their
#' genetic ancestry.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where results will be saved.
#' Per default, qcdir=indir. If do.evaluate_[analysis] is set to TRUE and
#' do.run_[analysis] is FALSE, \code{perIndividualQC} expects the
#' analysis-specific plink output files in qcdir i.e. do.check_sex expects
#' name.sexcheck, do.evaluate_check_het_and_miss expects name.het and name.imiss,
#' do.evaluate_check_relatedness expects name.genome and name.imiss and
#' Setting do.run_[analysis] TRUE will execute the checks and create the
#' required files. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param dont.check_sex [logical] If TRUE, no sex check will be conducted;
#' short for do.run_check_sex=FALSE and do.evaluate_check_sex=FALSE.
#' Takes precedence over do.run_check_sex and do.evaluate_check_sex.
#' @param dont.check_het_and_miss [logical] If TRUE, no heterozygosity and
#' missingness check will be conducted; short for
#' do.run_check_heterozygosity=FALSE, do.run_check_missingness=FALSE and
#' do.evaluate_check_het_and_miss=FALSE. Takes precedence over
#' do.run_check_heterozygosity, do.run_check_missingness and
#' do.evaluate_check_het_and_miss.
#' @param dont.check_relatedness [logical] If TRUE, no relatedness check will be
#' conducted; short for do.run_check_relatedness=FALSE and
#' do.evaluate_check_relatedness=FALSE. Takes precedence over
#' do.run_check_relatedness and do.evaluate_check_relatedness.
#' @param do.run_check_sex [logical] If TRUE, run \code{\link{run_check_sex}}
#' @param do.run_check_het_and_miss [logical] If TRUE, run
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}}
#' @param do.run_check_relatedness [logical] If TRUE, run
#' \code{\link{run_check_relatedness}}.
#' @param do.evaluate_check_sex [logical] If TRUE, run
#' \code{\link{evaluate_check_sex}}
#' @param do.evaluate_check_het_and_miss [logical] If TRUE, run
#' \code{\link{evaluate_check_het_and_miss}}.
#' @param do.run_superpop_classification [logical] If TRUE, run 
#' \code{\link{superpop_classification}}
#' @param do.evaluate_check_relatedness [logical] If TRUE, run
#' \code{\link{evaluate_check_relatedness}}.
#' @param subplot_label_size [integer] Size of the subplot labeling.
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sampleQC) via ggplot2::ggsave(p=p_sampleQC,
#' other_arguments) or pdf(outfile) print(p_sampleQC) dev.off().
#' If TRUE, i) depicts the X-chromosomal heterozygosity (SNPSEX) of the samples
#' split by their PEDSEX (if do.evaluate_check_sex is TRUE), ii) creates a
#' scatter plot
#' with samples' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis (if do.evaluate_check_het_and_miss is TRUE), and iii) depicts
#' all pair-wise IBD-estimates as histogram (if do.evaluate_check_relatedness is
#' TRUE)
#' .
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @inheritParams checkPlink
#' @inheritParams checkFiltering
#' @inheritParams check_sex
#' @inheritParams check_relatedness
#' @inheritParams check_het_and_miss
#' @inheritParams run_ancestry_format
#' @return Named [list] with i) fail_list, a named [list] with 1.
#' sample_missingness containing a [vector] with sample IIDs failing the
#' missingness threshold imissTh, 2. highIBD containing a [vector] with sample
#' IIDs failing the relatedness threshold highIBDTh, 3. outlying_heterozygosity
#' containing a [vector] with sample IIDs failing the heterozygosity threshold
#' hetTh, 4. mismatched_sex containing a [vector] with the sample IIDs failing
#' the sexcheck based on SNPSEX and femaleTh/maleTh and 5. ancestry containing
#' a dataframe of sample ids and ancestry probablities predicted by a classifier
#' ii) p_sampleQC, a ggplot2-object 'containing' a sub-paneled plot with the
#' QC-plots of \code{\link{check_sex}},
#' \code{\link{check_het_and_miss}}, and
#' \code{\link{check_relatedness}}, which can
#' be shown by print(p_sampleQC).
#' List entries contain NULL if that specific check was not chosen.
#' @details perIndividualQC wraps around the individual QC functions
#' \code{\link{check_sex}}, \code{\link{check_het_and_miss}}, and
#' \code{\link{check_relatedness}}. For details
#' on the parameters and outputs, check these function documentations. For
#' detailed output for fail IIDs (instead of simple IID lists), run each
#' function individually.
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data"
#' 
#' # All quality control checks
#' \dontrun{
#' # whole dataset
#' fail_individuals <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
#' do.run_check_sex=FALSE)
#'
#' # Only check sex and missingness/heterozygosity
#' fail_sex_het_miss <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' dont.check_relatedness=TRUE,
#' interactive=FALSE, verbose=FALSE)
#'
#' # subset of dataset with sample highlighting
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' individual_qc <- perIndividualQC(indir=indir, qcdir=qcdir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' path2plink=path2plink,
#' remove_individuals=remove_individuals_file,
#' highlight_samples=highlight_samples[,2],
#' highlight_type = c("text", "color"), highlight_color="goldenrod")
#' }
perIndividualQC <- function(indir, name, qcdir=indir,
                            dont.check_sex=FALSE,
                            do.run_check_sex=TRUE, do.evaluate_check_sex=TRUE,
                            maleTh=0.8, femaleTh=0.2,
                            externalSex=NULL, externalMale="M",
                            externalSexSex="Sex", externalSexID="IID",
                            externalFemale="F", fixMixup=FALSE,
                            dont.check_het_and_miss=FALSE,
                            do.run_check_het_and_miss=TRUE,
                            do.evaluate_check_het_and_miss=TRUE,
                            imissTh=0.03, hetTh=3,
                            dont.check_relatedness=FALSE,
                            do.run_check_relatedness=TRUE,
                            do.evaluate_check_relatedness=TRUE,
                            highIBDTh=0.1875,
                            mafThRelatedness=0.1,
                            filter_high_ldregion=TRUE,
                            high_ldregion_file=NULL,
                            genomebuild='hg38',
                            label_fail=TRUE,
                            highlight_samples = NULL,
                            highlight_type =
                                c("text", "label", "color", "shape"),
                            highlight_text_size = 3,
                            highlight_color = "#c51b8a",
                            highlight_shape = 17,
                            highlight_legend = FALSE,
                            interactive=FALSE, verbose=TRUE,
                            keep_individuals=NULL,
                            remove_individuals=NULL,
                            exclude_markers=NULL,
                            extract_markers=NULL,
                            legend_text_size = 5,
                            legend_title_size = 7,
                            axis_text_size = 5,
                            axis_title_size = 7,
                            subplot_label_size = 9,
                            title_size = 9,
                            path2plink=NULL, showPlinkOutput=TRUE,
                            path2plink2=NULL, 
                            do.run_superpop_classification=TRUE,
                            path2load_mat = NULL, 
                            plink2format=FALSE,
                            var_format=FALSE
                            ) {

    missing_genotype <- NULL
    highIBD <- NULL
    outlying_heterozygosity <- NULL
    mismatched_sex <- NULL
    ancestry_pred_maj <- NULL
    ancestry_pred_prob <- NULL

    p_sexcheck <- NULL
    p_het_imiss <- NULL
    p_relatedness <- NULL
    p_ancestry <- NULL

    out <- makepath(qcdir, name)

    if (!dont.check_sex) {
        if (do.run_check_sex) {
            run <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
                                 path2plink=path2plink,
                                 showPlinkOutput=showPlinkOutput,
                                 keep_individuals=keep_individuals,
                                 remove_individuals=remove_individuals,
                                 exclude_markers=exclude_markers,
                                 extract_markers=extract_markers,
                                 verbose=verbose)
        }
        if (do.evaluate_check_sex) {
            if (verbose) {
                message("Identification of individuals with discordant sex ",
                        "information")
            }
            fail_sex <- evaluate_check_sex(qcdir=qcdir, indir=indir, name=name,
                                           maleTh=maleTh, femaleTh=femaleTh,
                                           externalSex=externalSex,
                                           externalMale=externalMale,
                                           externalFemale=externalFemale,
                                           externalSexSex=externalSexSex,
                                           externalSexID=externalSexID,
                                           verbose=verbose,
                                           path2plink=path2plink,
                                           showPlinkOutput=showPlinkOutput,
                                           fixMixup=fixMixup,
                                           label_fail=label_fail,
                                           highlight_samples =
                                               highlight_samples,
                                           highlight_type = highlight_type,
                                           highlight_text_size =
                                               highlight_text_size,
                                           highlight_color = highlight_color,
                                           highlight_shape = highlight_shape,
                                           highlight_legend = highlight_legend,
                                           legend_text_size = legend_text_size,
                                           legend_title_size =
                                               legend_title_size,
                                           axis_text_size = axis_text_size,
                                           axis_title_size = axis_title_size,
                                           title_size = title_size,
                                           interactive=FALSE)
            write.table(fail_sex$fail_sex[,1:2],
                        file=paste(out, ".fail-sexcheck.IDs",
                                   sep=""),
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
            if (!is.null(fail_sex$fail_sex) &&
                nrow(fail_sex$fail_sex) != 0) {
                mismatched_sex<- select(fail_sex$fail_sex,
                                        .data$FID, .data$IID)
            }
            if (!is.null(fail_sex$mixup)) {
                write.table(fail_sex$mixup[,1:2],
                            file=paste(out, ".sexcheck_mixup.IDs",
                                       sep=""),
                            quote=FALSE, row.names=FALSE, col.names=FALSE)
            }
            p_sexcheck <- fail_sex$p_sexcheck
        }
    }
    if (!dont.check_het_and_miss) {
        if (do.run_check_het_and_miss) {
            run_miss <- run_check_missingness(qcdir=qcdir, indir=indir,
                                              name=name,
                                              path2plink=path2plink,
                                              showPlinkOutput=showPlinkOutput,
                                              keep_individuals=keep_individuals,
                                              remove_individuals=remove_individuals,
                                              exclude_markers=exclude_markers,
                                              extract_markers=extract_markers,
                                              verbose=verbose)
            run_het <- run_check_heterozygosity(qcdir=qcdir, indir=indir,
                                                name=name,
                                                path2plink=path2plink,
                                                showPlinkOutput=showPlinkOutput,
                                                keep_individuals=keep_individuals,
                                                remove_individuals=remove_individuals,
                                                exclude_markers=exclude_markers,
                                                extract_markers=extract_markers,
                                                verbose=verbose)
        }
        if (do.evaluate_check_het_and_miss) {
            if (verbose) {
                message("Identification of individuals with outlying missing ",
                        "genotype or heterozygosity rates")
            }
            fail_het_imiss <-
                evaluate_check_het_and_miss(qcdir=qcdir, name=name,
                                            imissTh=imissTh,
                                            hetTh=hetTh, label_fail=label_fail,
                                            highlight_samples =
                                                highlight_samples,
                                            highlight_type = highlight_type,
                                            highlight_text_size =
                                                highlight_text_size,
                                            highlight_color = highlight_color,
                                            highlight_shape = highlight_shape,
                                            highlight_legend = highlight_legend,
                                            legend_text_size = legend_text_size,
                                            legend_title_size =
                                                legend_title_size,
                                            axis_text_size = axis_text_size,
                                            axis_title_size = axis_title_size,
                                            title_size = title_size,
                                            interactive=FALSE)
            write.table(fail_het_imiss$fail_imiss[,1:2],
                        file=paste(out, ".fail-imiss.IDs",
                                   sep=""),
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
            if (!is.null(fail_het_imiss$fail_imiss) &&
                nrow(fail_het_imiss$fail_imiss) != 0) {
                missing_genotype <- select(fail_het_imiss$fail_imiss,
                                            .data$FID, .data$IID)
            }

            write.table(fail_het_imiss$fail_het[,1:2],
                        file=paste(out, ".fail-het.IDs",
                                   sep=""),
                        quote=FALSE, row.names=FALSE, col.names=FALSE)
            if (!is.null(fail_het_imiss$fail_het) &&
                nrow(fail_het_imiss$fail_het) != 0) {
                outlying_heterozygosity <- select(fail_het_imiss$fail_het,
                                                  .data$FID, .data$IID)
            } else {
                outlying_heterozygosity <- NULL
            }
            p_het_imiss <- fail_het_imiss$p_het_imiss
        }
    }
    
    
    if (!dont.check_relatedness) {
      if (do.run_check_relatedness) {
        run <- run_check_relatedness(qcdir=qcdir, indir=indir, name=name,
                                     path2plink=path2plink,
                                     mafThRelatedness=mafThRelatedness,
                                     filter_high_ldregion=
                                       filter_high_ldregion,
                                     high_ldregion_file=high_ldregion_file,
                                     genomebuild=genomebuild,
                                     showPlinkOutput=showPlinkOutput,
                                     keep_individuals=keep_individuals,
                                     remove_individuals=remove_individuals,
                                     exclude_markers=exclude_markers,
                                     extract_markers=extract_markers,
                                     verbose=verbose)
      }
      if (do.evaluate_check_relatedness) {
        if (verbose) message("Identification of related individuals")
        fail_relatedness <- evaluate_check_relatedness(qcdir=qcdir,
                                                       name=name,
                                                       imissTh=imissTh,
                                                       highIBDTh=highIBDTh,
                                                       legend_text_size =
                                                         legend_text_size,
                                                       legend_title_size =
                                                         legend_title_size,
                                                       axis_text_size =
                                                         axis_text_size,
                                                       axis_title_size =
                                                         axis_title_size,
                                                       title_size =
                                                         title_size,
                                                       interactive=FALSE)
        write.table(fail_relatedness$failIDs,
                    file=paste(out, ".fail-IBD.IDs", sep=""),
                    row.names=FALSE, quote=FALSE, col.names=FALSE,
                    sep="\t")
        if (!is.null(fail_relatedness$failIDs)  &&
            nrow(fail_relatedness$failIDs) != 0) {
          highIBD <- select(fail_relatedness$failIDs,
                            .data$FID, .data$IID)
        }
        p_relatedness <- fail_relatedness$p_IBD
      }
    }
    
    if (do.run_superpop_classification) {
      if ((plink2format == FALSE) | (var_format == FALSE)) { 
        ancestry_predictions <- run_ancestry_format(indir=indir, name=name, qcdir=qcdir, 
                            verbose=verbose,
                            path2plink2=path2plink2,
                            keep_individuals=keep_individuals,
                            remove_individuals=remove_individuals,
                            exclude_markers=exclude_markers,
                            extract_markers=extract_markers,
                            showPlinkOutput=showPlinkOutput,
                            format = "@:#[hg38]",
                            plink2format = plink2format,
                            var_format = var_format,
                            path2load_mat=path2load_mat, 
                            legend_text_size=legend_text_size,
                            legend_title_size=legend_title_size, 
                            axis_text_size=axis_text_size,
                            axis_title_size=axis_title_size, 
                            title_size=title_size,
                            legend_position = "bottom")
      }
      else { 
        ancestry_predictions <- superpop_classification(indir=indir, name=name, qcdir=qcdir, 
                                                      verbose=verbose,
                                                      path2plink2=path2plink2,
                                                      keep_individuals=keep_individuals,
                                                      remove_individuals=remove_individuals,
                                                      exclude_markers=exclude_markers,
                                                      extract_markers=extract_markers,
                                                      showPlinkOutput=showPlinkOutput,
                                                      path2load_mat=path2load_mat, 
                                                      legend_text_size=legend_text_size,
                                                      legend_title_size=legend_title_size, 
                                                      axis_text_size=axis_text_size,
                                                      axis_title_size=axis_title_size, 
                                                      title_size=title_size,
                                                      legend_position = "bottom")
      }
      ancestry_pred_maj <- ancestry_predictions$prediction_majority
      ancestry_pred_prob <- ancestry_predictions$prediction_prob
      p_ancestry <- ancestry_predictions$p_ancestry
    }

    fail_list <- list(missing_genotype=missing_genotype,
                      highIBD=highIBD,
                      outlying_heterozygosity=outlying_heterozygosity,
                      mismatched_sex=mismatched_sex)

    if(verbose) message(paste("Combine fail IDs into ", out, ".fail.IDs",
                              sep=""))

    uniqueFails <- do.call(rbind, fail_list)
    uniqueFails <- uniqueFails[!duplicated(uniqueFails$IID),]

    write.table(uniqueFails, file=paste(out, ".fail.IDs",sep=""),
                quote=FALSE, row.names=FALSE, col.names=FALSE)

    plots_sampleQC <- list(p_sexcheck=p_sexcheck,
                           p_het_imiss=p_het_imiss,
                           p_relatedness=p_relatedness,
                           p_ancestry=p_ancestry
                           )
    plots_sampleQC <- plots_sampleQC[sapply(plots_sampleQC,
                                            function(x) !is.null(x))]
    subplotLabels <- LETTERS[1:length(plots_sampleQC)]
    
    if (!is.null(p_ancestry)) {
      ancestry_legend <- cowplot::get_legend(p_ancestry)
      plots_sampleQC$p_ancestry <-  plots_sampleQC$p_ancestry +
        theme(legend.position = "None")
      plots_sampleQC$ancestry_legend <- ancestry_legend
      subplotLabels <- c(subplotLabels, "")
    }
    
    if (!is.null(p_sexcheck) && !is.null(p_het_imiss)) {
      first_plots <- cowplot::plot_grid(plotlist=plots_sampleQC[1:2],
                                        nrow=2,
                                        align = "v",
                                        axis = "lr",
                                        labels=subplotLabels[1:2],
                                        label_size = subplot_label_size
      )
      if (!is.null(p_ancestry)) {
        if (!is.null(p_relatedness)) {
          rel_heights <- c(2, 1, 1, 0.5)
          plots_sampleQC <- list(first_plots,
                                 plots_sampleQC[[3]],
                                 plots_sampleQC[[4]],
                                 plots_sampleQC[[5]]
          )
          subplotLabels <- c("", subplotLabels[3:5])
        } else {
          rel_heights <- c(2, 1, 0.5)
          plots_sampleQC <- list(first_plots,
                                 plots_sampleQC[[3]],
                                 plots_sampleQC[[4]]
          )
          subplotLabels <- c("", subplotLabels[3:4])
        }
      } else {
        if (!is.null(p_relatedness)) {
          rel_heights <- c(2, 1)
          plots_sampleQC <- list(first_plots, plots_sampleQC[[3]])
          subplotLabels <- c("", subplotLabels[3])
        } else {
          rel_heights <- 1
          plots_sampleQC <- first_plots
          subplotLabels <- ""
        }
      }
    } else {
      if (!is.null(p_ancestry)) {
        rel_heights <- c(rep(1, length(plots_sampleQC) -1), 0.5)
      } else {
        rel_heights <- c(rep(1, length(plots_sampleQC)))
      }
    }
    p_sampleQC <- cowplot::plot_grid(plotlist=plots_sampleQC,
                                     nrow=length(plots_sampleQC),
                                     labels=subplotLabels,
                                     label_size = subplot_label_size,
                                     rel_heights=rel_heights)
    if (interactive) {
      print(p_sampleQC)
    }
    return(list(fail_list=fail_list, p_sampleQC=p_sampleQC, 
                ancestry_pred_prob=ancestry_pred_prob, 
                ancestry_pred_maj=ancestry_pred_maj))
}

#' Overview of per sample QC
#'
#' \code{overviewPerIndividualQC} depicts results of
#' \code{\link{perIndividualQC}} as intersection plots
#' (via \code{\link[UpSetR]{upset}}) and returns dataframes indicating
#' which QC checks individuals failed or passed.
#'
#'
#' @param results_perIndividualQC [list] Output of \code{\link{perIndividualQC}}
#' i.e. named [list]
#' with i) sample_missingness containing a [vector] with sample IIDs failing
#' the selected missingness threshold imissTh, ii) highIBD containing
#' a [vector] with sample IIDs failing the selected relatedness threshold
#' highIBDTh, iii) outlying_heterozygosity containing a [vector] with sample
#' IIDs failing selected the heterozygosity threshold hetTh, iv) mismatched_sex
#' containing a [vector] with the sample IIDs failing the sexcheck based on
#' SNPSEX and selected femaleTh/maleTh, and
#' v) p_sampleQC, a ggplot2-object 'containing' a sub-paneled plot with the
#' QC-plots of \code{\link{check_sex}},
#' \code{\link{check_het_and_miss}}, and
#' \code{\link{check_relatedness}}.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_overview) via ggplot2::ggsave(p=p_overview,
#' other_arguments) or pdf(outfile) print(p_overview) dev.off().
#' @return Named [list] with i) nr_fail_samples: total number of samples
#' [integer] failing perIndividualQC, ii) fail_QC containing a [data.frame] with
#' samples that failed QC steps (excluding ancestry) with IID, FID,
#' all QC steps applied by perIndividualQC (max=4), with entries=0
#' if passing the QC and entries=1 if failing that particular QC and iii)
#' fail_QC_and_ancestry containing a [data.frame] with samples that failed
#' ancestry and QC checks with IID, FID, QC_fail and
#' Ancestry_fail, with entries=0 if passing and entries=1 if failing that check,
#' iii) p_overview, a ggplot2-object 'containing' a sub-paneled plot with the
#' QC-plots.
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data"
#' \dontrun{
#' fail_individuals <- perIndividualQC(qcdir=qcdir, indir=indir, name=name,
#' refSamplesFile=paste(qcdir, "/HapMap_ID2Pop.txt",sep=""),
#' refColorsFile=paste(qcdir, "/HapMap_PopColors.txt", sep=""),
#' prefixMergedDataset="data.HapMapIII", interactive=FALSE, verbose=FALSE,
#' do.run_check_het_and_miss=FALSE, do.run_check_relatedness=FALSE,
#' do.run_check_sex=FALSE)
#'
#' overview <- overviewPerIndividualQC(fail_individuals)
#' }

overviewPerIndividualQC <- function(results_perIndividualQC,
                                    interactive=FALSE) {
    if (length(perIndividualQC) == 2 &&
        !all(names(results_perIndividualQC) == c("fail_list", "p_sampleQC"))) {
        stop("results_perIndividualQC not direct output of perIndividualQC")
    }
    fail_list <- results_perIndividualQC$fail_list
    samples_fail_all <- do.call(rbind, fail_list)

    # a) overview QC fails independent of ethnicity
    fail_list_wo_ancestry <- fail_list[!names(fail_list) == "ancestry"]
    id_list_wo_ancestry <- sapply(fail_list_wo_ancestry, function(x) x$IID)
    unique_samples_fail_wo_ancestry <- unique(unlist(id_list_wo_ancestry))
    fail_counts_wo_ancestry <- UpSetR::fromList(id_list_wo_ancestry)
    rownames(fail_counts_wo_ancestry) <- unique_samples_fail_wo_ancestry

    p <- UpSetR::upset(fail_counts_wo_ancestry,
                       order.by = "freq",
                       empty.intersections = "on", text.scale=1.2,
                       # Include when UpSetR v1.4.1 is released
                       #title="Overview quality control failures",
                       mainbar.y.label="Samples failing multiple QC checks",
                       sets.x.label="Sample fails per QC check",
                       main.bar.color="#1b9e77", matrix.color="#1b9e77",
                       sets.bar.color="#d95f02")
    p_qc <- cowplot::plot_grid(NULL, p$Main_bar, p$Sizes, p$Matrix,
                               nrow=2, align='v', rel_heights = c(3,1),
                               rel_widths = c(2,3))
    if (interactive) {
        print(p_qc)
    }

    fail_counts_wo_ancestry <- merge(samples_fail_all, fail_counts_wo_ancestry,
                                     by.x=2, by.y=0)
    if ("ancestry" %in% names(fail_list) && !is.null(fail_list$ancestry)) {
        # b) overview of QC and ancestry fails
        fail_all <- list(QC_fail=unique_samples_fail_wo_ancestry,
                         Ancestry_fail=fail_list$ancestry$IID)

        unique_samples_fail_all <- unique(unlist(fail_all))
        fail_counts_all <- UpSetR::fromList(fail_all)
        rownames(fail_counts_all) <- unique_samples_fail_all

        m <- UpSetR::upset(fail_counts_all,
                           order.by = "freq",
                           # Include when UpSetR v1.4.1 is released
                           # title=
                           # "Intersection between QC and ancestry failures",
                           mainbar.y.label=
                               "Samples failing QC and ancestry checks",
                           sets.x.label="Sample fails per QC check",
                           empty.intersections = "on", text.scale=1.2,
                           main.bar.color="#7570b3", matrix.color="#7570b3",
                           sets.bar.color="#e7298a" )
        p_all <- cowplot::plot_grid(NULL, m$Main_bar, m$Sizes, m$Matrix,
                                    nrow=2, align='v', rel_heights = c(3,1),
                                    rel_widths = c(2,3))
        if (interactive) {
            print(p_all)
        }
        fail_counts_all <- merge(samples_fail_all, fail_counts_all,
                                 by.x=2, by.y=0)

        p_overview <- cowplot::plot_grid(NULL, p$Main_bar, p$Sizes, p$Matrix,
                                         NULL, m$Main_bar, m$Sizes, m$Matrix,
                                         nrow=4, align='v',
                                         rel_heights = c(3,1,3,1),
                                         rel_widths = c(2,3))
    } else {
        p_overview <- p_qc
        fail_counts_all <- NULL
    }
    nr_fail_samples <- length(unique(samples_fail_all$IID))
    return(list(nr_fail_samples=nr_fail_samples,
                fail_QC=fail_counts_wo_ancestry,
                fail_QC_and_ancestry=fail_counts_all,
                p_overview=p_overview))
}

#' Identification of individuals with discordant sex information
#'
#' Runs and evaluates results from plink --check-sex.
#' \code{\link{check_sex}} returns IIDs for individuals whose SNPSEX != PEDSEX
#' (where the SNPSEX is determined by the heterozygosity rate across
#' X-chromosomal variants).
#' Mismatching SNPSEX and PEDSEX IDs can indicate plating errors, sample-mixup
#' or generally samples with poor genotyping. In the latter case, these IDs are
#' likely to fail other QC steps as well.
#' Optionally, an extra data.frame (externalSex) with sample IDs and sex can be
#' provided to double check if external and PEDSEX data (often processed at
#' different centers) match. If a mismatch between PEDSEX and SNPSEX was
#' detected, while SNPSEX == Sex, PEDSEX of these individuals can optionally be
#' updated (fixMixup=TRUE).
#' \code{\link{check_sex}} depicts the X-chromosomal heterozygosity (SNPSEX) of
#' the individuals split by their (PEDSEX).
#'
#' \code{\link{check_sex}} wraps around \code{\link{run_check_sex}}  and
#' \code{\link{evaluate_check_sex}}. If run.check_sex is TRUE,
#' \code{\link{run_check_sex}} is executed ; otherwise it is assumed that plink
#' --check-sex has been run externally and qcdir/name.sexcheck exists.
#' \code{\link{check_sex}}  will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_sex, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#sexcheck}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.sexcheck as returned
#' by plink --check-sex. Per default qcdir=indir. If run.check_sex is FALSE, it
#' is assumed that plink --check-sex has been run and qcdir/name.sexcheck is
#' present. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam and name.sexcheck.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param run.check_sex [logical] Should plink --check-sex be run? if set to
#' FALSE, it is assumed that plink --check-sex has been run and
#' qcdir/name.sexcheck is present; \code{\link{check_sex}}  will fail with
#' missing file error otherwise.
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex while Sex==SNPSEX automatically corrected: this will directly
#' change the name.bim/.bed/.fam files!
#' @param externalSex [data.frame, optional] Dataframe with sample IDs
#' [externalSexID] and sex [externalSexSex] to double check if external and
#' PEDSEX data (often processed at different centers) match.
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @param externalSexSex [character] Column identifier for column containing sex
#' information in externalSex.
#' @param externalSexID [character] Column identifier for column containing ID
#' information in externalSex.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sexcheck) via ggplot2::ggsave(p=p_sexcheck,
#' other_arguments) or pdf(outfile) print(p_sexcheck) dev.off().
#' @inheritParams checkPlink
#' @inheritParams run_check_sex
#' @inheritParams evaluate_check_sex
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named list with i) fail_sex: [data.frame] with FID, IID, PEDSEX,
#' SNPSEX and Sex (if externalSex was provided) of individuals failing sex
#' check, ii) mixup: dataframe with FID, IID, PEDSEX, SNPSEX and Sex (if
#' externalSex was provided) of individuals whose PEDSEX != Sex and Sex ==
#' SNPSEX and iii) p_sexcheck, a ggplot2-object 'containing' a scatter plot of
#' the X-chromosomal heterozygosity (SNPSEX) of the sample split by their
#' (PEDSEX), which can be shown by print(p_sexcheck).
#' @export
#' @examples
#'  \dontrun{
#' indir <- system.file("extdata", package="plinkQC")
#' name <- "data"
#'
#' # whole dataset
#' fail_sex <- check_sex(indir=indir, name=name, run.check_sex=FALSE,
#' interactive=FALSE, verbose=FALSE)
#'
#' # subset of dataset with sample highlighting
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' fail_sex <- check_sex(indir=indir, name=name,
#' interactive=FALSE, path2plink=path2plink,
#' remove_individuals=remove_individuals_file,
#' highlight_samples=highlight_samples[,2], highlight_type = c("text", "shape"))
#' }
check_sex <- function(indir, name, qcdir=indir, maleTh=0.8, femaleTh=0.2,
                      run.check_sex=TRUE,
                      externalSex=NULL,
                      externalFemale="F", externalMale="M",
                      externalSexSex="Sex", externalSexID="IID",
                      fixMixup=FALSE,
                      interactive=FALSE, verbose=FALSE,
                      label_fail=TRUE,
                      highlight_samples = NULL,
                      highlight_type =
                          c("text", "label", "color", "shape"),
                      highlight_text_size = 3,
                      highlight_color = "#c51b8a",
                      highlight_shape = 17,
                      highlight_legend = FALSE,
                      path2plink=NULL,
                      keep_individuals=NULL,
                      remove_individuals=NULL,
                      exclude_markers=NULL,
                      extract_markers=NULL,
                      legend_text_size = 5,
                      legend_title_size = 7,
                      axis_text_size = 5,
                      axis_title_size = 7,
                      title_size = 9,
                      showPlinkOutput=TRUE) {
    if (run.check_sex) {
        run_sex <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
                                 verbose=verbose,
                                 path2plink=path2plink,
                                 keep_individuals=keep_individuals,
                                 remove_individuals=remove_individuals,
                                 exclude_markers=exclude_markers,
                                 extract_markers=extract_markers,
                                 showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_sex(qcdir=qcdir, name=name, externalSex=externalSex,
                               maleTh=maleTh, femaleTh=femaleTh,
                               interactive=interactive, fixMixup=fixMixup,
                               indir=indir,
                               label_fail=label_fail,
                               externalFemale=externalFemale,
                               externalMale=externalMale,
                               externalSexSex=externalSexSex,
                               externalSexID=externalSexID,
                               verbose=verbose, path2plink=path2plink,
                               highlight_samples = highlight_samples,
                               highlight_type = highlight_type,
                               highlight_text_size = highlight_text_size,
                               highlight_color = highlight_color,
                               highlight_shape = highlight_shape,
                               highlight_legend = highlight_legend,
                               keep_individuals=keep_individuals,
                               remove_individuals=remove_individuals,
                               exclude_markers=exclude_markers,
                               extract_markers=extract_markers,
                               legend_text_size = legend_text_size,
                               legend_title_size = legend_title_size,
                               axis_text_size = axis_text_size,
                               axis_title_size = axis_title_size,
                               title_size = title_size,
                               showPlinkOutput=showPlinkOutput)
    return(fail)
}


#' Identification of individuals with outlying missing genotype or
#' heterozygosity rates
#'
#' Runs and evaluates results from plink --missing (missing genotype rates
#' per individual) and plink --het (heterozygosity rates per individual).
#' Non-systematic failures in genotyping and outlying heterozygosity (hz) rates
#' per individual are often proxies for DNA sample quality. Larger than expected
#' heterozygosity can indicate possible DNA contamination.
#' The mean heterozygosity in PLINK is computed as hz_mean = (N-O)/N, where
#' N: number of non-missing genotypes and O:observed number of homozygous
#' genotypes for a given individual.
#' Mean heterozygosity can differ between populations and SNP genotyping panels.
#' Within a population and genotyping panel, a reduced heterozygosity rate can
#' indicate inbreeding - these individuals will then likely be returned by
#' \code{\link{check_relatedness}} as individuals that fail the relatedness
#' filters. \code{check_het_and_miss} creates a scatter plot with the
#' individuals' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis.
#'
#' \code{\link{check_het_and_miss}} wraps around
#' \code{\link{run_check_missingness}},
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{evaluate_check_het_and_miss}}.
#' If run.check_het_and_miss is TRUE, \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}} are executed ; otherwise it is assumed
#' that plink --missing and plink --het have been run externally and
#' qcdir/name.het and qcdir/name.imiss exist. \code{\link{check_het_and_miss}}
#' will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_imiss and fail_het, check the
#' original description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#imiss} and
#' \url{https://www.cog-genomics.org/plink/1.9/formats#het}
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory where name.het as returned by
#' plink --het and name.imiss as returned by plink --missing will be saved. Per
#' default qcdir=indir. If run.check_het_and_miss is FALSE, it is assumed that
#' plink --missing and plink --het have been run and qcdir/name.imiss and
#' qcdir/name.het are present. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.het and name.imiss.
#' @param run.check_het_and_miss [logical] Should plink --missing and plink
#' --het be run to determine genotype missingness and heterozygosity rates; if
#' FALSE, it is assumed that plink --missing and plink --het have been run and
#' qcdir/name.imiss and qcdir/name.het are present;
#' \code{\link{check_het_and_miss}} will fail with missing file error otherwise.
#' @param imissTh [double] Threshold for acceptable missing genotype rate per
#' individual; has to be proportion between (0,1)
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity per individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check; has to be
#' larger than 0.
#' @inheritParams checkPlink
#' @inheritParams evaluate_check_het_and_miss
#' @inheritParams run_check_heterozygosity
#' @inheritParams run_check_missingness
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_het_imiss) via ggplot2::ggsave(p=p_het_imiss
#' , other_arguments) or pdf(outfile) print(p_het_imiss) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_imiss [data.frame] containing FID (Family
#' ID), IID (Within-family ID), MISS_PHENO (Phenotype missing? (Y/N)), N_MISS
#' (Number of missing genotype call(s), not including obligatory missings),
#' N_GENO (Number of potentially valid call(s)), F_MISS (Missing call rate) of
#' individuals failing missing genotype check and ii) fail_het [data.frame]
#' containing FID (Family ID), IID (Within-family ID), O(HOM) (Observed number
#' of homozygotes), E(HOM) (Expected number of homozygotes), N(NM) (Number of
#' non-missing autosomal genotypes), F (Method-of-moments F coefficient
#' estimate) of individuals failing  outlying heterozygosity check and iii)
#' p_het_imiss, a ggplot2-object 'containing' a scatter plot with the samples'
#' missingness rates on x-axis and their heterozygosity rates on the y-axis,
#' which can be shown by print(p_het_imiss).
#' @export
#' @examples
#'  \dontrun{
#' indir <- system.file("extdata", package="plinkQC")
#' name <- "data"
#' path2plink <- "path/to/plink"
#'
#' # whole dataset
#' fail_het_miss <- check_het_and_miss(indir=indir, name=name,
#' interactive=FALSE,path2plink=path2plink)
#'
#' # subset of dataset with sample highlighting
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' fail_het_miss <- check_het_and_miss(indir=indir, name=name,
#' interactive=FALSE,path2plink=path2plink,
#' remove_individuals=remove_individuals_file,
#' highlight_samples=highlight_samples[,2], highlight_type = c("text", "shape"))
#' }
check_het_and_miss <- function(indir, name, qcdir=indir, imissTh=0.03, hetTh=3,
                               run.check_het_and_miss=TRUE,
                               label_fail=TRUE,
                               highlight_samples = NULL,
                               highlight_type =
                                   c("text", "label", "color", "shape"),
                               highlight_text_size = 3,
                               highlight_color = "#c51b8a",
                               highlight_shape = 17,
                               highlight_legend = FALSE,
                               interactive=FALSE, verbose=FALSE,
                               keep_individuals=NULL,
                               remove_individuals=NULL,
                               exclude_markers=NULL,
                               extract_markers=NULL,
                               legend_text_size = 5,
                               legend_title_size = 7,
                               axis_text_size = 5,
                               axis_title_size = 7,
                               title_size = 9,
                               path2plink=NULL, showPlinkOutput=TRUE) {
    if (run.check_het_and_miss) {
        run_het <- run_check_heterozygosity(indir=indir,qcdir=qcdir, name=name,
                                            verbose=verbose,
                                            path2plink=path2plink,
                                            keep_individuals=keep_individuals,
                                            remove_individuals=remove_individuals,
                                            exclude_markers=exclude_markers,
                                            extract_markers=extract_markers,
                                            showPlinkOutput=showPlinkOutput)
        run_miss <- run_check_missingness(indir=indir, qcdir=qcdir, name=name,
                                          verbose=verbose,
                                          path2plink=path2plink,
                                          keep_individuals=keep_individuals,
                                          remove_individuals=remove_individuals,
                                          exclude_markers=exclude_markers,
                                          extract_markers=extract_markers,
                                          showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,  hetTh=hetTh,
                                        imissTh=imissTh, interactive=interactive,
                                        label_fail=label_fail,
                                        highlight_samples = highlight_samples,
                                        highlight_type = highlight_type,
                                        highlight_text_size = highlight_text_size,
                                        highlight_color = highlight_color,
                                        highlight_shape = highlight_shape,
                                        highlight_legend = highlight_legend,
                                        legend_text_size = legend_text_size,
                                        legend_title_size = legend_title_size,
                                        title_size = title_size,
                                        axis_text_size = axis_text_size,
                                        axis_title_size = axis_title_size)
    return(fail)
}

#' Identification of related individuals
#'
#' Runs and evaluates results from plink --genome.
#' plink --genome calculates identity by state (IBS) for each pair of
#' individuals based on the average proportion of alleles shared at genotyped
#' SNPs. The degree of recent shared ancestry, i.e. the identity by descent
#' (IBD) can be estimated from the genome-wide IBS. The proportion of IBD
#' between two individuals is returned by plink --genome as PI_HAT.
#' check_relatedness finds pairs of samples whose proportion of IBD is larger
#' than the specified highIBDTh. Subsequently, for pairs of individuals that do
#' not have additional relatives in the dataset, the individual with the greater
#' genotype missingness rate is selected and returned as the individual failing
#' the relatedness check. For more complex family structures, the unrelated
#' individuals per family are selected (e.g. in a parents-offspring trio, the
#' offspring will be marked as fail, while the parents will be kept in the
#' analysis).
#' \code{check_relatedness} depicts all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT.
#'
#' \code{\link{check_relatedness}} wraps around
#' \code{\link{run_check_relatedness}} and
#' \code{\link{evaluate_check_relatedness}}. If run.check_relatedness is TRUE,
#' \code{\link{run_check_relatedness}} is executed ; otherwise it is assumed that
#' plink --genome has been run externally and qcdir/name.genome exists.
#' \code{\link{check_relatedness}}  will fail with missing file error otherwise.
#'
#' For details on the output data.frame fail_high_IBD, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#genome}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to where name.genome as returned
#' by plink --genome will be saved.  Per default qcdir=indir. If
#' run.check_relatedness is FALSE, it is assumed that plink
#' --missing and plink --genome have been run and qcdir/name.imiss and
#' qcdir/name.genome exist. User needs writing permission to qcdir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.genome and name.imiss.
#' @param run.check_relatedness [logical] Should plink --genome be run to
#' determine pairwise IBD of individuals; if FALSE, it is assumed that
#' plink --genome and plink --missing have been run and qcdir/name.imiss and
#' qcdir/name.genome are present;
#' \code{\link{check_relatedness}} will fail with missing file error otherwise.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @inheritParams checkPlink
#' @inheritParams run_check_relatedness
#' @inheritParams evaluate_check_relatedness
#' @inheritParams checkPlink
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_IBD() via ggplot2::ggsave(p=p_IBD,
#' other_arguments) or pdf(outfile) print(p_IBD) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Named [list] with i) fail_high_IBD containing a [data.frame] of
#' IIDs and FIDs of individuals who fail the IBDTh in columns
#' FID1 and IID1. In addition, the following columns are returned (as originally
#' obtained by plink --genome):
#' FID2 (Family ID for second sample), IID2 (Individual ID for second sample),
#' RT (Relationship type inferred from .fam/.ped file), EZ (IBD sharing expected
#' value, based on just .fam/.ped relationship), Z0 (P(IBD=0)), Z1 (P(IBD=1)),
#' Z2 (P(IBD=2)), PI_HAT (Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)), PHE
#' (Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)),
#' DST (IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)), PPC (IBS
#' binomial test), RATIO (HETHET : IBS0 SNP ratio (expected value 2)).
#' and ii) failIDs containing a [data.frame] with individual IDs [IID] and
#' family IDs [FID] of individuals failing the highIBDTh iii) p_IBD, a
#' ggplot2-object 'containing' all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT, which can be
#' shown by print(p_IBD).
#' @export
#' @examples
#' \dontrun{
#' indir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' path2plink <- "path/to/plink"
#'
#' # whole dataset
#' relatednessQC <- check_relatedness(indir=indir, name=name, interactive=FALSE,
#' run.check_relatedness=FALSE, path2plink=path2plink)
#'
#' # subset of dataset
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' fail_relatedness <- check_relatedness(indir=qcdir, name=name,
#' remove_individuals=remove_individuals_file, path2plink=path2plink)
#' }
check_relatedness <- function(indir, name, qcdir=indir, highIBDTh=0.1875,
                              filter_high_ldregion=TRUE,
                              high_ldregion_file=NULL,
                              genomebuild='hg19', imissTh=0.03,
                              run.check_relatedness=TRUE,
                              interactive=FALSE, verbose=FALSE,
                              mafThRelatedness=0.1, path2plink=NULL,
                              keep_individuals=NULL,
                              remove_individuals=NULL,
                              exclude_markers=NULL,
                              extract_markers=NULL,
                              legend_text_size = 5,
                              legend_title_size = 7,
                              axis_text_size = 5,
                              axis_title_size = 7,
                              title_size = 9,
                              showPlinkOutput=TRUE) {
    if (run.check_relatedness) {
        run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name,
                                     verbose=verbose,
                                     mafThRelatedness=mafThRelatedness,
                                     path2plink=path2plink,
                                     highIBDTh=highIBDTh,
                                     genomebuild=genomebuild,
                                     keep_individuals=keep_individuals,
                                     remove_individuals=remove_individuals,
                                     exclude_markers=exclude_markers,
                                     extract_markers=extract_markers,
                                     showPlinkOutput=showPlinkOutput)
    }
    fail <- evaluate_check_relatedness(qcdir=qcdir, name=name,
                                       highIBDTh=highIBDTh,
                                       imissTh=imissTh, interactive=interactive,
                                       legend_text_size = legend_text_size,
                                       legend_title_size = legend_title_size,
                                       axis_text_size = axis_text_size,
                                       axis_title_size = axis_title_size,
                                       title_size = title_size,
                                       verbose=verbose)
    return(fail)
}


#' Run PLINK sexcheck
#'
#' Run  plink --sexcheck to calculate the heterozygosity rate
#' across X-chromosomal variants.
#'
#' Both \code{\link{run_check_sex}} and its evaluation
#' \code{\link{evaluate_check_sex}} can simply be invoked by
#' \code{\link{check_sex}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.sexcheck as returned
#' by plink --check-sex. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @inheritParams checkPlink
#' @inheritParams checkFiltering
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' qcdir <- tempdir()
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # simple sexcheck on all individuals in dataset
#' run <- run_check_sex(indir=indir, qcdir=qcdir, name=name)
#'
#' # sexcheck on subset of dataset
#' keep_individuals_file <- system.file("extdata", "keep_individuals",
#' package="plinkQC")
#' run <- run_check_sex(indir=indir, qcdir=qcdir, name=name,
#' keep_individuals=keep_individuals_file, path2plink=path2plink)
#' }
run_check_sex <- function(indir, name, qcdir=indir, verbose=FALSE,
                          path2plink=NULL,
                          keep_individuals=NULL,
                          remove_individuals=NULL,
                          exclude_markers=NULL,
                          extract_markers=NULL,
                          showPlinkOutput=TRUE) {

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    checkFormat(prefix)
    path2plink <- checkPlink(path2plink)
    args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                  extract_markers, exclude_markers)

    if (verbose) message("Run check_sex via plink --check-sex")
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, "--check-sex", "--out", out,
                          args_filter),
                   std_out=showPlinkOutput, std_err=showPlinkOutput)
}

#' Evaluate results from PLINK sex check.
#'
#' Evaluates and depicts results from plink --check-sex (via
#' \code{\link{run_check_sex}} or externally conducted sex check).
#' Takes file qcdir/name.sexcheck and returns IIDs for samples whose
#' SNPSEX != PEDSEX (where the SNPSEX is determined by the heterozygosity rate
#' across X-chromosomal variants).
#' Mismatching SNPSEX and PEDSEX IDs can indicate plating errors, sample-mixup
#' or generally samples with poor genotyping. In the latter case, these IDs are
#' likely to fail other QC steps as well.
#' Optionally, an extra data.frame (externalSex) with sample IDs and sex can be
#' provided to double check if external and PEDSEX data (often processed at
#' different centers) match. If a mismatch between PEDSEX and SNPSEX was
#' detected while SNPSEX == Sex, PEDSEX of these individuals can optionally be
#' updated (fixMixup=TRUE).
#' \code{evaluate_check_sex} depicts the X-chromosomal heterozygosity (SNPSEX)
#' of the samples split by their (PEDSEX).
#'
#' Both \code{\link{run_check_sex}} and \code{\link{evaluate_check_sex}} can
#' simply be invoked by \code{\link{check_sex}}.
#'
#' For details on the output data.frame fail_sex, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#sexcheck}.
#'
#' @param qcdir [character] /path/to/directory containing name.sexcheck as
#' returned by plink --check-sex.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam and name.sexcheck.
#' @param maleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' males.
#' @param femaleTh [double] Threshold of X-chromosomal heterozygosity rate for
#' females.
#' @param fixMixup [logical] Should PEDSEX of individuals with mismatch between
#' PEDSEX and Sex, with Sex==SNPSEX automatically corrected: this will directly
#' change the name.bim/.bed/.fam files!
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files; only required of fixMixup==TRUE.
#' User needs writing permission to indir.
#' @param externalSex [data.frame, optional] with sample IDs [externalSexID] and
#' sex [externalSexSex] to double check if external and PEDSEX data (often
#' processed at different centers) match.
#' @param externalFemale [integer/character] Identifier for 'female' in
#' externalSex.
#' @param externalMale [integer/character] Identifier for 'male' in externalSex.
#' @param externalSexSex [character] Column identifier for column containing sex
#' information in externalSex.
#' @param externalSexID [character] Column identifier for column containing ID
#' information in externalSex.
#' @param label_fail [logical] Set TRUE, to add fail IDs as text labels in
#' scatter plot.
#' @param highlight_samples [character vector] Vector of sample IIDs to
#' highlight in the plot (p_sexcheck); all highlight_samples IIDs have to
#' be present in the IIDs of the name.fam file.
#' @param highlight_type [character]  Type of sample highlight, labeling
#' by IID ("text"/"label") and/or highlighting data points in different "color"
#' and/or "shape". "text" and "label" use ggrepel for minimal
#' overlap of text labels ("text) or label boxes ("label").  Only one of "text"
#' and "label" can be specified. Text/Label size can be specified with
#' highlight_text_size, highlight color with highlight_color,
#' or highlight shape with highlight_shape.
#' @param highlight_text_size [integer] Text/Label size for samples specified to
#' be highlighted (highlight_samples) by "text" or "label" (highlight_type).
#' @param highlight_color [character] Color for samples specified to
#' be highlighted (highlight_samples) by "color" (highlight_type).
#' @param highlight_shape [integer] Shape for samples specified to
#' be highlighted (highlight_samples) by "shape" (highlight_type). Possible
#' shapes and their encoding can be found at:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#sec:shape-spec}
#' @param highlight_legend [logical] Should a separate legend for the
#' highlighted samples be provided; only relevant for highlight_type == "color"
#' or highlight_type == "shape".
#' @param axis_text_size [integer] Size for axis text.
#' @param axis_title_size [integer] Size for axis title.
#' @param legend_text_size [integer] Size for legend text.
#' @param legend_title_size [integer] Size for legend title.
#' @param title_size [integer] Size for plot title.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_sexcheck) via ggplot2::ggsave(p=p_sexcheck,
#' other_arguments) or pdf(outfile) print(p_sexcheck) dev.off().
#' @inheritParams checkPlink
#' @inheritParams checkFiltering
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return named list with i) fail_sex: dataframe with FID, IID, PEDSEX, SNPSEX
#' and Sex (if externalSex was provided) of individuals failing sex check;
#' ii) mixup: dataframe with FID, IID, PEDSEX, SNPSEX and Sex (if externalSex
#' was provided) of individuals whose PEDSEX != Sex and Sex == SNPSEX; iii)
#' p_sexcheck, a ggplot2-object 'containing' a scatter plot of the X-chromosomal
#' heterozygosity (SNPSEX) of the individuals split by their (PEDSEX), which can
#' be shown by print(p_sexcheck) and iv) plot_data, a data.frame with the data
#' visualised in p_sexcheck (iii).
#' @export
#' @examples
#' qcdir <- system.file("extdata", package="plinkQC")
#' name <- "data"
#' path2plink <- '/path/to/plink'
#' \dontrun{
#' fail_sex <- evaluate_check_sex(qcdir=qcdir, name=name, interactive=FALSE,
#' verbose=FALSE, path2plink=path2plink)
#'
#' # highlight samples
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' fail_sex <- evaluate_check_sex(qcdir=qcdir, name=name, interactive=FALSE,
#' verbose=FALSE, path2plink=path2plink,
#' highlight_samples = highlight_samples[,2],
#' highlight_type = c("label", "color"), highlight_color = "darkgreen")
#' }
evaluate_check_sex <- function(qcdir, name, maleTh=0.8,
                               femaleTh=0.2, externalSex=NULL,
                               fixMixup=FALSE, indir=qcdir,
                               externalFemale="F", externalMale="M",
                               externalSexSex="Sex", externalSexID="IID",
                               verbose=FALSE, label_fail=TRUE,
                               highlight_samples = NULL,
                               highlight_type =
                                   c("text", "label", "color", "shape"),
                               highlight_text_size = 3,
                               highlight_color = "#c51b8a",
                               highlight_shape = 17,
                               highlight_legend = FALSE,
                               legend_text_size = 5,
                               legend_title_size = 7,
                               axis_text_size = 5,
                               axis_title_size = 7,
                               title_size = 9,
                               path2plink=NULL,
                               keep_individuals=NULL,
                               remove_individuals=NULL,
                               exclude_markers=NULL,
                               extract_markers=NULL,
                               showPlinkOutput=TRUE,
                               interactive=FALSE) {

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    if (!file.exists(paste(out, ".sexcheck", sep=""))){
        stop("plink --check-sex results file: ", out,
             ".sexcheck does not exist.")
    }
    testNumerics(numbers=c(maleTh, femaleTh), positives=c(maleTh, femaleTh),
                 proportions=c(maleTh, femaleTh))
    sexcheck <- read.table(paste(out, ".sexcheck",sep=""),
                           header=TRUE, stringsAsFactors=FALSE)

    names_sexcheck <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")
    if (!all(names_sexcheck == names(sexcheck))) {
        stop("Header of", out, ".sexcheck is not correct. Was your
             file generated with plink --check-sex?")
    }
    if (is.null(externalSex)) {
        fail_sex <- sexcheck[sexcheck$STATUS == "PROBLEM",]
        if (nrow(fail_sex) == 0) fail_sex <- NULL
        mixup_geno_pheno <- NULL
    } else {
        if (!(externalSexSex %in% names(externalSex))) {
            stop("Column ", externalSexSex, " not found in externalSex!")
        }
        if (!(externalSexID %in% names(externalSex))) {
            stop("Column ", externalSexID, " not found in externalSex!")
        }
        names(externalSex)[names(externalSex) == externalSexSex] <- "Sex"
        names(externalSex)[names(externalSex) == externalSexID] <- "IID"

        sexcheck_fuse <- merge(sexcheck, externalSex, by="IID")
        sex_mismatch <-
            apply(dplyr::select(sexcheck_fuse, .data$Sex, .data$PEDSEX,
                                .data$SNPSEX), 1,
                  function(ind) {
                      # Mismatch Sex in pheno and geno files
                      if (ind[1] == externalFemale && ind[2] %in% c(0, 1)) {
                          return(ifelse(ind[3] == 2, FALSE, TRUE))
                      }
                      if (ind[1] == externalMale && ind[2] %in% c(0, 2)) {
                          return(ifelse(ind[3] == 1, FALSE, TRUE))
                      }
                      # Mismatch Sex assigned and estimated
                      if (ind[1] == externalFemale && ind[2] == 2) {
                          return(ifelse(ind[3] == 1, TRUE, NA))
                      }
                      if (ind[1] == externalMale && ind[2] == 1) {
                          return(ifelse(ind[3] == 2, TRUE, NA))
                      }
                  })
        # SNPSEX != (Sex in pheno file, PEDSEX)
        fail_sex <-
            dplyr::select(sexcheck_fuse, .data$FID, .data$IID, .data$Sex,
                          .data$PEDSEX, .data$SNPSEX, .data$F)[which(sex_mismatch),]
        if (nrow(fail_sex) == 0) {
            fail_sex <- NULL
            mixup_geno_pheno <- NULL
        } else {
            # SNPSEX == Sex in pheno file and PEDSEX != (SNPSEX, Sex)
            mixup_geno_pheno <-
                dplyr::select(sexcheck_fuse, .data$FID, .data$IID, .data$Sex,
                              .data$PEDSEX, .data$SNPSEX, .data$F)[which(!sex_mismatch),]
            # Fix mismatch between PEDSEX and sex
            if (fixMixup) {
                checkFormat(prefix)
                path2plink <- checkPlink(path2plink)
                args_filter <- checkFiltering(keep_individuals,
                                              remove_individuals,
                                              extract_markers, exclude_markers)
                if (nrow(mixup_geno_pheno) != 0) {
                    file_mixup <- paste(out, ".mismatched_sex_geno_pheno",
                                        sep="")
                    write.table(dplyr::select(mixup_geno_pheno, .data$FID,
                                              .data$IID, .data$SNPSEX),
                                file=file_mixup,
                                row.names=FALSE, quote=FALSE, col.names=FALSE)
                    sys::exec_wait(path2plink, args=c("--bfile", prefix,
                                                      "--update-sex", file_mixup,
                                                      "--make-bed",
                                                      "--out", prefix,
                                                      args_filter),
                                   std_out=showPlinkOutput, std_err=showPlinkOutput)
                } else {
                    if (verbose) {
                        message("All assigned genotype sexes (PEDSEX) match",
                                " external sex assignment (Sex)")
                    }
                    mixup_geno_pheno <- NULL
                }
            } else {
                # Append mismatch to fail_sex, ie treat as fail IDs
                if (nrow(mixup_geno_pheno) != 0) {
                    fail_sex <- rbind(fail_sex, mixup_geno_pheno)
                } else {
                    mixup_geno_pheno <- NULL
                }
            }
        }
    }

    sexcheck$LABELSEX <- "Unassigned"
    sexcheck$LABELSEX[sexcheck$PEDSEX == 1] <- "Male"
    sexcheck$LABELSEX[sexcheck$PEDSEX == 2] <- "Female"

    colors <- c("#999999", "#377eb8", "#e41a1c")
    names(colors) <- c("Unassigned", "Male", "Female")

    sexcheck$shape <- "general"
    shape_guide <- FALSE

    if(!is.null(highlight_samples)) {
        if (!all(highlight_samples %in% sexcheck$IID)) {
            stop("Not all samples to be highlighted are present in the",
                 "name.fam file")
        }
        highlight_type <- match.arg(highlight_type, several.ok = TRUE)
        if (all(c("text", "label") %in% highlight_type)) {
            stop("Only one of text or label highlighting possible; either ",
                 "can be combined with shape and color highlighting")
        }

        if ("shape" %in% highlight_type) {
            sexcheck$shape[sexcheck$IID %in% highlight_samples] <- "highlight"
            shape_guide <- highlight_legend
        }
        if ("color"  %in% highlight_type && highlight_legend) {
            sexcheck$LABELSEX[sexcheck$IID %in% highlight_samples] <- "highlight"
            colors <- c(colors,  highlight_color)
            names(colors)[length(colors)] <- "highlight"
        }
    }


    sexcheck$LABELSEX <- factor(sexcheck$LABELSEX, levels=names(colors))
    sexcheck$PEDSEX <- as.factor(sexcheck$PEDSEX)
    sexcheck$shape <- as.factor(sexcheck$shape)

    p_sexcheck <- ggplot()
    p_sexcheck <- p_sexcheck + geom_point(data=sexcheck,
                                          aes_string(x='PEDSEX', y='F',
                                                     color='LABELSEX',
                                                     shape='shape')) +
        scale_shape_manual(values=c(16, highlight_shape), guide="none") +
        scale_color_manual(values=colors, name="Sex") +
        labs(title="Check assigned sex versus SNP sex",
             x="Reported Sex (PEDSEX)",
             y="ChrX heterozygosity") +
        geom_segment(data=data.frame(x=0.8, xend=1.2, y=maleTh,
                                     yend=maleTh),
                     aes_string(x='x', xend='xend', y='y', yend='yend'), lty=2,
                     color="#e7298a") +
        geom_segment(data=data.frame(x=1.8, xend=2.2, y=femaleTh,
                                     yend=femaleTh), lty=2,
                     aes_string(x='x', xend='xend', y='y', yend='yend'),
                     color="#e7298a")
    if (!is.null(fail_sex) && label_fail) {
        p_sexcheck <- p_sexcheck +
            ggrepel::geom_label_repel(
                data=dplyr::filter(sexcheck, .data$IID %in% fail_sex$IID),
                aes_string(x='PEDSEX',
                           y='F',
                           label='IID'),
                size=highlight_text_size)
    }

    if (!is.null(highlight_samples)) {
        highlight_data <- dplyr::filter(sexcheck, .data$IID %in% highlight_samples)
        if ("text"  %in% highlight_type) {
            p_sexcheck <- p_sexcheck +
                ggrepel::geom_text_repel(data=highlight_data,
                                         aes_string(x='PEDSEX', y='F',
                                                    label="IID"),
                                         size=highlight_text_size)
        }
        if ("label"  %in% highlight_type) {
            p_sexcheck <- p_sexcheck +
                ggrepel::geom_label_repel(data=highlight_data,
                                          aes_string(x='PEDSEX', y='F',
                                                     label='IID'),
                                          size=highlight_text_size)
        }

        if ("color"  %in% highlight_type && !highlight_legend) {
            p_sexcheck <- p_sexcheck +
                geom_point(data=highlight_data,
                           aes_string(x='PEDSEX', y='F', shape='shape'),
                           color=highlight_color)

        }
        if ("shape"  %in% highlight_type && highlight_legend) {
            p_sexcheck <- p_sexcheck +
                guides(shape = "legend") +
                labs(shape = "Individual")
        }
    }
    p_sexcheck <- p_sexcheck +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))

    if (interactive) print(p_sexcheck)
    return(list(fail_sex=fail_sex, mixup=mixup_geno_pheno,
                p_sexcheck=p_sexcheck, plot_data=sexcheck))
}

#' Run PLINK heterozygosity rate calculation
#'
#' Run plink --het to calculate heterozygosity rates per individual.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and their evaluation by
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.het as returned by
#' plink --het. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @inheritParams checkPlink
#' @inheritParams checkFiltering
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' qcdir <- tempdir()
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # heterozygosity check on all individuals in dataset
#' run <- run_check_heterozygosity(indir=indir, qcdir=qcdir, name=name,
#' path2plink=path2plink)
#'
#' #' # heterozygosity on subset of dataset
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' run <- run_check_heterozygosity(indir=indir, qcdir=qcdir, name=name,
#' remove_individuals=remove_individuals_file,path2plink=path2plink)

#' }
run_check_heterozygosity <- function(indir, name, qcdir=indir, verbose=FALSE,
                                     path2plink=NULL,
                                     keep_individuals=NULL,
                                     remove_individuals=NULL,
                                     exclude_markers=NULL,
                                     extract_markers=NULL,
                                     showPlinkOutput=TRUE) {

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    checkFormat(prefix)
    path2plink <- checkPlink(path2plink)
    args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                  extract_markers, exclude_markers)

    if (verbose) message("Run check_heterozygosity via plink --het")
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, "--het", "--out", out,
                          args_filter),
                   std_out=showPlinkOutput, std_err=showPlinkOutput)
}

#' Run PLINK missingness rate calculation
#'
#' Run  plink --missing to calculate missing genotype rates
#' per individual.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and their evaluation by
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.imiss as returned
#' by plink --missing. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @inheritParams checkPlink
#' @inheritParams checkFiltering
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' qcdir <- tempdir()
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # missingness check on all individuals in dataset
#' run <- run_check_missingness(indir=indir, qcdir=qcdir, name=name,
#' path2plink=path2plink)
#'
#' # missingness on subset of dataset
#' remove_individuals_file <- system.file("extdata", "remove_individuals",
#' package="plinkQC")
#' run <- run_check_missingness(indir=indir, qcdir=qcdir, name=name,
#' remove_individuals=remove_individuals_file, path2plink=path2plink)
#' }
run_check_missingness <- function(indir, name, qcdir=indir, verbose=FALSE,
                                  path2plink=NULL,
                                  keep_individuals=NULL,
                                  remove_individuals=NULL,
                                  exclude_markers=NULL,
                                  extract_markers=NULL,
                                  showPlinkOutput=TRUE) {

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    checkFormat(prefix)
    path2plink <- checkPlink(path2plink)
    args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                  extract_markers, exclude_markers)

    if (verbose) message("Run check_missingness via plink --missing")
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, "--missing", "--out", out,
                          args_filter),
                   std_out=showPlinkOutput, std_err=showPlinkOutput)
}

#' Evaluate results from PLINK missing genotype and heterozygosity rate check.
#'
#' Evaluates and depicts results from plink --missing (missing genotype rates
#' per individual) and plink --het (heterozygosity rates per individual) via
#' \code{\link{run_check_heterozygosity}} and
#' \code{\link{run_check_missingness}} or externally conducted check.)
#' Non-systematic failures in genotyping and outlying heterozygosity (hz) rates
#' per individual are often proxies for DNA sample quality. Larger than expected
#' heterozygosity can indicate possible DNA contamination.
#' The mean heterozygosity in PLINK is computed as hz_mean = (N-O)/N, where
#' N: number of non-missing genotypes and O:observed number of homozygous
#' genotypes for a given individual.
#' Mean heterozygosity can differ between populations and SNP genotyping panels.
#' Within a population and genotyping panel, a reduced heterozygosity rate can
#' indicate inbreeding - these individuals will then be returned by
#' check_relatedness as individuals that fail the relatedness filters.
#' \code{evaluate_check_het_and_miss} creates a scatter plot with the
#' individuals' missingness rates on x-axis and their heterozygosity rates on
#' the y-axis.
#'
#' All, \code{\link{run_check_heterozygosity}},
#' \code{\link{run_check_missingness}} and
#' \code{\link{evaluate_check_het_and_miss}} can simply be invoked by
#' \code{\link{check_het_and_miss}}.
#'
#' For details on the output data.frame fail_imiss and fail_het, check the
#' original description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#imiss} and
#' \url{https://www.cog-genomics.org/plink/1.9/formats#het}
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' name.imiss and name.het results as returned by plink --missing and plink
#' --het.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.het and name.imiss.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @param hetTh [double] Threshold for acceptable deviation from mean
#' heterozygosity in any individual. Expressed as multiples of standard
#' deviation of heterozygosity (het), i.e. individuals outside mean(het) +/-
#' hetTh*sd(het) will be returned as failing heterozygosity check; has to be
#' larger than 0.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_het_imiss) via ggplot2::ggsave(p=p_het_imiss
#' , other_arguments) or pdf(outfile) print(p_het_imiss) dev.off().
#' @param label_fail [logical] Set TRUE, to add fail IDs as text labels in
#' scatter plot.
#' @param highlight_samples [character vector] Vector of sample IIDs to
#' highlight in the plot (p_het_imiss); all highlight_samples IIDs have to
#' be present in the IIDs of the name.fam file.
#' @param highlight_type [character]  Type of sample highlight, labeling
#' by IID ("text"/"label") and/or highlighting data points in different "color"
#' and/or "shape". "text" and "label" use ggrepel for minimal
#' overlap of text labels ("text) or label boxes ("label").  Only one of "text"
#' and "label" can be specified.Text/Label size can be specified with
#' highlight_text_size, highlight color with highlight_color,
#' or highlight shape with highlight_shape.
#' @param highlight_text_size [integer] Text/Label size for samples specified to
#' be highlighted (highlight_samples) by "text" or "label" (highlight_type).
#' @param highlight_color [character] Color for samples specified to
#' be highlighted (highlight_samples) by "color" (highlight_type).
#' @param highlight_shape [integer] Shape for samples specified to
#' be highlighted (highlight_samples) by "shape" (highlight_type). Possible
#' shapes and their encoding can be found at:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#sec:shape-spec}
#' @param highlight_legend [logical] Should a separate legend for the
#' highlighted samples be provided; only relevant for highlight_type == "color"
#' or highlight_type == "shape".
#' @param axis_text_size [integer] Size for axis text.
#' @param axis_title_size [integer] Size for axis title.
#' @param legend_text_size [integer] Size for legend text.
#' @param legend_title_size [integer] Size for legend title.
#' @param title_size [integer] Size for plot title.
#' @return named [list] with i) fail_imiss dataframe containing FID (Family ID),
#' IID (Within-family ID), MISS_PHENO (Phenotype missing? (Y/N)), N_MISS (Number
#' of missing genotype call(s), not including obligatory missings), N_GENO (
#' Number of potentially valid call(s)), F_MISS (Missing call rate) of
#' individuals failing missing genotype check and ii) fail_het dataframe
#' containing FID (Family ID), IID (Within-family ID), O(HOM) (Observed number
#' of homozygotes), E(HOM) (Expected number of homozygotes), N(NM) (Number of
#' non-missing autosomal genotypes), F (Method-of-moments F coefficient
#' estimate) of individuals failing  outlying heterozygosity check; iii)
#' p_het_imiss, a ggplot2-object 'containing' a scatter plot with the samples'
#' missingness rates on x-axis and their heterozygosity rates on the y-axis,
#' which can be shown by print(p_het_imiss) and iv) plot_data, a data.frame with
#' the data visualised in p_het_imiss (iii).
#' @export
#' @examples
#' qcdir <- system.file("extdata", package="plinkQC")
#' name <- "data"
#' \dontrun{
#' fail_het_miss <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,
#' interactive=FALSE)
#'
#' #' # highlight samples
#' highlight_samples <- read.table(system.file("extdata", "keep_individuals",
#' package="plinkQC"))
#' fail_het_miss <- evaluate_check_het_and_miss(qcdir=qcdir, name=name,
#' interactive=FALSE, highlight_samples = highlight_samples[,2],
#' highlight_type = c("text", "color"))
#' }
evaluate_check_het_and_miss <- function(qcdir, name, imissTh=0.03,
                                        hetTh=3, label_fail=TRUE,
                                        highlight_samples = NULL,
                                        highlight_type =
                                            c("text", "label", "color", "shape"),
                                        highlight_text_size = 3,
                                        highlight_color = "#c51b8a",
                                        highlight_shape = 17,
                                        legend_text_size = 5,
                                        legend_title_size = 7,
                                        axis_text_size = 5,
                                        axis_title_size = 7,
                                        title_size = 9,
                                        highlight_legend = FALSE,
                                        interactive=FALSE) {

    prefix <- makepath(qcdir, name)

    if (!file.exists(paste(prefix, ".imiss",sep=""))){
        stop("plink --missing output file: ", prefix,
             ".imiss does not exist.")
    }
    if (!file.exists(paste(prefix, ".het",sep=""))){
        stop("plink --het output file: ", prefix,
             ".het does not exist.")
    }
    testNumerics(numbers=c(imissTh, hetTh), positives=c(imissTh, hetTh),
                 proportions=imissTh)
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)
    if (!all(names_imiss == names(imiss))) {
        stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    fail_imiss <- imiss[imiss$F_MISS > imissTh,]

    names_het <- c("FID", "IID", "O.HOM.", "E.HOM.", "N.NM.", "F")
    het <- read.table(paste(prefix, ".het", sep=""), header=TRUE, as.is=TRUE)
    if (!all(names_het == names(het))) {
        stop("Header of ", prefix, ".het is not correct. Was your
             file generated with plink --het?")
    }
    fail_het <- het[het$F < (mean(het$F)  - hetTh*sd(het$F)) |
                        het$F > (mean(het$F) + hetTh*sd(het$F)),]


    nr_samples <- nrow(imiss)
    imiss$logF_MISS <- log10(imiss$F_MISS)
    het_imiss <- merge(imiss, het, by="IID")
    fail_het_imiss <- het_imiss[which(het_imiss$IID %in%
                                          union(fail_het$IID, fail_imiss$IID)),]
    if (nrow(fail_het_imiss) == 0) {
        fail_het_imiss <- NULL
    }
    het_imiss$type <- "pass"
    het_imiss$type[het_imiss$IID %in% fail_het$IID] <- "fail het"
    het_imiss$type[het_imiss$IID %in% fail_imiss$IID] <- "fail miss"
    het_imiss$type[het_imiss$IID %in%
                       intersect(fail_het$IID, fail_imiss$IID)] <- "fail het + miss"

    minus_sd <- mean(het_imiss$F) - 1:5*(sd(het_imiss$F))
    plus_sd <- mean(het_imiss$F) + 1:5*(sd(het_imiss$F))

    colors <- c("#666666", "#1b9e77", "#d95f02", "#7570b3")
    names(colors) <- c("pass", "fail het", "fail miss", "fail het + miss" )

    het_imiss$shape <- "general"
    shape_guide <- FALSE

    if(!is.null(highlight_samples)) {
        if (!all(highlight_samples %in% het_imiss$IID)) {
            stop("Not all samples to be highlighted are present in the",
                 "prefixMergedDataset")
        }
        highlight_type <- match.arg(highlight_type, several.ok = TRUE)
        if (all(c("text", "label") %in% highlight_type)) {
            stop("Only one of text or label highlighting possible; either ",
                 "can be combined with shape and color highlighting")
        }

        if ("shape" %in% highlight_type) {
            het_imiss$shape[het_imiss$IID %in% highlight_samples] <- "highlight"
            shape_guide <- highlight_legend
        }
        if ("color"  %in% highlight_type && highlight_legend) {
            het_imiss$type[het_imiss$IID %in% highlight_samples] <- "highlight"
            colors <- c(colors,  highlight_color)
            names(colors)[length(colors)] <- "highlight"
        }
    }

    het_imiss$type <- factor(het_imiss$type, levels=names(colors))
    het_imiss$shape <- as.factor(het_imiss$shape)

    p_het_imiss <- ggplot()
    p_het_imiss <- p_het_imiss + geom_point(data=het_imiss,
                                            aes_string(x='logF_MISS', y='F',
                                                       color='type',
                                                       shape="shape")) +
        scale_shape_manual(values=c(16, highlight_shape), guide="none") +
        scale_color_manual(values=colors) +
        labs(x = "Proportion of missing SNPs",
             y = "heterozygosity rate (and sd)",
             color = "Marker",
             title = "heterozygosity by Missingness across samples") +
        geom_hline(yintercept=c(minus_sd[1:3], plus_sd[1:3]), lty=2,
                   col="azure4") +
        scale_y_continuous(labels=c("-5", "-4", "-3" ,"+3", "+4", "+5"),
                           breaks=c(minus_sd[3:5], plus_sd[3:5])) +
        scale_x_continuous(labels=c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1),
                           breaks=c(-4,-3,-2, log10(0.03), log10(0.05),-1,0)) +
        geom_hline(yintercept=mean(het_imiss$F) - (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_hline(yintercept=mean(het_imiss$F) + (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_vline(xintercept=log10(imissTh), col="#e7298a", lty=2)

    if (!is.null(fail_het_imiss) && label_fail) {
        p_het_imiss <-
            p_het_imiss + ggrepel::geom_label_repel(
                data=data.frame(x=fail_het_imiss$logF_MISS,
                                y=fail_het_imiss$F,
                                label=fail_het_imiss$IID),
                aes_string(x='x', y='y', label='label'),
                size=highlight_text_size)
    }
    highlight_data <- dplyr::filter(het_imiss, .data$IID %in% highlight_samples)
    if (!is.null(highlight_samples)) {
        if ("text"  %in% highlight_type) {
            p_het_imiss <- p_het_imiss +
                ggrepel::geom_text_repel(data=highlight_data,
                                         aes_string(x='logF_MISS', y='F',
                                                    label="IID"),
                                         size=highlight_text_size)
        }
        if ("label"  %in% highlight_type) {
            p_het_imiss <- p_het_imiss +
                ggrepel::geom_label_repel(data=highlight_data,
                                          aes_string(x='logF_MISS', y='F',
                                                     label="IID"),
                                          size=highlight_text_size)
        }

        if ("color"  %in% highlight_type && !highlight_legend) {
            p_het_imiss <- p_het_imiss +
                geom_point(data=highlight_data,
                           aes_string(x='logF_MISS', y='F', shape='shape'),
                           color=highlight_color,
                           show.legend=highlight_legend)
        }
        if ("shape"  %in% highlight_type && highlight_legend) {
            p_het_imiss <- p_het_imiss +
                guides(shape = "legend")
        }
    }
    p_het_imiss <- p_het_imiss +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              title = element_text(size = title_size),
              axis.title = element_text(size = axis_title_size))
    if (interactive) print(p_het_imiss)
    return(list(fail_imiss=fail_imiss, fail_het=fail_het,
                p_het_imiss=p_het_imiss,
                plot_data = het_imiss))
}

#' Run PLINK IBD estimation
#'
#' Run LD pruning on dataset with plink --exclude range highldfile
#' --indep-pairwise 50 5 0.2, where highldfile contains regions of high LD as
#' provided by Anderson et (2010) Nature Protocols. Subsequently, plink
#' --genome is run on the LD pruned, maf-filtered data. plink --genome
#' calculates identity by state (IBS) for each pair of individuals based on the
#' average proportion of alleles shared at genotyped SNPs. The degree of recent
#' shared ancestry,i.e. the identity by descent (IBD) can be estimated from the
#' genome-wide IBS. The proportion of IBD between two individuals is returned by
#' --genome as PI_HAT.
#'
#' Both \code{\link{run_check_relatedness}} and its evaluation via
#' \code{\link{evaluate_check_relatedness}} can simply be invoked by
#' \code{\link{check_relatedness}}.
#'
#' @param indir [character] /path/to/directory containing the basic PLINK data
#' files name.bim, name.bed, name.fam files.
#' @param qcdir [character] /path/to/directory to save name.genome as returned
#' by plink --genome. User needs writing permission to qcdir. Per default
#' qcdir=indir.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals; only pairwise relationship estimates larger than this
#' threshold will be recorded.
#' @param mafThRelatedness [double] Threshold of minor allele frequency filter
#' for selecting variants for IBD estimation.
#' @param filter_high_ldregion [logical] Should high LD regions be filtered
#'   before IBD estimation; carried out per default with high LD regions for
#'   hg19 provided as default via \code{genomebuild}. For alternative genome
#'   builds not provided or non-human data, high LD regions files can be
#'   provided via \code{high_ldregion_file}.
#' @param high_ldregion_file [character] Path to file with high LD regions used
#' for filtering before IBD estimation if \code{filter_high_ldregion} == TRUE,
#' otherwise ignored; for human genome data, high LD region files are provided
#' and can simply be chosen via \code{genomebuild}. Files have to be
#' space-delimited, no column names with the following columns: chromosome,
#' region-start, region-end, region number. Chromosomes are specified without
#' 'chr' prefix. For instance:
#' 1 48000000 52000000 1
#' 2 86000000 100500000 2
#' @param genomebuild [character] Name of the genome build of the PLINK file
#' annotations, ie mappings in the name.bim file. Will be used to remove
#' high-LD regions based on the coordinates of the respective build. Options
#' are hg18, hg19 and hg38. See @details.
#' @inheritParams checkFiltering
#' @inheritParams checkPlink
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @details The IBD estimation is conducted on LD pruned data and in a first
#' step, high LD regions are excluded. The regions were derived from the
#' high-LD-regions file provided by Anderson et (2010) Nature Protocols. These
#' regions are in NCBI36 (hg18) coordinates and were lifted to GRCh37 (hg19)
#' and GRC38 (hg38) coordinates using the liftOver tool available here:
#' \url{https://genome.ucsc.edu/cgi-bin/hgLiftOver}. The 'Minimum ratio of bases
#' that must remap' which was set to 0.5 and the 'Allow multiple output regions'
#' box ticked; for all other parameters, the default options were selected.
#' LiftOver files were generated on July 9,2019. The commands for formatting
#' the files are provided in system.file("extdata", 'liftOver.cmd',
#' package="plinkQC").
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' qcdir <- tempdir()
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' # Relatedness estimation based in all markers in dataset
#' run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name,
#' path2plink=path2plink)
#'
#' # relatedness estimation on subset of dataset
#' keep_individuals_file <- system.file("extdata", "keep_individuals",
#' package="plinkQC")
#' run <- run_check_relatedness(indir=indir, qcdir=qcdir, name=name,
#' keep_individuals=keep_individuals_file, path2plink=path2plink)
#' }
run_check_relatedness <- function(indir, name, qcdir=indir, highIBDTh=0.185,
                                  mafThRelatedness=0.1,
                                  path2plink=NULL,
                                  filter_high_ldregion=TRUE,
                                  high_ldregion_file=NULL,
                                  genomebuild='hg19',
                                  showPlinkOutput=TRUE,
                                  keep_individuals=NULL,
                                  remove_individuals=NULL,
                                  exclude_markers=NULL,
                                  extract_markers=NULL,
                                  verbose=FALSE) {

    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)

    checkFormat(prefix)
    path2plink <- checkPlink(path2plink)
    args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                  extract_markers, exclude_markers)
    if (filter_high_ldregion) {
        if (!is.null(high_ldregion_file)) {
            if (!file.exists(high_ldregion_file)) {
                stop("high_ldregion_file (", high_ldregion_file ,
                     ") cannot be read")
            }
            highld <- data.table::fread(high_ldregion_file, sep=" ",
                                        header=FALSE, data.table=FALSE)
            if(ncol(highld) != 4) {
                stop("high_ldregion_file (", high_ldregion_file ,
                    ") is incorrectly formated: ",
                    "contains more/less than 4 columns")
            }
            if(any(grepl("chr", highld[,1]))) {
                stop("high_ldregion_file (", high_ldregion_file ,
                     ") is incorrectly formated: ",
                     "chromosome specification in first column",
                     "cannot contain 'chr'")
            }
            if (verbose) message(paste("Using", high_ldregion_file, "coordinates for pruning of",
                                       "high-ld regions"))
        } else {
            if (tolower(genomebuild) == 'hg18' || tolower(genomebuild) == 'NCBI36') {
                high_ldregion_file <- system.file("extdata", 'high-LD-regions-hg18-NCBI36.txt',
                                      package="plinkQC")
            } else if (tolower(genomebuild) == 'hg19' ||
                       tolower(genomebuild) == 'grch37') {
                high_ldregion_file <- system.file("extdata", 'high-LD-regions-hg19-GRCh37.txt',
                                      package="plinkQC")
            } else if (tolower(genomebuild) == 'hg38' ||
                       tolower(genomebuild) == 'grch38') {
                high_ldregion_file <- system.file("extdata", 'high-LD-regions-hg38-GRCh38.txt',
                                      package="plinkQC")
            } else {
                stop(genomebuild, "is not a known/provided human genome build.",
                     "Options are: hg18, hg19, and hg38")
            }
            if (verbose) message(paste("Use", genomebuild, "coordinates for pruning of",
                                       "high-ld regions"))
        }
        if (verbose) message(paste("Prune", prefix, "for relatedness estimation"))
        sys::exec_wait(path2plink,
                       args=c("--bfile", prefix,
                              "--exclude", "range", high_ldregion_file,
                              "--indep-pairwise", 50, 5, 0.2,
                              "--out", out,
                              args_filter),
                       std_out=showPlinkOutput, std_err=showPlinkOutput)
    } else {
        if (verbose) message("No pruning of high-ld regions")
        if (verbose) message(paste("Prune", prefix,
                                   "for relatedness estimation"))
        sys::exec_wait(path2plink,
                       args=c("--bfile", prefix,
                              "--indep-pairwise", 50, 5, 0.2, "--out", out,
                              args_filter),
                       std_out=showPlinkOutput, std_err=showPlinkOutput)
    }

    if (verbose) message("Run check_relatedness via plink --genome")
    if (!is.null(mafThRelatedness)) {
        maf <- c("--maf", mafThRelatedness)
    } else {
        maf <- NULL
    }
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, "--extract",
                          paste(out, ".prune.in", sep=""),
                          maf, "--genome",
                          # "--min", highIBDTh,
                          "--out", out,
                          args_filter),
                   std_out=showPlinkOutput, std_err=showPlinkOutput)
    if (!file.exists(paste(prefix, ".imiss", sep=""))) {
        sys::exec_wait(path2plink,
                       args=c("--bfile", prefix, "--missing", "--out", out,
                              args_filter
                              ),
                       std_out=showPlinkOutput, std_err=showPlinkOutput)
    }
}

#' Evaluate results from PLINK IBD estimation.
#'
#' Evaluates and depicts results from plink --genome on the LD pruned dataset
#' (via \code{\link{run_check_relatedness}} or externally conducted IBD
#' estimation). plink --genome calculates identity by state (IBS) for each pair
#' of individuals based on the average proportion of alleles shared at genotyped
#' SNPs. The degree of recent shared ancestry, i.e. the identity by descent
#' (IBD) can be estimated from the genome-wide IBS. The proportion of IBD
#' between two individuals is returned by --genome as PI_HAT.
#' \code{evaluate_check_relatedness} finds pairs of samples whose proportion of
#' IBD is larger than the specified highIBDTh. Subsequently, for pairs of
#' individual that do not have additional relatives in the dataset, the
#' individual with the greater genotype missingness rate is selected and
#' returned as the individual failing the relatedness check. For more complex
#' family structures, the unrelated individuals per family are selected (e.g. in
#' a parents-offspring trio, the offspring will be marked as fail, while the
#' parents will be kept in the analysis).
#' \code{evaluate_check_relatedness} depicts all pair-wise IBD-estimates as
#' histograms stratified by value of PI_HAT.
#'
#' Both \code{\link{run_check_relatedness}} and
#' \code{\link{evaluate_check_relatedness}} can simply be invoked by
#' \code{\link{check_relatedness}}.
#'
#' For details on the output data.frame fail_high_IBD, check the original
#' description on the PLINK output format page:
#' \url{https://www.cog-genomics.org/plink/1.9/formats#genome}.
#'
#' @param qcdir [character] path/to/directory/with/QC/results containing
#' name.imiss and name.genome results as returned by plink --missing and plink
#' --genome.
#' @param name [character] Prefix of PLINK files, i.e. name.bed, name.bim,
#' name.fam, name.genome and name.imiss.
#' @param highIBDTh [double] Threshold for acceptable proportion of IBD between
#' pair of individuals.
#' @param imissTh [double] Threshold for acceptable missing genotype rate in any
#' individual; has to be proportion between (0,1)
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_IBD() via ggplot2::ggsave(p=p_IBD,
#' other_arguments) or pdf(outfile) print(p_IBD) dev.off().
#' @param axis_text_size [integer] Size for axis text.
#' @param axis_title_size [integer] Size for axis title.
#' @param legend_text_size [integer] Size for legend text.
#' @param legend_title_size [integer] Size for legend title.
#' @param title_size [integer] Size for plot title.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return a named [list] with i) fail_high_IBD containing a [data.frame] of
#' IIDs and FIDs of individuals who fail the IBDTh in columns
#' FID1 and IID1. In addition, the following columns are returned (as originally
#' obtained by plink --genome):
#' FID2 (Family ID for second sample), IID2 (Individual ID for second sample),
#' RT (Relationship type inferred from .fam/.ped file), EZ (IBD sharing expected
#' value, based on just .fam/.ped relationship), Z0 (P(IBD=0)), Z1 (P(IBD=1)),
#' Z2 (P(IBD=2)), PI_HAT (Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)), PHE
#' (Pairwise phenotypic code (1, 0, -1 = AA, AU, and UU pairs, respectively)),
#' DST (IBS distance, i.e. (IBS2 + 0.5*IBS1) / (IBS0 + IBS1 + IBS2)), PPC (IBS
#' binomial test), RATIO (HETHET : IBS0 SNP ratio (expected value 2)).
#' and ii) failIDs containing a [data.frame] with individual IDs [IID] and
#' family IDs [FID] of individuals failing the highIBDTh; iii) p_IBD, a
#' ggplot2-object 'containing' all pair-wise IBD-estimates as histograms
#' stratified by value of PI_HAT, which can be shown by print(p_IBD and iv)
#' plot_data, a data.frame with the data visualised in p_IBD (iii).
#' @export
#' @examples
#' qcdir <- system.file("extdata", package="plinkQC")
#' name <- 'data'
#' \dontrun{
#' relatednessQC <- evaluate_check_relatedness(qcdir=qcdir, name=name,
#' interactive=FALSE)
#' }
evaluate_check_relatedness <- function(qcdir, name, highIBDTh=0.1875,
                                       imissTh=0.03, interactive=FALSE,
                                       legend_text_size = 5,
                                       legend_title_size = 7,
                                       axis_text_size = 5,
                                       axis_title_size = 7,
                                       title_size = 9,
                                       verbose=FALSE) {

    prefix <- makepath(qcdir, name)

    if (!file.exists(paste(prefix, ".imiss", sep=""))){
        stop("plink --missing output file: ", prefix,
             ".imiss does not exist.")
    }
    if (!file.exists(paste(prefix, ".genome",sep=""))){
        stop("plink --genome output file: ", prefix,
             ".genome does not exist.")
    }
    testNumerics(numbers=highIBDTh, positives=highIBDTh, proportions=highIBDTh)
    names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
    imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE, stringsAsFactors=FALSE)
    if (!all(names_imiss == names(imiss))) {
        stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
    names_genome <- c("FID1", "IID1", "FID2", "IID2", "RT", "EZ", "Z0", "Z1",
                      "Z2", "PI_HAT", "PHE", "DST", "PPC", "RATIO")
    genome <- read.table(paste(prefix, ".genome", sep=""), header=TRUE,
                         as.is=TRUE, stringsAsFactors=FALSE)
    if (!all(names_genome == names(genome))) {
        stop("Header of ", prefix, ".genome is not correct. Was your
             file generated with plink --genome?")
    }

    fail_highIBD <- relatednessFilter(relatedness=genome,  otherCriterion=imiss,
                                      relatednessTh=highIBDTh,
                                      relatednessFID1="FID1",
                                      relatednessFID2="FID2",
                                      otherCriterionTh=imissTh,
                                      otherCriterionThDirection="gt",
                                      otherCriterionMeasure="F_MISS" )

    genome$PI_HAT_bin <- ifelse(genome$PI_HAT > 0.05, 0, 1)
    p_allPI_HAT <- ggplot(genome, aes_string('PI_HAT'))
    p_allPI_HAT <- p_allPI_HAT + geom_histogram(binwidth = 0.005,
                                                fill="#66a61e") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for all sample pairs") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_text_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_highPI_HAT <- ggplot(dplyr::filter(genome, .data$PI_HAT_bin == 0),
                           aes_string('PI_HAT'))
    p_highPI_HAT <- p_highPI_HAT + geom_histogram(binwidth = 0.005,
                                                  fill="#e6ab02") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for sample pairs with PI_HAT >0.1") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_text_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_histo <- cowplot::plot_grid(p_allPI_HAT, p_highPI_HAT)
    title <- cowplot::ggdraw() +
        cowplot::draw_label("Relatedness estimated as pairwise IBD (PI_HAT)",
                            size=title_size)
    p_IBD <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                rel_heights = c(0.1, 1))
    if (interactive) print(p_IBD)
    return(list(fail_highIBD=fail_highIBD$relatednessFails,
                failIDs=fail_highIBD$failIDs, p_IBD=p_IBD,
                plot_data = genome))
}



