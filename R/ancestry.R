#' Running functions to format data for ancestry identification
#' 
#' This function runs convert_to_plink2 and rename_variant_identifiers to format 
#' the data for the ancestry identification with superpop_classification
#' 
#' @inheritParams convert_to_plink2
#' @inheritParams rename_variant_identifiers
#' @inheritParams checkPlink2
#' @inheritParams superpop_classification
#' @param plink2format [logical] If TRUE, data is in plink2 format already and 
#' convert_to_plink2 will not be run
#' @param var_format [logical] If TRUE, variant identifiers are in correct 
#' format already and rename_variant_identifiers will not be run
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Creates plink 2.0 datafiles
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' run_ancestry_format(indir=indir, qcdir=qcdir, 
#'   name=name, path2plink2 = path2plink2)
#' }
run_ancestry_format <- function(indir, name, qcdir=indir, verbose=FALSE,
                                path2plink2=NULL,
                                keep_individuals=NULL,
                                remove_individuals=NULL,
                                exclude_markers=NULL,
                                extract_markers=NULL,
                                showPlinkOutput=TRUE,
                                format = "@:#[hg38]",
                                plink2format = FALSE,
                                var_format = FALSE,
                                path2load_mat, legend_text_size,
                                legend_title_size, axis_text_size,
                                axis_title_size, title_size,
                                legend_position = "bottom") {
  
  path2plink2 <- checkPlink2(path2plink2)
  if (plink2format==FALSE) {
    convert_to_plink2(indir=indir, name=name, qcdir=qcdir, verbose=verbose,
                      path2plink2=path2plink2,
                      keep_individuals=keep_individuals,
                      remove_individuals=remove_individuals,
                      exclude_markers=exclude_markers,
                      extract_markers=extract_markers,
                      showPlinkOutput=showPlinkOutput)
    indir=qcdir
  }
  
  if (var_format==FALSE) {
    rename_variant_identifiers(indir=indir, name=name, qcdir=qcdir, verbose=verbose,
                      path2plink2=path2plink2,
                      format=format,
                      showPlinkOutput=showPlinkOutput)
    name = paste0(name, ".renamed")
  }
  
  results <- superpop_classification(indir=indir, name=name, qcdir=qcdir, 
                                     verbose=verbose, interactive=interactive,
                          path2plink2=path2plink2,
                          path2load_mat=path2load_mat,
                          keep_individuals=keep_individuals,
                          remove_individuals=remove_individuals,
                          extract_markers=extract_markers,
                          exclude_markers=exclude_markers,
                          legend_text_size=legend_text_size,
                          legend_title_size=legend_title_size,
                          axis_text_size=axis_text_size,
                          axis_title_size=axis_title_size,
                          title_size=title_size,
                          showPlinkOutput=showPlinkOutput, 
                          legend_position=legend_position)
  return(results)
}



#' Converting PLINK v1.9 data files into PLINK v2.0 data files
#' 
#' This converts files in the PLINK v1.9 format (i.e. name.bim, name.fam, and 
#' name.bed) into PLINK v2.0 format (i.e. name.pvar, name.psam, and name.pgen)
#' 
#' @param indir [character] /path/to/directory containing the basic PLINK 1.9 data
#' file name.bim, name.fam, name.bed
#' @param qcdir [character] /path/to/directory where the plink2 data formations
#' as returned by plink2 --make-pgen will be saved to. User needs writing 
#' permission to qcdir. Per default is qcdir=indir.  
#' @param name [character] Prefix of PLINK 1.9 files, i.e. name.bim, name.fam, 
#' name.bed
#' @inheritParams checkPlink2
#' @inheritParams checkFiltering
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Creates plink 2.0 datafiles
#' @export
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' convert_to_plink2(indir=indir, qcdir=qcdir, name=name, path2plink2 = path2plink2)
#' }
convert_to_plink2 <- function(indir, name, qcdir=indir, verbose=FALSE,
                              path2plink2=NULL,
                              keep_individuals=NULL,
                              remove_individuals=NULL,
                              exclude_markers=NULL,
                              extract_markers=NULL,
                              showPlinkOutput=TRUE) {
  
  prefix <- makepath(indir, name)
  out <- makepath(qcdir, name) 
  
  checkFormatPlink2(prefix)
  path2plink2 <- checkPlink2(path2plink2)
  
  if (showPlinkOutput) {
    showPlinkOutput = ""}
  else {
    showPlinkOutput = FALSE}
  
  args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                extract_markers, exclude_markers) 
  
  if (verbose) message("Converting to plink2 data types")
  system2(path2plink2, 
          args=c("--bfile", prefix,
                 args_filter,
                 "--make-pgen",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)
}

#' Predicting sample superpopulation ancestry 
#' 
#' Predicts the ancestry of inputted samples using plink2. This program utilizes a random 
#' forest algorithm that was trained on the 1000 genomes dataset to predict the 
#' ancestry of the samples. There are five possible ancestral superpopulation
#' categories: AFR, AMR, EAS, EUR, and SAS. Genomic data version hg38 with variant 
#' identifiers in the format of 1:12345[hg38] is needed for the function to work 
#' 
#' @param indir [character] /path/to/directory containing the basic PLINK 2.0 data
#' file name.pgen, name.pvar, name.psam
#' @param qcdir [character] /path/to/directory where name.sscore as returned by 
#' plink2 --score will be saved to. User needs writing permission to qcdir. Per 
#' default is qcdir=indir.  
#' @param name [character] Prefix of PLINK 2.0 files, i.e. name.pgen, name.pvar, 
#' name.psam
#' @param path2load_mat [character] /path/to/directory where loading matrices are 
#' kept. This can be downloaded from: https://github.com/meyer-lab-cshl/plinkQCAncestryData
#' @inheritParams checkPlink2
#' @inheritParams checkFiltering
#' @param legend_text_size [integer] Size for legend text.
#' @param legend_title_size [integer] Size for legend title.
#' @param axis_text_size [integer] Size for axis text.
#' @param axis_title_size [integer] Size for axis title.
#' @param title_size [integer] Size for plot title.
#' @param legend_position [character] Legend position for the plot. 
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param interactive [logical] Should plots be shown interactively? When
#' choosing this option, make sure you have X-forwarding/graphical interface
#' available for interactive plotting. Alternatively, set interactive=FALSE and
#' save the returned plot object (p_ancestry) via ggplot2::ggsave(p=p_ancestry,
#' other_arguments) or pdf(outfile) print(p_ancestry) dev.off().
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return A named [factor] with five different levels AFR, AMR, EAS, EUR, and 
#' SAS that correspond to the ancestry that a sample comes from and the names 
#' corresponds to the IDs.
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data.hg38"
#' path2plink <- '/path/to/plink'
#' path2load_mat <- '/path/to/load_mat/merged_chrs.postQC.train.pca'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' superpop_classification(indir=indir, qcdir=qcdir, name=name, 
#' path2plink2 = path2plink2, path2load_mat = path2load_mat)
#' }
#'@export 
superpop_classification <- function(indir, name, qcdir=indir, verbose=FALSE,
                                    interactive = FALSE,
                                    path2plink2=NULL,
                                    path2load_mat=NULL,
                                    keep_individuals=NULL,
                                    remove_individuals=NULL,
                                    extract_markers=NULL,
                                    exclude_markers=NULL,
                                    legend_text_size = 5,
                                    legend_title_size = 7,
                                    axis_text_size = 5,
                                    axis_title_size = 7,
                                    title_size = 9,
                                    showPlinkOutput=TRUE,
                                    legend_position="right") {
  
  prefix <- makepath(indir, name)
  out <- makepath(qcdir, name) 

  checkFormatPlink2(prefix)
  path2plink2 <- checkPlink2(path2plink2)
  #checkLoadingMat(path2load_mat)

  if (showPlinkOutput) {
    showPlinkOutput = ""
  }
  else {
    showPlinkOutput = FALSE 
  }
  
  args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                extract_markers, exclude_markers) 
  
  if (verbose) message("Projecting data on 1000G PC space via Plink2 --score")
  system2(path2plink2,
          args=c("--pfile", prefix,
                 args_filter,
                 "--snps-only",
                 "--max-alleles 2",
                 "--read-freq", paste0(path2load_mat, ".acount"),
                 "--score", paste0(path2load_mat, ".eigenvec.allele"),
                 "2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-26",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)

  proj <- read.csv(paste0(out,".sscore"), 
                   sep='\t', header = TRUE)
  
  #seeing if FID is part of the projection
  if ("X.FID" %in% names(proj)) {
    colnames(proj) <- c("FID", "IID", "Allele_Count", "Allele_Dosage", paste0("PC", 1:20))
  }
  else {
    colnames(proj) <- c("IID", "Allele_Count", "Allele_Dosage", paste0("PC", 1:20))
  }


  rf_path <- system.file("extdata", 'final_model.RDS',
                                    package="plinkQC")
  superpop <- readRDS(rf_path)
  prediction_prob <- predict(superpop, proj, type = "prob")
  prediction_prob <- data.frame(prediction_prob)
  prediction_prob <- data.frame(IID = proj["IID"], predictions = prediction_prob)
  
  prediction_majority <- predict(superpop, proj)
  prediction_majority <- data.frame(prediction_majority)
  prediction_majority <- data.frame(IID = proj["IID"], predictions = prediction_majority)
  
  colnames(prediction_prob) <- sub("^predictions.", "", colnames(prediction_prob))
  prediction_prob_long <- pivot_longer(prediction_prob, 
                                       cols = 'Africa':'Middle_East',
                                       names_to = "Ancestry",
                                       values_to = "predictions")
  p_ancestry <- 
     ggplot(prediction_prob_long, aes(x = as.factor(.data$IID), 
                                      y = .data$predictions, 
                                      fill = .data$Ancestry)) + 
     geom_bar(stat = "identity") +
     xlab("Samples") + ylab("Ancestral Prediction") + 
     theme_bw() + 
     theme(legend.position = legend_position,
           legend.text = element_text(size = legend_text_size),
           legend.title = element_text(size = legend_title_size),
           title = element_text(size = title_size),
           axis.text = element_text(size = axis_text_size),
           axis.title = element_text(size = axis_title_size)) +
    scale_fill_brewer(palette = "Set2") 
  
   return(list(prediction_prob=prediction_prob, 
               prediction_majority = prediction_majority, 
               p_ancestry = p_ancestry))
}


#' Renaming variants  
#' 
#' Changes the format of the variant identifier. The default is in the format of
#' chr1:12345[hg38].
#' 
#' @param indir [character] /path/to/directory containing the basic PLINK 2.0 data
#' file name.pgen, name.pvar, name.psam
#' @param qcdir [character] /path/to/directory where name.sscore as returned by 
#' plink2 --score will be saved to. User needs writing permission to qcdir. Per 
#' default is qcdir=indir.  
#' @param name [character] Prefix of PLINK 2.0 files, i.e. name.pgen, name.pvar, 
#' name.psam
#' @param format [character] This gives the template to rewrite the variant identifier.
#' A '@' represents the chromosome code, and a '#' represents the base-pair position.
#' @inheritParams checkPlink2
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return Files with a .renamed in them that have the renamed variants
#' @export 
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data.hg38"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' rename_variant_identifiers(indir=indir, qcdir=qcdir, name=name, path2plink2 = path2plink2)
#' }
rename_variant_identifiers <- function(indir, name, qcdir=indir, verbose=FALSE,
                                    path2plink2=NULL,
                                    format = "@:#[hg38]",
                                    showPlinkOutput=TRUE) {
  
  prefix <- makepath(indir, name)
  out <- makepath(qcdir, paste0(name, ".renamed")) 
  
  checkFormatPlink2(prefix)
  checkPlink2(path2plink2)
  
  if (showPlinkOutput) {
    showPlinkOutput = ""
  } 
  else {
    showPlinkOutput = FALSE 
  }
   
  system2(path2plink2, 
          args=c("--pfile", prefix,
                 "--snps-only",
                 "--max-alleles 2",
                 "--chr 1-22",
                 "--set-all-var-ids", paste0('"', format, '"'),
                 "--rm-dup exclude-all",
                 "--make-pgen",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)
}



checkLoadingMat <- function(path2load_mat) {
  if (!file.exists(makepath(path2load_mat, "unrel_hg38.maf.pruned.pca.acount"))) {
    stop("The file unrel_hg38.maf.pruned.pca.acount is not found in the path 
         given for path2load_mat. Please check that the file path is correct")
  }
  if (!file.exists(makepath(path2load_mat, "unrel_hg38.maf.pruned.pca.eigenvec.allele"))) {
    stop("The file unrel_hg38.maf.pruned.pca.eigenvec.allele is not found in the 
         path given for path2load_mat. Please check that the file path is correct")
  }
}

