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
  
  checkFormat(prefix)
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

#' Predicting sample superpopualation ancestry 
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
#' @param axis_text_size [integer] Size for axis text.
#' @param axis_title_size [integer] Size for axis title.
#' @param title_size [integer] Size for plot title.
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return A named [factor] with five different levels AFR, AMR, EAS, EUR, and SAS that
#' correspond to the ancestry that a sample comes from and the names corresponds
#' to the IDs.
#' @export 
#' @examples
#' indir <- system.file("extdata", package="plinkQC")
#' qcdir <- tempdir()
#' name <- "data.hg38"
#' path2plink <- '/path/to/plink'
#' # the following code is not run on package build, as the path2plink on the
#' # user system is not known.
#' \dontrun{
#' superpop_classification(indir=indir, qcdir=qcdir, name=name, path2plink2 = path2plink2)
#' }
superpop_classification <- function(indir, name, qcdir=indir, verbose=FALSE,
                                    path2plink2=NULL,
                                    path2load_mat=NULL,
                                    keep_individuals=NULL,
                                    remove_individuals=NULL,
                                    exclude_markers=NULL,
                                    extract_markers=NULL,
                                    axis_text_size = 5,
                                    axis_title_size = 7,
                                    title_size = 9,
                                    showPlinkOutput=TRUE) {
  
  prefix <- makepath(indir, name)
  out <- makepath(qcdir, name) 

  checkFormatPlink2(prefix)
  path2plink2 <- checkPlink2(path2plink2)
  checkLoadingMat(path2load_mat)

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
                 "--read-freq", makepath(path2load_mat, "all_hg38.pca.acount"),
                 "--score", makepath(path2load_mat, "all_hg38.pca.eigenvec.allele"),
                 "2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-26",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)

  proj <- read.csv(paste0(out,".sscore"), 
                   sep='\t', header = TRUE)

  colnames(proj) <- c("IID", "Allele_Count", "Allele_Dosage", paste0("PC", 1:20))
  
  #load RF
  rf_path <- system.file("extdata", 'superpop_rf_0909.RDS',
              package="plinkQC")
  superpop <- readRDS(rf_path)
  predictions <- predict(superpop, proj)
  names(predictions) <- proj$IID
  
  predictions <- data.frame(predictions)
  
  p_ancestry <- 
    ggplot(predictions, aes(x = predictions, fill = predictions)) + geom_bar() + 
    xlab("Ancestral Prediction") + ylab("Count") + 
    theme_bw() + 
    theme(legend.position = "none",
          title = element_text(size = title_size),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size))
  return(list(predictions=predictions, p_ancestry = p_ancestry))
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
                 "--set-all-var-ids", paste0('"', format, '"'),
                 "--rm-dup exclude-all",
                 "--make-pgen",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)
}



checkLoadingMat <- function(path2load_mat) {
  if (!file.exists(makepath(path2load_mat, "all_hg38.pca.acount"))) {
    stop("The file all_hg38.pca.acount is not found in the path given. Please 
         check that the file path is correct")
  }
  if (!file.exists(makepath(path2load_mat, "all_hg38.pca.eigenvec.allele"))) {
    stop("The file all_hg38.pca.eigenvec.allele is not found in the path given. 
    Please check that the file path is correct")
  }
}

