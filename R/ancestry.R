#' Identifying ancestries of samples from plink data samples
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
#' @return 
#' 
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
  
  if (verbose) message("Projecting data on 1000G PC space via Plink2 --score")
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
#' categories: AFR, AMR, EAS, EUR, and SAS. 
#' 
#' @param indir [character] /path/to/directory containing the basic PLINK 2.0 data
#' file name.pgen, name.pvar, name.psam
#' @param qcdir [character] /path/to/directory where name.sscore as returned by 
#' plink2 --score will be saved to. User needs writing permission to qcdir. Per 
#' default is qcdir=indir.  
#' @param name [character] Prefix of PLINK 2.0 files, i.e. name.pgen, name.pvar, 
#' name.psam
#' @inheritParams checkPlink2
#' @param showPlinkOutput [logical] If TRUE, plink log and error messages are
#' printed to standard out.
#' @param verbose [logical] If TRUE, progress info is printed to standard out.
#' @return A named [factor] with five different levels AFR, AMR, EAS, EUR, and SAS that
#' correspond to the ancestry that a sample comes from and the names corresponds
#' to the IDs.
#' @export 
superpop_classification <- function(indir, name, qcdir=indir, verbose=FALSE,
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

  ref_pca_acount <- system.file("extdata", 'all_hg38_maf_first_pca.acount',
                                    package="plinkQC")
  ref_pca_eigen_allele <- system.file("extdata", 'all_hg38_maf_first_pca.eigenvec.allele',
                                package="plinkQC")
  
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
                 "--read-freq", ref_pca_acount,
                 "--score", ref_pca_eigen_allele,
                 "2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-26",
                 "--out", out),
          stdout = showPlinkOutput, stderr = showPlinkOutput)

  proj <- read.csv(paste0(out,".sscore"), 
                   sep='\t', header = TRUE)
  proj <- proj %>%
    select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
  colnames(proj) <- c("ID", paste0("PC", 1:20))
  
  #load RF
  load("R/sysdata.rda")
  predictions = predict(proj_superpop, proj)
  return(setNames(predictions, proj$ID))
  #delete files??
}


  
  
#putting it here as an example
run_check_heterozygosity <- function(indir, name, qcdir=indir, verbose=FALSE,
                                       path2plink=NULL,
                                       keep_individuals=NULL,
                                       remove_individuals=NULL,
                                       exclude_markers=NULL,
                                       extract_markers=NULL,
                                       showPlinkOutput=TRUE) {
    
    prefix <- makepath(indir, name)
    out <- makepath(qcdir, name)
    
    checkFormat(prefix) #ok so making sure prefix is right format
    path2plink <- checkPlink(path2plink)
    args_filter <- checkFiltering(keep_individuals, remove_individuals,
                                  extract_markers, exclude_markers) #this is just if you want to filter individuals
    
    if (verbose) message("Run check_heterozygosity via plink --het")
    sys::exec_wait(path2plink,
                   args=c("--bfile", prefix, "--het", "--out", out,
                          args_filter),
                   std_out=showPlinkOutput, std_err=showPlinkOutput)
  }
  
  

  
  