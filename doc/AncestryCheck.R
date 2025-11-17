## ----setup knitr, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setting up, eval = FALSE-------------------------------------------------
# package.dir <- find.package('plinkQC')
# indir <- file.path(package.dir, 'extdata')
# qcdir <- tempdir()
# name <- 'data.hg38'
# path2plink2 <- "/Users/syed/bin/plink2"
# path2load_mat <- "/Users/syed/Documents/loading_matrix"

## ----eval = FALSE-------------------------------------------------------------
# convert_to_plink2(indir=indir, qcdir=qcdir, name=name,
#                   path2plink2 = path2plink2)
# rename_variant_identifiers(indir=qcdir, qcdir=qcdir, name=name,
#                            path2plink2 = path2plink2)
# name <- paste0(name, ".renamed")

## ----ancestry identification, eval = FALSE------------------------------------
# ancestries <- superpop_classification(indir=qcdir, qcdir=qcdir, name=name,
#                                       path2plink2 = path2plink2,
#                                       path2load_mat = path2load_mat)
# ancestries$p_ancestry

## ----load ancestry, out.width = "500px", echo=FALSE---------------------------
par(mfrow=c(2,1), las=1)
knitr::include_graphics("ancestral_prediction.pdf")

## ----renaming, eval = FALSE---------------------------------------------------
# library(plinkQC)
# name <- "data.no_ac_gt_snps"
# refname <- "all_hg38.no_ac_gt_snps"
# path2plink2 <- "/Users/syed/bin/plink2"
# rename_variant_identifiers(indir=qcdir, qcdir=qcdir, name=name,
#                            path2plink2 = path2plink2)

## ----marker_sample_QC, eval=FALSE---------------------------------------------
# refname <- "all_hg38.renamed.studysnps"
# 
# fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name=name,
#                             path2plink=path2plink,
#                             verbose=TRUE, interactive=TRUE,
#                             showPlinkOutput=FALSE)
# marker_ids <- cleanData(indir = indir, qcdir = qcdir, name = name,
# 				path2plink = path2plink, filterAncestry = FALSE, filterRelated = FALSE, macTh = NULL, mafTh = 0.05, verbose = TRUE,
# 				filterSex = FALSE, filterHeterozygosity = FALSE,
# 				filterSNPMissingness = TRUE, filterSampleMissingness = FALSE,
# 				filterMAF = TRUE, filterHWE = TRUE)
# name = paste(name, ".clean", sep = "")
# 
# pruning_ld(indir = indir, qcdir = qcdir, name = name,
# 	  path2plink = path2plink, genomebuild="hg38")
# name = paste(name, ".pruned", sep = "")
# 
# fail_samples <- perIndividualQC(indir=indir, qcdir=qcdir, name=name, dont.check_sex = TRUE,
# 								dont.check_relatedness = TRUE,
#                             path2plink=path2plink, dont.check_ancestry = TRUE,
#                             interactive=TRUE, verbose=TRUE)
# 
# sample_ids <- cleanData(indir = indir, qcdir = qcdir, name = name,
# 			 path2plink = path2plink, filterAncestry = FALSE, filterRelated = FALSE,verbose = TRUE,
# 			filterSex = FALSE, filterHeterozygosity = TRUE,
# 			filterSNPMissingness = FALSE, filterSampleMissingness = TRUE,
# 			filterMAF = FALSE, filterHWE = FALSE)
# 
# 

## ----eval = FALSE-------------------------------------------------------------
# library(tidyverse)
# library(randomForest)
# 
# filepath <- 'insert path to sscore file here'
# proj <- read.csv(
#   file= filepath,
#   sep='\t', header = TRUE)
# 
# package.dir <- find.package('plinkQC')
# ancestry_info <-
#   read_delim(file.path(package.dir,"extdata/Genomes1000_ID2Pop.txt"))
# superpop <-
#   read_delim(file.path(package.dir,"extdata/AncestryGeoLocations.csv"))

## ----eval = FALSE-------------------------------------------------------------
# proj <- proj %>%
#   select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
# colnames(proj) <- c("IID", paste0("PC", 1:20))
# 
# labeled_proj <- merge(proj, superpop)
# 
# n_individuals <- 448
# set.seed(123)
# idx <- split(seq(nrow(labeled_proj)), labeled_proj$Ancestry)
# 
# train.ids <- lapply(idx, function(i) {
#   sample(i, n_individuals)
# })
# 
# ids_to_keeps<- unlist(train.ids, use.names = FALSE, recursive = FALSE)
# train_proj_1000g <- labeled_proj[ids_to_keeps,]
# 

## ----eval = FALSE-------------------------------------------------------------
# train_proj_1000g$Ancestry <- factor(train_proj_1000g$Ancestry)
# 
# ancestry_rf <- randomForest(Ancestry ~ .,
#                         data = train_proj_1000g[,-c(2,23)],
#                         method = "rf",
#                         ntree = 750,
#                         importance = TRUE)
# ancestry_rf
# 

## ----eval = FALSE-------------------------------------------------------------
# filepath <- "Insert path to newdata.sscore here"
# newdata <- read.csv(
#   file= filepath,
#   sep='\t', header = TRUE)
# 
# newdata <- newdata %>%
#   select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
# colnames(newdata) <- c("ID", paste0("PC", 1:20))
# 
# predictions <- predict(ancestry_rf, newdata)

## ----eval = FALSE-------------------------------------------------------------
# ancestry_rf

## ----eval = FALSE-------------------------------------------------------------
# my_rf <- ancestry_rf$confusion
# heatmap(my_rf[,-c(27)], scale = 'column', trace = 'none')

## ----eval = FALSE-------------------------------------------------------------
# train_results <- data.frame(mtry = c(0), ntrees = c(0), num_pc = c(0), acc = c(0))
# train_proj_noids <- train_proj_1000g[,-c(1)]
# for (PC_inc in c(1:20)) {
#   for (ntree in c(1, 5)) {
#     set.seed(123)
#     fit <- train(Pop~., data=train_proj_noids[,c(1:PC_inc,21)],
#                  method="rf", metric="Accuracy",
#                  trControl=control, ntree=ntree)
#     train_proj_1000g <- rbind(train_proj_1000g,
#                      data.frame(mtry = fit$results$mtry,
#                                 ntrees = rep(ntree, nrow(fit$results)),
#                                 num_pc = rep(PC_inc, nrow(fit$results)),
#                                 acc = fit$results$Accuracy))
#   }
# }
# 
# 
# 

## ----eval = FALSE-------------------------------------------------------------
# # This is a visual of the random forest being made.
# # https://github.com/araastat/reprtree
# reprtree:::plot.getTree(rf_mtry)
# # Two different plots
# reprtree:::plot.getTree(rf_mtry, k=1, d=5)
# reprtree:::plot.getTree(rf_mtry, k=1, d=6)

