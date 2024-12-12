## ----setup knitr, include = FALSE---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setting up, eval = FALSE-------------------------------------------------
#  package.dir <- find.package('plinkQC')
#  indir <- file.path(package.dir, 'extdata')
#  qcdir <- tempdir()
#  name <- 'data.hg38'
#  path2plink2 <- "/Users/syed/bin/plink2"
#  path2load_mat <- "/Users/syed/Documents/loading_matrix"

## ----eval = FALSE-------------------------------------------------------------
#  convert_to_plink2(indir=indir, qcdir=qcdir, name=name,
#                    path2plink2 = path2plink2)
#  rename_variant_identifiers(indir=qcdir, qcdir=qcdir, name=name,
#                             path2plink2 = path2plink2)
#  name <- paste0(name, ".renamed")

## ----ancestry identification, eval = FALSE------------------------------------
#  ancestries <- superpop_classification(indir=qcdir, qcdir=qcdir, name=name,
#                                        path2plink2 = path2plink2,
#                                        path2load_mat = path2load_mat)
#  head(ancestries$predictions)

## ----prune ref, eval = FALSE--------------------------------------------------
#  highld='high-LD-regions-hg38-GRCh38.txt'
#  
#  plink2 --pfile filtered_hg38.renamed\
#         --exclude range $qcdir/$highld \
#         --indep-pairwise 50 5 0.2 \
#         --out filtered_hg38.renamed
#  
#  plink2 --pfile filtered_hg38.renamed \
#         --extract filtered_hg38.renamed.prune.in \
#         --make-pgen \
#         --out filtered_hg38.pruned
#  
#  mv $qcdir/filtered_hg38.pruned $qcdir/plink_log

## ----eval = FALSE-------------------------------------------------------------
#  library(tidyverse)
#  library(randomForest)
#  
#  filepath <- '~/qcdir/filtered_hg38.projection.sscore'
#  proj <- read.csv(
#    file= filepath,
#    sep='\t', header = TRUE)
#  
#  package.dir <- find.package('plinkQC')
#  ancestry_info <-
#    read_delim(file.path(package.dir,"extdata/Genomes1000_ID2Pop.txt"))
#  superpop <-
#    read_delim(file.path(package.dir,"extdata/AncestryGeoLocations.csv"))

## ----eval = FALSE-------------------------------------------------------------
#  proj <- proj %>%
#    select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
#  colnames(proj) <- c("IID", paste0("PC", 1:20))
#  
#  labeled_proj <- merge(proj, superpop_id)
#  
#  n_individuals <- 448
#  set.seed(123)
#  idx <- split(seq(nrow(labeled_proj)), labeled_proj$Ancestry)
#  
#  train.ids <- lapply(idx, function(i) {
#    sample(i, n_individuals)
#  })
#  
#  ids_to_keeps<- unlist(train.ids, use.names = FALSE, recursive = FALSE)
#  train_proj_1000g <- labeled_proj[ids_to_keeps,]
#  

## ----eval = FALSE-------------------------------------------------------------
#  train_proj_1000g$Ancestry <- factor(train_proj_1000g$Ancestry)
#  
#  ancestry_rf <- randomForest(Ancestry ~ .,
#                          data = train_proj_1000g[,-c(2,23)],
#                          method = "rf",
#                          ntree = 750,
#                          importance = TRUE)
#  ancestry_rf
#  

## ----eval = FALSE-------------------------------------------------------------
#  filepath <- "Insert path to newdata.sscore here"
#  newdata <- read.csv(
#    file= filepath,
#    sep='\t', header = TRUE)
#  
#  newdata <- newdata %>%
#    select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM))
#  colnames(newdata) <- c("ID", paste0("PC", 1:20))
#  
#  predictions <- predict(ancestry_rf, newdata)

## ----eval = FALSE-------------------------------------------------------------
#  ancestry_rf

## ----eval = FALSE-------------------------------------------------------------
#  my_rf <- ancestry_rf$confusion
#  heatmap(my_rf[,-c(27)], scale = 'column', trace = 'none')

## ----eval = FALSE-------------------------------------------------------------
#  rf_100_trees <- randomForest(Ancestry ~ .,
#                     data = train_proj_1000g[,-c(1,21,22)],
#                     ntree = 100,
#                     method = "rf",
#                     importance = TRUE)
#  
#  rf_500_trees <- randomForest(Ancestry ~ .,
#                     data = train_proj_1000g[,-c(1,21,22)],
#                     ntree = 500,
#                     method = "rf",
#                     importance = TRUE)
#  
#  rf_10000_trees <- randomForest(Ancestry ~ .,
#                     data = train_proj_1000g[,-c(1,21,22)],
#                     ntree = 10000,
#                     method = "rf",
#                     mtry = best_mtry,
#                     importance = TRUE)

## ----eval = FALSE-------------------------------------------------------------
#  # This first line will provide the best mtry
#  mtry<- tuneRF(x = train_proj_1000g[,-c(1,21,22)], y = train_proj_1000g$Ancestry,
#                stepFactor = 1.5, improve = 0.01, trace = TRUE, plot = TRUE)
#  best_mtry <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
#  
#  set.seed(123)
#  rf_mtry <- randomForest(Ancestry ~ .,
#                          data = train_proj_1000g[,-c(1,21,22)],
#                          mtry = best_mtry,
#                          ntree = 1500,
#                          method = "rf",
#                          importance = TRUE)
#  rf_mtry

## ----eval = FALSE-------------------------------------------------------------
#  library(caret)
#  
#  set.seed(123)
#  trainControl_5 <- trainControl(method = "cv",
#                            number = 5,
#                            search = "grid")
#  
#  forest_5 <- train(Ancestry ~ .,
#                 data=train_proj_1000g[,-c(1,21,22)],
#                 method = "rf",
#                 metric = "Accuracy",
#                 trControl = trainControl_5
#                 )
#  # Print out the model details
#  forest_5$finalModel
#  
#  # Now for 10-fold cross validation
#  trainControl_10 <- trainControl(method = "cv",
#                            number = 10,
#                            search = "grid")
#  forest_10 <- train(Ancestry ~ .,
#                 data=train_proj_1000g[,-c(1,21,22)],
#                 method = "rf",
#                 metric = "Accuracy",
#                 trControl = trainControl_10
#                 )
#  
#  # Print out model details
#  rf_10_fold_cv <- forest_10$finalModel
#  

## ----eval = FALSE-------------------------------------------------------------
#  # This is a visual of the random forest being made.
#  # https://github.com/araastat/reprtree
#  reprtree:::plot.getTree(rf_mtry)
#  # Two different plots
#  reprtree:::plot.getTree(rf_mtry, k=1, d=5)
#  reprtree:::plot.getTree(rf_mtry, k=1, d=6)

## ----check ancestry, eval=FALSE, fig.height=3, fig.width=5, fig.align='center'----
#  library(plinkQC)
#  indir <- system.file("extdata", package="plinkQC")
#  qcdir <- "~/qcdir"
#  name <- 'data'
#  refname <- 'all_hg38'
#  prefixMergedDataset <- paste(name, ".", refname, sep="")
#  
#  exclude_ancestry <-
#      evaluate_check_ancestry(indir=qcdir, name=name,
#                              prefixMergedDataset=prefixMergedDataset,
#                              refSamplesFile=paste(indir,
#                                                   "/Genomes1000_ID2Pop.txt",
#                                                   sep=""),
#                              refColorsFile=paste(indir,
#                                                  "/Genomes1000_PopColors.txt",
#                                                  sep=""),
#                              interactive=TRUE)

## ----load ancestry, out.width = "500px", echo=FALSE, fig.align='center'-------
knitr::include_graphics("checkAncestry.png")

