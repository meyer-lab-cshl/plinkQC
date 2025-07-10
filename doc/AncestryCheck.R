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
# head(ancestries$predictions)

## ----eval = FALSE-------------------------------------------------------------
# library(tidyverse)
# library(randomForest)
# 
# filepath <- '~/qcdir/filtered_hg38.projection.sscore'
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
# labeled_proj <- merge(proj, superpop_id)
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
# rf_100_trees <- randomForest(Ancestry ~ .,
#                    data = train_proj_1000g[,-c(1,21,22)],
#                    ntree = 100,
#                    method = "rf",
#                    importance = TRUE)
# 
# rf_500_trees <- randomForest(Ancestry ~ .,
#                    data = train_proj_1000g[,-c(1,21,22)],
#                    ntree = 500,
#                    method = "rf",
#                    importance = TRUE)
# 
# rf_10000_trees <- randomForest(Ancestry ~ .,
#                    data = train_proj_1000g[,-c(1,21,22)],
#                    ntree = 10000,
#                    method = "rf",
#                    mtry = best_mtry,
#                    importance = TRUE)

## ----eval = FALSE-------------------------------------------------------------
# # This first line will provide the best mtry
# mtry<- tuneRF(x = train_proj_1000g[,-c(1,21,22)], y = train_proj_1000g$Ancestry,
#               stepFactor = 1.5, improve = 0.01, trace = TRUE, plot = TRUE)
# best_mtry <- mtry[mtry[, 2] == min(mtry[, 2]), 1]
# 
# set.seed(123)
# rf_mtry <- randomForest(Ancestry ~ .,
#                         data = train_proj_1000g[,-c(1,21,22)],
#                         mtry = best_mtry,
#                         ntree = 1500,
#                         method = "rf",
#                         importance = TRUE)
# rf_mtry

## ----eval = FALSE-------------------------------------------------------------
# library(caret)
# 
# set.seed(123)
# trainControl_5 <- trainControl(method = "cv",
#                           number = 5,
#                           search = "grid")
# 
# forest_5 <- train(Ancestry ~ .,
#                data=train_proj_1000g[,-c(1,21,22)],
#                method = "rf",
#                metric = "Accuracy",
#                trControl = trainControl_5
#                )
# # Print out the model details
# forest_5$finalModel
# 
# # Now for 10-fold cross validation
# trainControl_10 <- trainControl(method = "cv",
#                           number = 10,
#                           search = "grid")
# forest_10 <- train(Ancestry ~ .,
#                data=train_proj_1000g[,-c(1,21,22)],
#                method = "rf",
#                metric = "Accuracy",
#                trControl = trainControl_10
#                )
# 
# # Print out model details
# rf_10_fold_cv <- forest_10$finalModel
# 

## ----eval = FALSE-------------------------------------------------------------
# # This is a visual of the random forest being made.
# # https://github.com/araastat/reprtree
# reprtree:::plot.getTree(rf_mtry)
# # Two different plots
# reprtree:::plot.getTree(rf_mtry, k=1, d=5)
# reprtree:::plot.getTree(rf_mtry, k=1, d=6)

