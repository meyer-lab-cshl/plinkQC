context('Test ancestry functions')
qcdir <- '.'
indir <- '.'
name <- 'data.hg38.renamed'
path2load_mat <- "../plinkQC_validation/ancestry_unrel_model/merged_chrs.postQC.train.pca"


context('Test run_ancestry_format')
test_that('run_ancestry_format throws file error',{
  mockery::stub(where = run_ancestry_format,
                what = "checkPlink2",
                how = TRUE)
  expect_error(run_ancestry_format(qcdir=qcdir, indir=indir, name="nodata",
                         verbose=FALSE),
               "does not exist")
})

test_that("run_ancestry_prediction works", {
  mockery::stub(where = run_ancestry_prediction,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = run_ancestry_prediction,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = run_ancestry_prediction,
                what = "system2",
                how = TRUE)
  mockery::stub(where = run_ancestry_prediction,
                what = "checkLoadingMat",
                how = TRUE)
  expect_no_error(run_ancestry_prediction(indir = indir, qcdir = qcdir, 
                                          name = name, path2load_mat = path2load_mat)
  )
})

context('Testing convert to plink2')
test_that("Throws file error if data is in wrong format for convert to plink", {
  mockery::stub(where = convert_to_plink2,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = convert_to_plink2,
                what = "system2",
                how = TRUE)
  expect_error(convert_to_plink2(indir = indir, qcdir = qcdir, name="nodata",
                                 verbose=FALSE),
               "does not exist")

})

test_that("Throws error if pathtoplink isn't right", {
  mockery::stub(where = convert_to_plink2,
                what = "system2",
                how = TRUE)
  expect_error(convert_to_plink2(indir = indir, qcdir = qcdir, name=name,
                                 verbose=FALSE),
               "PLINK 2.0 software required")
  
})

test_that("Testing that the loading matrix path will not run if the path is not correct", {
  expect_error(checkLoadingMat(path2load_mat = "Nothing"),
               "The loading matrix .acount file is not found in the path \n         given for path2load_mat. Please check that the file path is correct.\n         Note that that filepath requires for the filename before .acount to be \n         included.")
  
})



