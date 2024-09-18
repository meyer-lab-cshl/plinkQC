context('Test ancestry functions')
qcdir <- '.'
indir <- '.'
name <- 'data.hg38.renamed'
path2load_mat <- '../../../loading_matrix'


context('Test superpop_classifcation')
test_that('Superpop_classification throws file error',{
  expect_error(superpop_classification(qcdir=qcdir, indir=indir, name="nodata",
                         verbose=FALSE),
               "does not exist")
})

test_that("Superpop_classification works", {
  mockery::stub(where = superpop_classification,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = superpop_classification,
                what = "system2",
                how = TRUE)
  mockery::stub(where = superpop_classification,
                what = "checkLoadingMat",
                how = TRUE)
  expect_no_error(superpop_classification(indir = indir, qcdir = qcdir, 
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


