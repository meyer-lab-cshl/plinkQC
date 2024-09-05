context('Test ancestry functions')
qcdir <- '.'
indir <- '.'
name <- 'data.hg38.renamed'

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
  expect_no_error(superpop_classification(indir = indir, qcdir = qcdir, name = name)
  )
})

test_that("Superpop_classification works", {
  mockery::stub(where = superpop_classification,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = superpop_classification,
                what = "system2",
                how = TRUE)
  expect_no_error(superpop_classification(indir = indir, qcdir = qcdir, name = name)
  )
})

test_that("Throws file error if data is in wrong format for convert to plink", {
  mockery::stub(where = superpop_classification,
                what = "checkPlink2",
                how = TRUE)
  mockery::stub(where = superpop_classification,
                what = "system2",
                how = TRUE)
  expect_no_error(convert_to_plink2(indir = indir, qcdir = qcdir, 
                                          name = "data.hg38.renamed")
  )
})

