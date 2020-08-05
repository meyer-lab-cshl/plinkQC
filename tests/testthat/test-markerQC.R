context('Test markerQC functions')
qcdir <- '.'
indir <- '.'
name <- 'data'


context('Test check_snp_missingness')
test_that('check_snp_missingness throws file error',{
    expect_error(check_snp_missingness(qcdir=qcdir, indir=indir, name="nodata",
                                       verbose=FALSE),
                 "plink family file")
})
test_that('check_snp_missingness runs with no perIndividualQC fail IDs', {
    mockery::stub(where = check_snp_missingness,
                  what = "checkPlink",
                  how = TRUE)
    mockery::stub(where = check_snp_missingness,
                  what = "sys::exec_wait",
                  how = TRUE)
    expect_message(check_snp_missingness(qcdir=qcdir, indir=indir,
                                         name="data.clean",
                             verbose=FALSE),
                   "No individuals failed perIndividualQC")
})

context('Test check_hwe')
test_that('check_hwe throws file error',{
    expect_error(check_hwe(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})
test_that('check_hwe runs with no perIndividualQC fail IDs', {
    mockery::stub(where = check_hwe,
                  what = "checkPlink",
                  how = TRUE)
    mockery::stub(where = check_hwe,
                  what = "sys::exec_wait",
                  how = TRUE)
    expect_message(check_hwe(qcdir=qcdir, indir=indir, name="data.clean",
                       verbose=FALSE),
             "No individuals failed perIndividualQC")
})

context('Test check_maf')
test_that('check_maf throws file error',{
    expect_error(check_maf(qcdir=qcdir, indir=indir, name="nodata",
                           verbose=FALSE),
                 "plink family file")
})
test_that('check_maf runs with no perIndividualQC fail IDs', {
    mockery::stub(where = check_maf,
                  what = "checkPlink",
                  how = TRUE)
    mockery::stub(where = check_maf,
                  what = "sys::exec_wait",
                  how = TRUE)
    expect_message(check_maf(qcdir=qcdir, indir=indir, name="data.clean",
                             verbose=FALSE),
                   "No individuals failed perIndividualQC")
})

#fail_markers <- perMarkerQC(indir=indir, qcdir=qcdir, name="data.clean", path2plink = path2plink)
