# plinkQC 0.2.3
## major changes
* Enable return of overview plots as ggplot object
  (fixes #6 in ab2e840f0f22ccdeb5317475698aa58e7eecd345 and
  101e74e318752b2f8038330245d875f6afd5b57b)
* Relationship filter now deals with more complicated relationship scenarios as
  observed in plant genotype sets (fixed #11)

## minor changes
* give user option to choose maf threshold for relatedness filtering (relates to
  #3)

## bug fixes
* Include check for zero related individuals fixing #3 in
  1445a88d41d985e73ffdd161144229683f8352cd
* Include check in case all samples fail perIndividual QC in
   894acc1fa03dadfe0ad2028888142171bcc641eb and
   04642246d18ed4eaac5b9d9a6931d1ecb08308e8)
* Include checks for diagonal derived relationship estimates, and estimate data
  containing only related individuals; fixes #11
* Fix command for genotype conversion in 1000Genomes vignette, addressing issue
  #10
* fix missing rownames error for overviewPerIndividualQC, when relatedness check
  was included (issue #16, fc7a38b1f2b345d9c6c5d69f5dcf0bc57a857f62)


# plinkQC 0.2.2
## minor changes
* Fix IDs written to fail.IDs file: previous versions wrote IID,IID, now fixed
to FID, IID (fixes #2).

# plinkQC 0.2.1
## minor changes
* Fix path check bug in checkPlink
* Include test data in build!

# plinkQC 0.2.0
## major changes
* All system calls to plink are conducted with sys::exec_wait - this should
solve platform dependent issues, mainly comptatibility with windows.
* Make path construction compatible with windows.
* path2plink now requires full path to plink executable, no tilde expansion or
simple pointer to directory supported.
* Fix bug in return of cleanData function: list now contains keep and fail IDs.
* Fix bug in return of maf computation: if fail.IDs does not exist, set
  fail_samples to zero.

## minor changes
* IBD-fail.IDs now saved without column names to be consistent with other
xxx-fail.IDs files.
* Include additional progress messages in cleanData()
* Remove default double-specification of mafTh and macTh
* use checkPlink to return correct path2plink and export to make checkPlink
directly accesible to user.


# plinkQC 0.1.1
## major changes
* run_check_relatedness will only save IBD estimates of individuals whose
estimates are higher than the threshold.

## minor changes
* Fix examples in vignettes 1000 Genomes and HapMap III reference.
* Change file access in function examples
* Add additional checks in check_ancestry and fix missing refSamplesFile test


# plinkQC 0.1.0

* Added a `NEWS.md` file to track changes to the package.



