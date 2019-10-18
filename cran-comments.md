## Resubmission
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


## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note


* This is a new release

## Comments
* move of code repository
* new urls
* major version update which should have been done with the previous release

