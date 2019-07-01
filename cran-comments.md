## Test environments
* local OS X install, R 3.6.0
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Comments
Minor change to version 0.2.1:
* Fix IDs written to fail.IDs file: previous versions wrote IID,IID, now fixed
to FID, IID (see [here](https://github.com/HannahVMeyer/plinkQC/issues/2)).
