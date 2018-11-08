# plinkQC 0.2.0
## major changes
* all system calls to plink are conducted with sys::exec_wait - this should
solve platform dependent issues, mainly comptatibility with windows.
* path2plink now requires full path to plink executable, no tilde expansion or
simple pointer to directory supported.
* fix bug in return of cleanData function: list now contains keep and fail IDs.
* fix bug in return of maf computation: if fail.IDs does not exist, set
  fail_samples to zero.

## minor changes
* IBD-fail.IDs now saved without column names to  be consistent with other
*-fail.IDs files.

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



