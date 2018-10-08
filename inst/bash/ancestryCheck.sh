###############################################################################
###                                                                         ###
### Sample ancestry check based on genotypes                                ###
###                                                                         ###
###     * numbers in QC indicate steps described in the protocol            ###
###      "Anderson et al. (2010) Data quality control in genetic            ###
###       case-control association studies." Nature protocols               ###
###       5(9):1564-73                                                      ###
###                                                                         ###
###                                                                         ###
###     by Hannah Meyer                                                     ###
###                                                                         ###
###############################################################################

####################################
### 1. input parameters and data ###
####################################

### a) parameters
alg=$1
qcdir=$2
highld=$3
hapmap=$4
hapmapdir=$5

mkdir -p $qcdir/plink_log

##########################
### 2. Ancestry check  ###
##########################

# 15. merge study genotypes to HapMap3 data
# Filter hapmap3 data for the same SNPset as in study
plink --bfile  $hapmapdir/$hapmap --extract $qcdir/$alg.prune.in \
      --make-bed --out $qcdir/$hapmap.pruned
mv  $qcdir/$hapmap.pruned.log $qcdir/plink_log/$hapmap.pruned

# check if SNPs which are in the HapMapIII (position matched to UK10K1KG)
# reference sets are concordant with HapMapIII info
# perfect match between SNPs in batch and HapMap
awk 'FNR==NR {{a[$1$2$4$5$6]; next}} $1$2$4$5$6 in a {print $2}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/$hapmap.match

# allele order mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/$hapmap.toUpdate

# possible allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$alg.pruned.bim $qcdir/$hapmap.pruned.bim > \
    $qcdir/$hapmap.toFlip

# all snps considered in first match to reference
cat $qcdir/$hapmap.match $qcdir/$hapmap.toFlip \
    $qcdir/$hapmap.toUpdate > $qcdir/$hapmap.keep

plink --bfile $qcdir/$hapmap.pruned \
      --update-alleles $qcdir/${hapmap}.toUpdate \
      --extract $qcdir/${hapmap}.keep \
      --flip $qcdir/${hapmap}.toFlip \
      --make-bed \
      --out $qcdir/${hapmap}.flipped
mv $qcdir/$hapmap.flipped.log \
    $qcdir/plink_log/$hapmap.flipped.log

plink --bfile $qcdir/${alg}.pruned \
      --extract $qcdi/${hapmap}.keep \
      --make-bed \
      --out $qcdir/${alg}.pruned.keep
mv $qcdir/$alg.pruned.keep.log \
    $qcdir/plink_log/$alg.pruned.keep.log

# allele order mismatch of flipped alleles
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] == $6$5)  {print $2,$5,$6,$6,$5}' \
    $qcdir/$alg.pruned.keep.bim $qcdir/$hapmap.flipped.bim >\
    $qcdir/${hapmap}.toUpdate_flipped

plink --bfile $qcdir/$hapmap.flipped \
      --update-alleles $qcdir/$hapmap.toUpdate_flipped \
      --make-bed \
      --out $qcdir/$hapmap.updated_flipped
mv $qcdir/$hapmap.updated_flipped.log \
    $qcdir/plink_log/$hapmap.updated_flipped.log

# snps that do not match between dataset and reference
awk 'BEGIN {OFS="\t"} FNR==NR {{a[$1$2$4]=$5$6; next}} \
    ($1$2$4 in a && a[$1$2$4] != $5$6)  {print $2}' \
    $qcdir/$alg.pruned.keep.bim \
    $qcdir/$hapmap.updated_flipped.bim > \
    $qcdir/$hapmap.mismatch

plink --bfile $qcdir/$hapmap.updated_flipped \
      --exclude $qcdir/$hapmap.mismatch \
      --make-bed \
      --out $qcdir/$hapmap.clean
mv $qcdir/$hapmap.clean.log \
    $qcdir/plink_log/$hapmap.clean.log

plink --bfile $qcdir/$alg.pruned.keep \
      --exclude $qcdir/$hapmap.mismatch \
      --make-bed \
      --out $qcdir/$alg.hapmap.clean
mv $qcdir/$alg.hapmap.clean.log \
    $qcdir/plink_log/$alg.hapmap.clean.log

# merge HapMap-mapped genotypes to genotypes of HapMap populations
plink --bfile $qcdir/$hapmap.clean  \
      --bmerge $qcdir/$alg.hapmap.clean.bed $qcdir/$alg.hapmap.clean.bim \
         $qcdir/$alg.hapmap.clean.fam  \
      --make-bed --out $qcdir/$alg.merge.$hapmap
mv $qcdir/$alg.merge.$hapmap.log $qcdir/plink_log

# 17. Conduct PCA on the merged data
plink --bfile $qcdir/$alg.merge.$hapmap --pca \
      --out $qcdir/$alg.HapMapIII.pruned.pca
mv $qcdir/$alg.HapMapIII.pruned.pca.log $qcdir/plink_log
