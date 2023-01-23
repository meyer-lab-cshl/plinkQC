# Format anderson et al provided, NCBI36-annotated high-LD file into UCSCF bed format
awk '{print "chr" $1, $2 -1, $3}' high-LD-regions-hg18-NCBI36.txt | \
    sed 's/chr23/chrX/' | sed 's/chr24/chrY/' \
    > high-LD-regions-hg18-NCBI36-tolift.txt

# Use the liftOVer tool at https://genome.ucsc.edu/cgi-bin/hgLiftOver with 'Minimum ratio of bases that must
#' remap' set to 0.5 , 'Allow multiple output regions' ticked and default options otherwise. Do liftover
# into hg39 and hg38

# Reformat UCSC bed output into plink 1.9 --make-set format:
awk '{print $1,$2,$3,FNR}' high-LD-regions-hg38-GRCh38.bed > high-LD-regions-hg38-GRCh38.txt
awk '{print $1,$2,$3,FNR}' high-LD-regions-hg19-GRCh37.bed > high-LD-regions-hg19-GRCh37.txt
