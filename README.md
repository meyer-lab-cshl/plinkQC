Genotyping arrays enable the direct measurement of an individuals
    genotype at thousands of markers. plinkQC facilitates genotype quality
    control for genetic association studies as described by Anderson et al
    (2010) Nature Protocols. It wraps around PLINK basic statistics (e.g.
    missing genotyping rates per individual, allele frequencies per genetic
    marker) and relationship functions and generates a per-individual and
    per-marker quality control report. Individuals and markers that fail the
    quality control can subsequently be removed to generate a new, clean
    dataset. Removal of individuals based on relationship status is optimised
    to retain as many individuals as possible in the study.
