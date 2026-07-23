#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# 00-extract_gwas_loci.sh
#
# Extract the genotype data falling within 1 Mb (+/- 1,000,000 bp) of every
# published critical COVID-19 GWAS lead SNP (GenOMICC PMID 37198478 and
# HGI PMID 37674002). The range list (ranges_to_extract.txt) is written by
# 01-genotype_processing.Rmd. This script uses it with PLINK
# --extract range and produces the locus-restricted BED/BIM/FAM used
# downstream by MatrixEQTL.
#
# This does not need to be re-run for every analysis: the extracted
# genotype set is a stable pipeline input.
# ============================================================================

# ---- CONFIGURABLE VARIABLES (edit these for your environment) --------------
# PLINK1.9/2.0 binary prefix for the full post-QC IMPACC genotype set,
# the range file produced in R, and the output directory.
GENO_PREFIX="/path/to/impacc_geno"
RESULTS_DIR="/path/to/hospital_admission_results"
RANGES="${RESULTS_DIR}/ranges_to_extract.txt"
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Reference / historical note on the published-SNP handling.
#
# The published list is 69 non-overlapping critical COVID-19 GWAS SNPs
# (snps_to_extract.txt). An earlier, single-position extraction against the
# imputed reference panel found 58 of the 69 SNPs in the IMPACC cohort, and
# the 11 unmatched SNPs were written out for follow-up. That step is retained
# here, commented, for provenance.
#
# # Source biogrids
# source /programs/biogrids.shrc
#
# genetics_dir="/path/to/"
# base_dir="/path/to/impacc_allGWAS_eQTLs_TG/"
#
# plink --bfile ${genetics_dir}/impacc_chrALL_imputed_rsq07_ancestry_maf05_38 \
#       --extract range ${base_dir}/published_gwas_snps/snps_to_extract.txt \
#       --make-bed \
#       --out ${base_dir}/genotype_data/extracted_snps
#
# # 58 of the 69 SNPs are found in our data
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Active extraction: pull all variants within the +/- 1 Mb ranges around the
# published SNPs. 299,380 variants were extracted across the 69 ranges.
# ----------------------------------------------------------------------------
plink2 --bfile "${GENO_PREFIX}" \
       --extract range "${RANGES}" \
       --make-bed \
       --out "${RESULTS_DIR}/extracted_ranges"
