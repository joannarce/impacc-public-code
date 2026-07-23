#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# 00-extract_gwas_loci.sh
#
# Severity-stratified, longitudinal cis-eQTL track (TG1-3 vs TG4-5 across
# nine visits).
#
# Extract the genotype data falling within 1 Mb (+/- 1,000,000 bp) of every
# published critical COVID-19 GWAS lead SNP (GenOMICC PMID 37198478 and
# HGI PMID 37674002). The range list (ranges_to_extract.txt) is written by
# 01-genotype_processing.Rmd. This script consumes it with PLINK
# --extract range and produces the locus-restricted BED/BIM/FAM used
# downstream by MatrixEQTL.
#
# This does not need to be re-run for every analysis: the extracted
# genotype set is a stable pipeline input.
# ============================================================================

# ---- CONFIGURABLE VARIABLES (edit these for your environment) --------------
# PLINK1.9/2.0 binary prefix for the full post-QC IMPACC genotype set, the
# range file produced in R, the single-position SNP list, and the output
# directory.
GENO_PREFIX="/path/to/impacc_geno"
RESULTS_DIR="/path/to/severity_longitudinal_results"
RANGES="${RESULTS_DIR}/ranges_to_extract.txt"
SNPS="${RESULTS_DIR}/snps_to_extract.txt"
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Reference / historical note on the published-SNP handling.
#
# The published list is 69 non-overlapping critical COVID-19 GWAS SNPs
# (snps_to_extract.txt). An earlier, single-position extraction against the
# imputed reference panel found 58 of the 69 SNPs in the IMPACC cohort, and
# the 11 unmatched SNPs were written out for follow-up. That step is retained
# here for provenance.
# ----------------------------------------------------------------------------
plink2 --bfile "${GENO_PREFIX}" \
       --extract range "${SNPS}" \
       --make-bed \
       --out "${RESULTS_DIR}/extracted_snps"
# 58 of the 69 SNPs are found in our data.

# ----------------------------------------------------------------------------
# Active extraction: pull all variants within the +/- 1 Mb ranges around the
# published SNPs. 299,380 variants were extracted across the 69 ranges.
# ----------------------------------------------------------------------------
plink2 --bfile "${GENO_PREFIX}" \
       --extract range "${RANGES}" \
       --make-bed \
       --out "${RESULTS_DIR}/extracted_ranges"

# ============================================================================
# Per-stratum genotype extraction (run AFTER expression processing).
#
# 02-expression_processing.Rmd writes one PLINK keep-list of post-QC sample
# IDs per severity-by-visit stratum (named "*_coreIDs.txt", under
# genotype_data/visit_*/{case,control}/). The loop below restricts the
# genotype set to each stratum's samples within the published GWAS ranges,
# producing the three PLINK outputs consumed downstream:
#   *_eQTL.traw : A-transpose genotype dosages for MatrixEQTL (03)
#   *_pca.bed   : binary set for PC-AiR genotype PCs (01)
#   *_ped.map   : map file used to build the snploc file (03)
#
# Original IMPACC server version sourced biogrids and used PLINK1.9 (plink);
# the variable names below replace the hardcoded paths.
# ----------------------------------------------------------------------------

# ---- ADDITIONAL CONFIGURABLE VARIABLES -------------------------------------
# Directory holding the per-stratum "*_coreIDs.txt" keep-lists written by R
# (genotype_data, with visit_*/{case,control} subfolders).
SAMPLES_DIR="${RESULTS_DIR}/genotype_data"
# ----------------------------------------------------------------------------

# Loop over each visit folder (e.g., "visit_1", "visit_2", etc.)
for visit_dir in "$SAMPLES_DIR"/visit_*; do
    # Loop over each group folder (e.g., "case" and "control") within the visit folder
    for group_dir in "$visit_dir"/*; do
        # Find the sample extraction file (assumed to match "*_coreIDs.txt")
        sample_file=$(find "$group_dir" -maxdepth 1 -type f -name "*_coreIDs.txt")

        # Extract visit and group basenames for naming output files
        visit_basename=$(basename "$visit_dir")  # e.g., visit_1
        group_basename=$(basename "$group_dir")  # e.g., case or control

        # Define output prefixes in the same group folder using a concise naming scheme
        out_prefix_eQTL="${group_dir}/${visit_basename}_${group_basename}_eQTL"
        out_prefix_pca="${group_dir}/${visit_basename}_${group_basename}_pca"
        out_prefix_ped="${group_dir}/${visit_basename}_${group_basename}_ped"

        # PLINK Run 1: Extract variants in A-transpose format
        plink --bfile "$GENO_PREFIX" \
              --keep "$sample_file" \
              --extract range "$RANGES" \
              --recode A-transpose \
              --out "$out_prefix_eQTL"

        # Check if the .traw file was generated
        if [ ! -f "${out_prefix_eQTL}.traw" ]; then
            echo "Error: ${out_prefix_eQTL}.traw not generated. Stopping execution."
            exit 1
        fi

        # PLINK Run 2: Create PLINK binary files for PCA (--make-bed)
        plink --bfile "$GENO_PREFIX" \
              --keep "$sample_file" \
              --make-bed \
              --out "$out_prefix_pca"

        # PLINK Run 3: Convert to ped/map format for downstream snpsloc file.
        plink --bfile "$GENO_PREFIX" \
              --keep "$sample_file" \
              --extract range "$RANGES" \
              --recode \
              --out "$out_prefix_ped"

        # Remove the .ped file to save space, leaving only the .map file.
        if [ -f "${out_prefix_ped}.ped" ]; then
            rm "${out_prefix_ped}.ped"
        fi

        # Optionally remove any .nosex files generated by PLINK
        for f in "${out_prefix_pca}.nosex" "${out_prefix_ped}.nosex" "${out_prefix_eQTL}.nosex"; do
          if [ -f "$f" ]; then
            rm "$f"
          fi
        done

    done
done

# ============================================================================
# Per-stratum allele-frequency extraction for colocalization (run during 05).
#
# 05-colocalization.Rmd writes ranges_for_coloc.txt (+/- 500 kb around the
# significant lead SNPs). The loop below restricts each stratum's genotypes to
# those ranges and computes minor-allele frequencies (--freq) plus a .bim, used
# to harmonize the eQTL summary statistics with the GenOMICC GWAS in 05.
# ----------------------------------------------------------------------------

# ---- ADDITIONAL CONFIGURABLE VARIABLES -------------------------------------
# Range file for coloc windows and the coloc output directory (per-stratum
# .frq / .bim are written under coloc_dir/visit_*/{case,control}/).
RANGES_COLOC="${RESULTS_DIR}/genotype_data/ranges_for_coloc.txt"
COLOC_DIR="${RESULTS_DIR}/coloc"
# ----------------------------------------------------------------------------

# Loop over each visit folder (e.g., "visit_1", "visit_2", etc.)
for visit_dir in "$SAMPLES_DIR"/visit_*; do
    visit_basename=$(basename "$visit_dir")  # e.g., visit_1

    # Create corresponding visit subdirectory in COLOC_DIR
    mkdir -p "$COLOC_DIR/$visit_basename"

    # Loop over each group folder ("case" and "control") within each visit folder
    for group_dir in "$visit_dir"/*; do
        group_basename=$(basename "$group_dir")  # e.g., case or control

        # Create corresponding group subdirectory in COLOC_DIR
        mkdir -p "$COLOC_DIR/$visit_basename/$group_basename"

        # Find the sample extraction file matching "*_coreIDs.txt"
        sample_file=$(find "$group_dir" -maxdepth 1 -type f -name "*_coreIDs.txt")

        # Define output prefix for allele frequency inside COLOC_DIR
        out_prefix_freq="$COLOC_DIR/$visit_basename/$group_basename/${visit_basename}_${group_basename}"

        # PLINK run: Calculate allele frequencies (--freq)
        plink --bfile "$GENO_PREFIX" \
              --keep "$sample_file" \
              --extract range "$RANGES_COLOC" \
              --freq \
              --make-bed \
              --out "$out_prefix_freq"

        # Check if .frq file was generated
        if [ ! -f "${out_prefix_freq}.frq" ]; then
            echo "Error: ${out_prefix_freq}.frq not generated. Stopping execution."
            exit 1
        fi

        # Remove unnecessary files generated by PLINK
        for file_ext in .nosex .bed .fam; do
          if [ -f "${out_prefix_freq}${file_ext}" ]; then
              rm "${out_prefix_freq}${file_ext}"
          fi
        done

    done
done
