#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# 04-run_eigenmt.sh
#
# Run eigenMT per chromosome for each severity-by-visit stratum, as a
# computationally efficient proxy for permutation-based multiple-testing
# correction. eigenMT estimates the effective number of independent tests per
# gene; the per-stratum, per-chromosome outputs are combined back into one
# merged eigenMT table in 04-eigenmt_correction.Rmd.
#
# Inputs are the eigenMT-formatted files written by 04-eigenmt_correction.Rmd
# (QTL, genotype, snploc, geneloc). For each chromosome the genotype, snploc
# and QTL files are subset by leading "chr:" key, eigenMT is run, and the
# per-chromosome intermediates are removed afterwards.
#
# Severity groups: case = TG4-5, control = TG1-3. visit_5 control is skipped.
# ============================================================================

# ---- CONFIGURABLE VARIABLES (edit these for your environment) --------------
# results_dir         : root of this track's results.
# EIGEN_FORMATTED_DIR : eigenMT-formatted inputs from 04-eigenmt_correction.Rmd.
# RESULTS_DIR         : where eigenMT per-chromosome outputs are written.
# EIGENMT_PATH        : path to eigenMT.py.
results_dir="/path/to/severity_longitudinal_results"
EIGEN_FORMATTED_DIR="${results_dir}/eigenMT_FDR/eigenMT_formatted_files"
RESULTS_DIR="${results_dir}/eigenMT_FDR/eigenMT_results"
EIGENMT_PATH="/path/to/eigenMT/eigenMT.py"
# ----------------------------------------------------------------------------

mkdir -p "${RESULTS_DIR}"

# Define visits (all 10) and groups (case and control)
visits=("visit_1" "visit_2" "visit_3" "visit_4" "visit_5" "visit_6" "visit_7" "visit_8" "visit_9" "visit_10")
groups=("case" "control")

# List of chromosomes to process (the 18 autosomes with extracted variants;
# X is excluded from eigenMT here).
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 16 17 19 20 21)

# Loop over visits and groups
for visit in "${visits[@]}"; do
  for gr in "${groups[@]}"; do
    # Skip visit_5 control
    if [[ "$visit" == "visit_5" && "$gr" == "control" ]]; then
      echo "Skipping ${visit} ${gr}..."
      continue
    fi

    echo "Processing ${visit} ${gr}..."

    # Input file paths from the eigenMT_formatted_files directory:
    QTL_FILE_in="${EIGEN_FORMATTED_DIR}/${visit}/${gr}/${visit}_${gr}_cis_eQTL_eigenMT.txt"
    GEN_FILE_in="${EIGEN_FORMATTED_DIR}/${visit}/${gr}/${visit}_${gr}_SNP_eigenMT.txt"
    GENPOS_FILE_in="${EIGEN_FORMATTED_DIR}/${visit}/${gr}/${visit}_${gr}_snploc_eigenMT.txt"
    PHEPOS_FILE_in="${EIGEN_FORMATTED_DIR}/${visit}/${gr}/${visit}_${gr}_geneloc_eigenMT.txt"

    # Create output directory for eigenMT results for this visit and group
    OUT_DIR="${RESULTS_DIR}/${visit}/${gr}"
    mkdir -p "${OUT_DIR}"

    # Loop over each chromosome
    for chr in "${chroms[@]}"; do
      echo "  Preparing files for chromosome ${chr}..."

      # Create chromosome-specific GEN file.
      GEN_FILE_chr="${OUT_DIR}/${visit}_${gr}_SNP_chr${chr}.txt"
      awk '$1 ~ /^'"${chr}"':/' ${GEN_FILE_in} > ${GEN_FILE_chr}

      # Create chromosome-specific GENPOS file.
      GENPOS_FILE_chr="${OUT_DIR}/${visit}_${gr}_snploc_chr${chr}.txt"
      awk 'NR==1 || $1 ~ /^'"${chr}"':/' ${GENPOS_FILE_in} > ${GENPOS_FILE_chr}

      # Create chromosome-specific QTL file.
      QTL_FILE_chr="${OUT_DIR}/${visit}_${gr}_cis_eQTL_chr${chr}.txt"
      awk 'NR==1 || $1 ~ /^'"${chr}"':/' ${QTL_FILE_in} > ${QTL_FILE_chr}

      # Set the output file for eigenMT results for this chromosome.
      out_file="${OUT_DIR}/chr${chr}_${visit}_${gr}_eigenMT.txt"
      echo "  Processing chromosome ${chr} for ${visit} ${gr}..."

      # Run eigenMT for the current chromosome.
      python -W ignore "${EIGENMT_PATH}" --CHROM ${chr} \
          --QTL ${QTL_FILE_chr} \
          --GEN ${GEN_FILE_chr} \
          --GENPOS ${GENPOS_FILE_chr} \
          --PHEPOS ${PHEPOS_FILE_in} \
          --OUT ${out_file}
    done

    # Once eigenMT files for all chromosomes are computed, remove the intermediate files.
    echo "Cleaning up intermediate files for ${visit} ${gr}..."
    rm ${OUT_DIR}/${visit}_${gr}_SNP_chr*.txt \
       ${OUT_DIR}/${visit}_${gr}_snploc_chr*.txt \
       ${OUT_DIR}/${visit}_${gr}_cis_eQTL_chr*.txt
  done
done
