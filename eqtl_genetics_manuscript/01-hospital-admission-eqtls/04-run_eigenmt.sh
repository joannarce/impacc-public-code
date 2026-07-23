#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# 04-run_eigenmt.sh
#
# Run eigenMT per chromosome as a computationally efficient proxy for
# permutation-based control. The MatrixEQTL outputs are first reformatted
# for eigenMT in 04-eigenmt_correction.Rmd; this script subsets each input by
# chromosome, runs eigenMT.py, then 04-eigenmt_correction.Rmd combines the
# per-chromosome results.
# ============================================================================

# ---- CONFIGURABLE VARIABLES (edit these for your environment) --------------
# Analysis results root, the eigenMT-formatted input directory, the per-visit
# output directory, and the path to the eigenMT.py script.
RESULTS_DIR="/path/to/hospital_admission_results"
EIGENMT_PY="/path/to/eigenMT/eigenMT.py"
FMT_DIR="${RESULTS_DIR}/eigenMT_formatted_files"
OUT_DIR="${RESULTS_DIR}/eigenMT_results/visit1"   # results for this visit
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# One-time environment setup (run once, then reuse the env).
# ----------------------------------------------------------------------------
# still inside the base env
# conda create -n eigenmt python=3.11 numpy scipy pandas Scikit-learn
# add any other eigenMT requirements (e.g. h5py) if needed

# source ~/anaconda3/etc/profile.d/conda.sh   # initialize conda for *this* subshell
conda activate eigenmt            # prompt becomes (eigenmt)

mkdir -p "$OUT_DIR"

# formatted input files (produced earlier in the R pipeline)
QTL_IN="$FMT_DIR/impacc_panEQTL_visit1_cis_eQTL_eigenMT.txt"
GEN_IN="$FMT_DIR/impacc_panEQTL_visit1_SNP_eigenMT.txt"
GENPOS_IN="$FMT_DIR/impacc_panEQTL_visit1_snploc_eigenMT.txt"
PHEPOS_IN="$FMT_DIR/impacc_panEQTL_visit1_geneloc_eigenMT.txt"

# chromosomes to process (non-standard chr15, 18, 22 were absent in the GWAS hits)
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 16 17 19 20 21)

# ------------------------------------------------------------------
# Loop over chromosomes and run eigenMT
# ------------------------------------------------------------------
for chr in "${chroms[@]}"; do
  echo "Preparing files for chromosome ${chr} …"

  # --- subset to this chromosome ---
  GEN_CHR="$OUT_DIR/SNP_chr${chr}.txt"
  awk -v chr="${chr}:" '$1 ~ "^"chr {print}' "$GEN_IN"  > "$GEN_CHR"

  GENPOS_CHR="$OUT_DIR/snploc_chr${chr}.txt"
  awk -v chr=$chr 'NR==1 || $2==chr {print}' "$GENPOS_IN" > "$GENPOS_CHR"

  QTL_CHR="$OUT_DIR/cis_eQTL_chr${chr}.txt"
  awk -v chr="${chr}:" 'NR==1 || $1 ~ "^"chr {print}' "$QTL_IN" > "$QTL_CHR"

  # --- run eigenMT ---
  OUT_CHR="$OUT_DIR/chr${chr}_eigenMT.txt"
  echo "Running eigenMT on chromosome ${chr}"
  python -W ignore "$EIGENMT_PY" \
        --CHROM  "$chr" \
        --QTL    "$QTL_CHR" \
        --GEN    "$GEN_CHR" \
        --GENPOS "$GENPOS_CHR" \
        --PHEPOS "$PHEPOS_IN" \
        --OUT    "$OUT_CHR"
done

# ------------------------------------------------------------------
# Clean-up intermediates
# ------------------------------------------------------------------
echo "Cleaning up intermediate files"
rm "$OUT_DIR"/SNP_chr*.txt "$OUT_DIR"/snploc_chr*.txt "$OUT_DIR"/cis_eQTL_chr*.txt
