set -euo pipefail

# ==============================
# Project root
# ==============================
# Set this to the full path of your repository on your system.
# Example (HPC): /exports/eddie/scratch/<user>/Barley_SymPath
PROJ="/path/to/your/repo"

# ==============================
# Threads
# ==============================
# Number of CPU threads used by FastQC/HISAT2/Trimmomatic where supported.
THREADS=16

# ==============================
# Sample IDs
# ==============================
# These should match the per-sample folder names under fastq_raw/ (or your raw FASTQ location).
# Example from your Eddie pipeline: P2 P3 P6 ...
SAMPLES=(
  P2 P3 P6 P7 P8 P9 P15 P16 P17 P19 P20 P21
)

# ==============================
# Directory structure (relative to PROJ)
# ==============================
# Keep these aligned with your repo conventions: data_raw/, figures/, scripts/.
# For sequencing data, it is common to keep FASTQs OUT of GitHub; you can still
# use these folders locally or symlink them.

RAW_DIR="${PROJ}/fastq_raw"          # raw FASTQs (NOT typically committed)
CLEAN_DIR="${PROJ}/fastq_clean"      # trimmed paired/unpaired FASTQs (NOT committed)
QC_DIR="${PROJ}/qc"                  # FastQC outputs (NOT committed)
ALIGN_DIR="${PROJ}/alignments"       # HISAT2 SAM/BAM outputs (NOT committed)
COUNT_DIR="${PROJ}/counts"           # HTSeq count outputs + merged count matrices (counts can be committed)

# Optional: keep final figures consistent with your repo
FIG_DIR="${PROJ}/figures"            # your existing figures/ folder

# ==============================
# Reference files
# ==============================
# You mentioned the FASTA is available as:
# Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa
#
# Put references in genome/ locally (often not committed) OR point to a shared reference path on HPC.

GENOME_DIR="${PROJ}/genome"

FASTA="${GENOME_DIR}/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel.fa"

# HTSeq requires an annotation file (GFF3 or GTF). Update to the exact filename you use.
# Example placeholder (edit this to your file):
GFF3="${GENOME_DIR}/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.59.gff3"

# HISAT2 index output directory + basename
HISAT2_INDEX_DIR="${GENOME_DIR}/hisat2_index"
HISAT2_INDEX_BASE="${HISAT2_INDEX_DIR}/MorexV3"

# ==============================
# Tool locations / environment
# ==============================
# If you are using an HPC module system, you can load modules inside each script
# (recommended) rather than here. Put only paths here that are system-specific.

# Trimmomatic v0.39
TRIMMOMATIC_JAR="/path/to/Trimmomatic-0.39/trimmomatic-0.39.jar"
TRIMMOMATIC_ADAPTERS="/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

# HTSeq environment
# If you use conda:
USE_CONDA_FOR_HTSEQ=1
CONDA_MODULE="roslin/conda/5.3.0"
CONDA_ENV_NAME="htseq"

# ==============================
# FASTQ naming patterns
# ==============================
# Your scripts likely assume paired-end FASTQs named with 1/2.
# If your naming differs (e.g. _R1_ / _R2_), adjust these patterns.
R1_GLOB="*1.fq.gz"
R2_REPLACE_FROM="1.fq.gz"
R2_REPLACE_TO="2.fq.gz"

# Alternative (Illumina-style) patterns (uncomment if needed)
# R1_GLOB="*_R1_*.fastq.gz"
# R2_REPLACE_FROM="_R1_"
# R2_REPLACE_TO="_R2_"

# ==============================
# HTSeq options (match your Methods)
# ==============================
# stranded: "yes" / "no" / "reverse"
HTSEQ_STRANDED="yes"
HTSEQ_FEATURE_TYPE="gene"

# ==============================
# Housekeeping
# ==============================
# Create key directories if they do not exist (safe; no error if they already exist)
mkdir -p "${CLEAN_DIR}" "${QC_DIR}" "${ALIGN_DIR}" "${COUNT_DIR}" "${HISAT2_INDEX_DIR}" "${FIG_DIR}"
