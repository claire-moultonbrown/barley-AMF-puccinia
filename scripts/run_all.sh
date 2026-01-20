#!/usr/bin/env bash
set -euo pipefail

bash scripts/01_fastqc_raw.sh
bash scripts/02_trimmomatic.sh
bash scripts/03_fastqc_clean.sh
bash scripts/04_hisat2_build_index.sh
bash scripts/05_hisat2_align.sh
bash scripts/06_htseq_count.sh
bash scripts/07_collect_counts.sh
