# barley-AMF-puccinia
Data and R code for Moulton-Brown, Brzezinska, Orosa Puenta & Helgason "AMF primes immune genes against Puccinia hordei (Brown rust) in Hordeum vulgare but does not reduce pathogen burden".

A preprint of the manuscript and supplementary material can be found on bioRxiv here: https://www.biorxiv.org/content/10.1101/2025.09.16.676302v1 

This repository contains reproducible scripts for processing barley (Hordeum vulgare, MorexV3) RNA-seq data and performing differential expression and downstream analyses consistent with the manuscript Methods.

The workflow follows:

1. Raw read QC with FastQC (v0.11.9)
2. Adapter and quality trimming with Trimmomatic (v0.39)
3. QC of cleaned reads with FastQC (v0.11.9)
4. Alignment to MorexV3 with HISAT2 (v2.1.0)
5. Gene-level counting with HTSeq (v2.0.5)
6. Differential expression with edgeR (v3.42.4), using CPM filtering and TMM normalisation
7. GO enrichment of DEGs with g:Profiler (FDR correction)
8. Hierarchical clustering of DEGs using Ward’s method

# Home directory
Contents:
- this README
- barley-amf-puccinia.Rproj
  - R project file. Download the entire repository and open this file in R studio to run all R scripts with the correct relative file paths.
-.gitignore
- LICESNSE  
-scripts
  - 00_config.example.sh
  - 01_fastqc_raw.sh
  - 02_trimmomatic.sh
  - 03_fastqc_clean.sh
  - 04_hisat2_build_index.sh
  - 05_hisat2_align.sh
  - 06_htseq_count.sh
  - 07_collect_counts.sh
  - run_all.sh
- r_scripts
  - 01_edger.dge.R
  - 02_go_profiler.R
  - figure_1_experimental_design.R
  - figure_2_endpoint_measurements.R
  - figure_3_qPCR_defence_genes.R
  - figure_4_transcriptomics.R
  - figure_5_wrkys_and_cdpks.R
  - figure_6_ubiquitination.R
  - figure_s1_plant_biomass_allocation.R
- data_raw
  - meta_data.csv
  - fig2_endpoint_measurements.csv
  - fig3_qpcr_defence_genes.csv
  - fig4a_transcriptomics_volcano.csv
  - fig4b_transcriptomics_amf_induced_go_terms.csv
  - fig4c_transcriptomics_expression_by_treatment.csv
  - fig4c_transcriptomics_cluster_annotation.csv
  - fig4d_transcriptomics_enriched_cluster_go_terms.csv
  - fig5a_transcriptomics_wrky_expression.csv
  - fig5b_transcriptomics_wrky_go_terms.csv
  - fig5c_transcriptomics_cdpk_expression.csv
  - fig6a_ubiquitination_western_data.csv
  - fig6b_ubiquitination_expression.csv
  - fig6b_ubiquitination_annotation.csv 
- results
  - counts
  - deg
  - go
- figures
  - fig1_experimental_plan_fig.png
  - fig2a_root_annotated.png
  - fig6a_ubiquitination_western.png

# Requirements

You will need the following software available on your system (via modules, conda, or PATH):

- FastQC v0.11.9
- Trimmomatic v0.39
- HISAT2 v2.1.0
- HTSeq v2.0.5
- R (>= 4.2 recommended) with:
  - edgeR (v3.42.4)
  - tidyverse
  - gprofiler2

