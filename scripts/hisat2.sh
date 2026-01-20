#!/bin/bash
# Initialise the environment modules
. /etc/profile.d/modules.sh

for sample in P2	P3	P6	P7	P8	P9	P15	P16	P17	P19	P20	P21
do
cd /exports/eddie/scratch/s1763420/RNAseq_data/X204SC22020675-Z01-F001/cleaned_data/$sample/
for R1 in *1_paired.fq.gz
do
R2=${R1//1_paired.fq.gz/2_paired.fq.gz}
module load igmm/apps/HISAT2/2.1.0
hisat2 -q -x /exports/eddie/scratch/s1763420/genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.toplevel -1 $R1 -2 $R2 -S ${sample}.sam --summary-file ${sample}__alignmentSummary.txt
done
done
