. /etc/profile.d/modules.sh
GFF_FILE=/exports/eddie/scratch/s1763420/Genome/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.59.gff3
module load roslin/conda/5.3.0
eval "$(conda shell.bash hook)"
conda activate htseq
for sample in P3	P6	P7	P8	P9	P15	P16	P17	P19	P20	P21
do
cd /exports/eddie/scratch/s1763420/cleaned_data/$sample/
htseq-count /exports/eddie/scratch/s1763420/cleaned_data/$sample/${sample}.sam -t gene -s yes $GFF_FILE > ${sample}_counts.txt
done
