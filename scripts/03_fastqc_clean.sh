. /etc/profile.d/modules.sh
# Load fastqc
module load igmm/apps/FastQC/0.11.9
for sample in P3	P9	P20	P21
do
cd /exports/eddie/scratch/s1763420/cleaned_data/$sample/
for file in $(ls *fq.gz)
do
mkdir fastqc_results_$sample
fastqc $file -o /exports/eddie/scratch/s1763420/cleaned_data/$sample/fastqc_results_$sample/ 2> /exports/eddie/scratch/s1763420/cleaned_data/$sample/fastqc_results_$sample/${sample}_fastqc_results.log -t 16
done
cd ..;
done
echo 'I am done'
