. /etc/profile.d/modules.sh
# Load fastqc
module load igmm/apps/FastQC/0.11.9
for sample in P2	P3	P6	P7	P8	P9	P15	P16	P17	P19	P20	P21	
do
cd /exports/eddie/scratch/s1763420/X204SC22020675-Z01-F003/01.RawData/$sample/
for file in $(ls *fq.gz)
do
mkdir fastqc_results_$sample
fastqc $file -o /exports/eddie/scratch/s1763420/X204SC22020675-Z01-F003/01.RawData/$sample/fastqc_results_$sample/ 2> /exports/eddie/scratch/s1763420/cleaned_data/$sample/fastqc_results_$sample/${sample}_fastqc_results.log
done
cd ..;
done
echo 'I am done'