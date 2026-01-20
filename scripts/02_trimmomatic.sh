for sample in P3	P6	P7	P8	P9	P15	P16	P17	P19	P20	P21
do
cd /exports/eddie/scratch/s1763420/X204SC22020675-Z01-F003/01.RawData/$sample/
for R1 in *1.fq.gz
do
R2=${R1//1.fq.gz/2.fq.gz}
R1paired=${R1//.fq.gz/_paired.fq.gz}
R1unpaired=${R1//.fq.gz/_unpaired.fq.gz}
R2paired=${R2//.fq.gz/_paired.fq.gz}
R2unpaired=${R2//.fq.gz/_unpaired.fq.gz}
mkdir /exports/eddie/scratch/s1763420/cleaned_data/$sample/
java -jar /exports/eddie/scratch/s1763420/Programmes/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $R1 $R2 /exports/eddie/scratch/s1763420/cleaned_data/$sample/$R1paired /exports/eddie/scratch/s1763420/cleaned_data/$sample/$R1unpaired \
/exports/eddie/scratch/s1763420/cleaned_data/$sample/$R2paired /exports/eddie/scratch/s1763420/cleaned_data/$sample/$R2unpaired ILLUMINACLIP:/exports/eddie/scratch/s1763420/Programmes/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36 \
2> /exports/eddie/scratch/s1763420/cleaned_data/$sample/${sample}.log
done
done
