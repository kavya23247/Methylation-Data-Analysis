!#/bin/bash
#for tumor samples
./bowtie2 -x bowtie2-2.5.1-linux-x86_64/hg19/hg19 \
     bowtie2-2.5.1-linux-x86_64/PA220KH-lib09-P19-Tumor_S2_L001_R1_001.fastq \
     bowtie2-2.5.1-linux-x86_64/PA220KH-lib09-P19-Tumor_S2_L001_R2_001.fastq \
     bowtie2-2.5.1-linux-x86_64/PA220KH-lib09-P19-Tumor_S2_aligned.sam
#for normal samples
./bowtie2 -x bowtie2-2.5.1-linux-x86_64/hg19/hg19 \
     bowtie2-2.5.1-linux-x86_64/PA221MH-lib09-P19-Norm_S1_L001_R1_001.fastq \
     bowtie2-2.5.1-linux-x86_64/PA221MH-lib09-P19-Norm_S1_L001_R2_001.fastq \
     bowtie2-2.5.1-linux-x86_64/PA221MH-lib09-P19-Norm_S1_aligned.sam
# Convert tumor SAM to BAM
samtools view -bS PA220KH-lib09-P19-Tumor_S2_aligned.sam > PA220KH-lib09-P19-Tumor_S2_aligned.bam

# Convert normal SAM to BAM
samtools view -bS PA221MH-lib09-P19-Norm_S1_aligned.sam > PA221MH-lib09-P19-Norm_S1_aligned.bam

# Sort the BAM files
samtools sort PA220KH-lib09-P19-Tumor_S2_aligned.bam -o PA220KH-lib09-P19-Tumor_S2_sorted.bam
samtools sort PA221MH-lib09-P19-Norm_S1_aligned.bam -o PA221MH-lib09-P19-Norm_S1_sorted.bam

# Index the BAM files
samtools index PA220KH-lib09-P19-Tumor_S2_sorted.bam
samtools index PA221MH-lib09-P19-Norm_S1_sorted.bam



