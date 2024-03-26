#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=j_02_hisat2-map
#SBATCH --ntasks=1
#SBATCH --mem=10gb 
#SBATCH --time=8:00:00 
#SBATCH --output=HISAT2.%j.out
#SBATCH --error=HISAT2.%j.err

###############################################
### Map RNA-seq to indexed reference genome ###
###############################################

ml HISAT2/2.2.1-gompi-2022a
ml SAMtools/1.16.1-GCC-11.3.0

ARRAY=(D12_fil_NOrRNAstd.fastq.gz E05_fil_NOrRNAstd.fastq.gz E07_fil_NOrRNAstd.fastq.gz F01_fil_NOrRNAstd.fastq.gz F04_fil_NOrRNAstd.fastq.gz F05_fil_NOrRNAstd.fastq.gz F07_fil_NOrRNAstd.fastq.gz F09_fil_NOrRNAstd.fastq.gz)

INDEX="DSS3mega"
NP="4"

mkdir map_Night

for FASTQ in "${ARRAY[@]}"
  do
     BASENAME=$(echo "$FASTQ" | cut -f 1 -d '.')
    time hisat2 -p $NP --dta -x $INDEX -U Night/$FASTQ -S map_Night/${BASENAME}.sam
    time samtools view -bS  -@ $NP map_Night/${BASENAME}.sam > map_Night/${BASENAME}.bam
    rm map_Night/${BASENAME}.sam 
    time samtools sort -@ $NP map_Night/${BASENAME}.bam map_Night/${BASENAME}_sorted
    rm map_Night/${BASENAME}.bam
  done

