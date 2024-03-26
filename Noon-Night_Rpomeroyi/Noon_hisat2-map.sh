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

ARRAY=(D08_fil_NOrRNAstd.fastq.gz E02_fil_NOrRNAstd.fastq.gz E06_fil_NOrRNAstd.fastq.gz E08_fil_NOrRNAstd.fastq.gz F03_fil_NOrRNAstd.fastq.gz F06_fil_NOrRNAstd.fastq.gz)

INDEX="DSS3mega"
NP="4"

mkdir map_Noon

for FASTQ in "${ARRAY[@]}"
  do
     BASENAME=$(echo "$FASTQ" | cut -f 1 -d '.')
    time hisat2 -p $NP --dta -x $INDEX -U Noon/$FASTQ -S map_Noon/${BASENAME}.sam
    time samtools view -bS  -@ $NP map_Noon/${BASENAME}.sam > map_Noon/${BASENAME}.bam
    rm map_Noon/${BASENAME}.sam 
    time samtools sort -@ $NP map_Noon/${BASENAME}.bam map_Noon/${BASENAME}_sorted
    rm map_Noon/${BASENAME}.bam
  done

