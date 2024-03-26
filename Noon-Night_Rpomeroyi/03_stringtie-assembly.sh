#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=j_03_stringtie-assembly
#SBATCH --ntasks=1
#SBATCH --mem=10gb 
#SBATCH --time=8:00:00 
#SBATCH --output=STRINGTIE.%j.out
#SBATCH --error=STRINGTIE.%j.err

###########################################################################
### Assemble the mapped reads into transcripts guided by the annotation ###
###########################################################################

ml StringTie/2.2.1-GCC-11.2.0

#Morning
#ARRAY1=(D10_fil_NOrRNAstd.fastq.gz D11_fil_NOrRNAstd.fastq.gz E01_fil_NOrRNAstd.fastq.gz E04_fil_NOrRNAstd.fastq.gz E10_fil_NOrRNAstd.fastq.gz F02_fil_NOrRNAstd.fastq.gz)

#Night
ARRAY2=(D12_fil_NOrRNAstd.fastq.gz E05_fil_NOrRNAstd.fastq.gz E07_fil_NOrRNAstd.fastq.gz F01_fil_NOrRNAstd.fastq.gz F04_fil_NOrRNAstd.fastq.gz F05_fil_NOrRNAstd.fastq.gz F07_fil_NOrRNAstd.fastq.gz F09_fil_NOrRNAstd.fastq.gz)

#Noon
ARRAY3=(D08_fil_NOrRNAstd.fastq.gz E02_fil_NOrRNAstd.fastq.gz E06_fil_NOrRNAstd.fastq.gz E08_fil_NOrRNAstd.fastq.gz F03_fil_NOrRNAstd.fastq.gz F06_fil_NOrRNAstd.fastq.gz)

#Dusk
#ARRAY4=(D09_fil_NOrRNAstd.fastq.gz E03_fil_NOrRNAstd.fastq.gz E09_fil_NOrRNAstd.fastq.gz E11_fil_NOrRNAstd.fastq.gz E12_fil_NOrRNAstd.fastq.gz F08_fil_NOrRNAstd.fastq.gz)

NP="4"
GFF="DSS3mega.gff"

mkdir assembly

#for BAM in "${ARRAY1[@]}"
#  do
#    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
#    time stringtie map_Morn/${BASENAME}_sorted.bam -l $BASENAME -p $NP -G $GFF -o assembly/${BASENAME}.gtf
#    echo assembly/${BASENAME}.gtf > mergelist.txt
#  done

for BAM in "${ARRAY2[@]}"
  do
    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
    time stringtie map_Night/${BASENAME}_sorted.bam -l $BASENAME -p $NP -G $GFF -o assembly/${BASENAME}.gtf
    echo assembly/${BASENAME}.gtf > mergelist.txt
  done

for BAM in "${ARRAY3[@]}"
  do
    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
    time stringtie map_Noon/${BASENAME}_sorted.bam -l $BASENAME -p $NP -G $GFF -o assembly/${BASENAME}.gtf
    echo assembly/${BASENAME}.gtf >> mergelist.txt
  done

#for BAM in "${ARRAY4[@]}"
#  do
#    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
#    time stringtie map_Dusk/${BASENAME}_sorted.bam -l $BASENAME -p $NP -G $GFF -o assembly/${BASENAME}.gtf
#    echo assembly/${BASENAME}.gtf >> mergelist.txt
#  done

time stringtie --merge -p $NP -G $GFF -o stringtie_merged.gtf mergelist.txt

