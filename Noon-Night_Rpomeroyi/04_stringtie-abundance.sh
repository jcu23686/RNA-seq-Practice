#!/bin/sh
#SBATCH --partition=batch
#SBATCH --job-name=j_03_stringtie-assembly
#SBATCH --ntasks=1
#SBATCH --mem=10gb 
#SBATCH --time=8:00:00 
#SBATCH --output=STRINGTIE.%j.out
#SBATCH --error=STRINGTIE.%j.err

######################################
### Estimate transcript abundances ###
######################################


ml StringTie/2.2.1-GCC-11.2.0

#Morning
#ARRAY1=(D10_fil_NOrRNAstd_sorted.bam D11_fil_NOrRNAstd_sorted.bam E01_fil_NOrRNAstd_sorted.bam E04_fil_NOrRNAstd_sorted.bam E10_fil_NOrRNAstd_sorted.bam F02_fil_NOrRNAstd_sorted.bam)

#Night
ARRAY2=(D12_fil_NOrRNAstd_sorted.bam E05_fil_NOrRNAstd_sorted.bam E07_fil_NOrRNAstd_sorted.bam F01_fil_NOrRNAstd_sorted.bam F04_fil_NOrRNAstd_sorted.bam F05_fil_NOrRNAstd_sorted.bam F07_fil_NOrRNAstd_sorted.bam F09_fil_NOrRNAstd_sorted.bam)

#Noon
ARRAY3=(D08_fil_NOrRNAstd_sorted.bam E02_fil_NOrRNAstd_sorted.bam E06_fil_NOrRNAstd_sorted.bam E08_fil_NOrRNAstd_sorted.bam F03_fil_NOrRNAstd_sorted.bam F06_fil_NOrRNAstd_sorted.bam)

#Dusk
#ARRAY4=(D09_fil_NOrRNAstd_sorted.bam E03_fil_NOrRNAstd_sorted.bam E09_fil_NOrRNAstd_sorted.bam E11_fil_NOrRNAstd_sorted.bam E12_fil_NOrRNAstd_sorted.bam F08_fil_NOrRNAstd_sorted.bam)


NP="4"
GTF="stringtie_merged.gtf"

#for BAM in "${ARRAY1[@]}"  
#  do
#    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
#    time stringtie -e -B -p $NP -G $GTF -o ballgown/${BASENAME}/${BASENAME}.gtf map_Morn/$BAM
#  done

for BAM in "${ARRAY2[@]}"
  do
    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
    time stringtie -e -B -p $NP -G $GTF -o ballgown/${BASENAME}/${BASENAME}.gtf map_Night/$BAM
  done

for BAM in "${ARRAY3[@]}"
  do
    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
    time stringtie -e -B -p $NP -G $GTF -o ballgown/${BASENAME}/${BASENAME}.gtf map_Noon/$BAM
  done

#for BAM in "${ARRAY4[@]}"
#  do
#    BASENAME=$(echo "$BAM" | cut -f 1 -d '.')
#    time stringtie -e -B -p $NP -G $GTF -o ballgown/${BASENAME}/${BASENAME}.gtf map_Dusk/$BAM
#  done

