#!/bin/bash
#SBATCH --job-name=j_HISAT2
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=HISAT2.%j.out
#SBATCH --error=HISAT2.%j.err
cd $SLURM_SUBMIT_DIR
ml HISAT2/2.2.1-gompi-2022a

REF="DSS3mega.fna" #genome reference fasta file
NP="4" #number of processors

BASENAME=$(echo "$REF" | cut -f 1 -d '.')

time hisat2-build \
-p $NP \
$REF \
$BASENAME
