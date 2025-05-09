#!/bin/bash
#SBATCH -a 1-382
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=depth_bam
#SBATCH -o /nfs/scratch/oostinto/stdout/depth_bam.%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/depth_bam.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#load modules
export PATH=~/bin/bcftools-1.18:$PATH
export PATH=~/bin/htslib-1.18:$PATH
export PATH=~/bin/samtools-1.18:$PATH

#variables
PROJECT=$1
SET=$2
N=${SLURM_ARRAY_TASK_ID}

#
SET_FILE=/nfs/scratch/oostinto/projects/$PROJECT/resources/sample_lists/$PROJECT'_'$SET'_inds.tsv'
dos2unix $SET_FILE
SAMPLE=$( head -n $N $SET_FILE | tail -n 1 )

BAM=/nfs/scratch/oostinto/projects/$PROJECT/data/bam/$SAMPLE'_1.0.nuclear.bam'
DPT=/nfs/scratch/oostinto/projects/$PROJECT/output/samtools/$SAMPLE'_depth.txt'

samtools depth -a $BAM | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > $DPT
