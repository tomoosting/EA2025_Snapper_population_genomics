#!/bin/bash
#SBATCH -a 1-211
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=bcf_stats
#SBATCH -o /nfs/scratch/oostinto/stdout/bcf_stats.%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/bcf_stats.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#load modules
module load htslib/1.9
module load samtools/1.10

#variables
PROJECT=$1
SET=$PROJECT'_'$2
i=${SLURM_ARRAY_TASK_ID}

#paths
BAMLIST=$SCRATCH/projects/$PROJECT/resources/sample_info/$SET'_bam.list'
OUT_DIR=$SCRATCH/projects/$PROJECT/output/$SET/bam_stats
mkdir -p $OUT_DIR

#bam
BAM=$( head -n $i $BAMLIST | tail -n 1 )
SAMPLE=$( basename $BAM | cut -d_ -f1 )

##generate samtools stats
samtools stats $BAM > $OUT_DIR/$SAMPLE'_stats.txt'
##depth
DP=$( samtools depth -a $BAM | awk '{sum+=$3} END { print sum/NR}' )
printf "$SAMPLE\t$DP\r\n" >> $OUT_DIR/depth_list.tsv
