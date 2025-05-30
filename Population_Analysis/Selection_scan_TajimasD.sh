#!/bin/bash
#SBATCH -a 1-24
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=parallel
#SBATCH --time=5-0:00
#SBATCH --job-name=VCFtools_TjD
#SBATCH -o /nfs/scratch/oostinto/stdout/VCFtools_TjD_%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/VCFtools_TjD_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

source activate /nfs/scratch/oostinto/conda/vcftools

###run input
PROJECT=$1
SET=$PROJECT'_'$2
CLUSTER=$3
LG=LG${SLURM_ARRAY_TASK_ID}

VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
OUT=$SCRATCH/projects/$PROJECT/output/$SET/vcftools/$LG'_'$SET'_'$CLUSTER
mkdir $SCRATCH/projects/$PROJECT/output/$SET/vcftools

vcftools 	--gzvcf $VCF'_qc.vcf.gz' 		\
			--keep ./$SET'_'$CLUSTER.tsv	\
			--chr $LG						\
			--out $OUT						\
			--TajimaD 5000