#!/bin/bash
#SBATCH --array 1-24
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=parallel
#SBATCH --time=3-0:00
#SBATCH --job-name=bcftools_mpileup
#SBATCH -o /nfs/scratch/oostinto/stdout/bcftools_mpileup.%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/bcftools_mpileup.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###!!!###									    ###!!!####
# determine the number of scaffolds you want to genotype #
#               change --array accordingly			     #
###!!!###									    ###!!!####

#export pahts
export PATH=~/bin/bcftools-1.18:$PATH
export PATH=~/bin/htslib-1.18:$PATH

#variables
SCAFFOLD=${SLURM_ARRAY_TASK_ID}
PROJECT=$1
SET_NEW=$PROJECT'_'$2
TMP_DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/tmp
mkdir -p $TMP_DIR

#genome
REF=$SCRATCH/projects/$PROJECT/resources/reference_genomes/nuclear/Chrysophrys_auratus.v.1.0.all.assembly.units.fasta

BAMLIST=$SCRATCH/projects/$PROJECT/resources/bam_lists/$SET_NEW'_bam.list'
REGION=$( head -n $SCAFFOLD $REF.fai | tail -n 1 | cut -f 1 )

#genotype
bcftools mpileup 	-Ov 																\
					-a 'FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD'	\
					-f $REF 															\
					-r $REGION															\
					-b $BAMLIST															|
bcftools call -Ov -mv > $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'

#update INFO fields
bcftools 	+fill-tags 	$TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'			\
			-Oz -o 		$TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz'		\
			-- -t AC,AF,AN,MAF,NS,AC_Hom,AC_Het			
bgzip --reindex $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz'

rm $TMP_DIR/$REGION'_'$SET_NEW'_raw_tmp1.vcf'
 