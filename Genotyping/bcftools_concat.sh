#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --partition=parallel
#SBATCH --time=2-0:00
#SBATCH --job-name=bcf_concat
#SBATCH -o /nfs/scratch/oostinto/stdout/bcf_concat.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/bcf_concat.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#load modules
module load htslib/1.9
module load vcftools/0.1.16
#module load R/4.0.2
module load bcftools/1.10.1

Rscript=/nfs/home/oostinto/bin/R-4.3.1/bin/Rscript
vcf2R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R/vcf2Rinput.R

#variables
PROJECT=$1
SET_NEW=$PROJECT'_'$2

VCF_NEW=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/$SET_NEW'_raw'
TMP_DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET_NEW/tmp

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF_NEW'.vcf.gz' $( ls -v $TMP_DIR/*'_raw_tmp2.vcf.gz' ) --threads 10
bgzip --reindex $VCF_NEW'.vcf.gz'
tabix -p vcf $VCF_NEW'.vcf.gz'

$Rscript $vcf2R 	--gzvcf 		$VCF_NEW.vcf.gz \
					--snprelate_out $VCF_NEW

#remove tmp folder when you don't need it anymore
#rm -r $TMP_DIR
