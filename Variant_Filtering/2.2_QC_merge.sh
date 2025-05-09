#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --partition=parallel
#SBATCH --time=0-5:00
#SBATCH --job-name=merge_vcf_bcftools
#SBATCH -o /nfs/scratch/oostinto/stdout/merge_vcf_bcftools.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/merge_vcf_bcftools.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load bcftools/1.10.1
#module load R/4.0.2
#module load GCC/10.3.0
#module load OpenMPI/4.1.1
#module load R/4.1.0
Rscript=/nfs/home/oostinto/bin/R-4.3.1/bin/Rscript


#script
vcf2R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R/vcf2Rinput.R
gds2plink=$SCRATCH/scripts/population_genomics/4_variant_filtering/R/gds2plink.R

#pop file
POP=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET'_qc'
DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF.vcf.gz $( ls -v $DIR/*'_tmp5_qc.vcf.gz' ) --threads 10
bgzip --reindex $VCF.vcf.gz
tabix -p vcf $VCF.vcf.gz

#create gds
$Rscript $vcf2R 	--gzvcf $VCF.vcf.gz 	\
					--snprelate_out $VCF	
#					--genlight_out $VCF

#create plink bed file
$Rscript $gds2plink 	--gds_file	$VCF.gds 	\
						--out 		$VCF		\
						--pop_file	$POP

#clean up
#rm -r $DIR


##3#SBATCH --cpus-per-task=10
##3#SBATCH --mem=20G
##3#SBATCH --partition=parallel