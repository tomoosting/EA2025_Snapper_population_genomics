#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=filter_neutral
#SBATCH -o /nfs/scratch/oostinto/stdout/filter_neutral.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/filter_neutral.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1
#module load R/4.0.2
module load sf/0.9-5-R-4.0.0-Python-3.8.2
module load plink/1.90


#Rscript=/nfs/home/oostinto/bin/R-4.3.1/bin/Rscript

#Rscript paths
R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R

###resrouces
pop_file=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
TMP=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp
PCADAPT=$SCRATCH/projects/$PROJECT/output/$SET/outlier_analyses/pcadapt

### create directories
mkdir $TMP
#parameters you used to for outlier detection using pcadapt to remove the correct outliers
QVAL=0.05
K=1
##################################### filter neutral loci  #############################
if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' ]
then
		tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_all_outlier_LOC.tsv'
		vcftools	--gzvcf $VCF'_qc.vcf.gz' 								\
					--out $TMP/$SET'_tmp6'									\
					--exclude-positions $PCADAPT/$SET'_all_outlier_LOC.tsv' \
					--hwe 0.05 --maf 0.05 --thin 5000						\
					--recode-INFO-all --recode
else
		vcftools	--gzvcf $VCF'_qc.vcf.gz' 			\
					--out $TMP/$SET'_tmp6'				\
					--hwe 0.05 --maf 0.05 --thin 5000	\
					--recode-INFO-all --recode		
fi
bgzip -ci $TMP/$SET'_tmp6.recode.vcf' > $TMP/$SET'_tmp6.vcf.gz'
#get vcf.gz
$Rscript $vcf2R 	--gzvcf $TMP/$SET'_tmp6.vcf.gz' 	\
					--snprelate_out $TMP/$SET'_tmp6'
#get gds
$Rscript $gds2plink 	--gds_file	$TMP/$SET'_tmp6.gds'	\
						--out 		$TMP/$SET'_tmp6'		\
						--pop_file	$pop_file
#filter linked sites
plink 	--bfile $TMP/$SET'_tmp6' 	\
		--indep-pairwise 50 5 0.2 	\
		--out $TMP/$SET'_tmp6'
sed -i 's/:/\t/g' $TMP/$SET'_tmp6.prune.in'
# extract independant SNPs
vcftools 	--gzvcf $TMP/$SET'_tmp6.vcf.gz' 		\
			--out $VCF'_neutral'					\
			--positions $TMP/$SET'_tmp6.prune.in' 	\
			--recode-INFO-all --recode
#convert to vcf.gz
mv $VCF'_neutral.recode.vcf' $VCF'_neutral.vcf'
bgzip -i $VCF'_neutral.vcf'
#convert to gds
$Rscript $vcf2R 		--gzvcf $VCF'_neutral.vcf.gz' 	\
						--snprelate_out $VCF'_neutral'	\
						--genlight_out  $VCF'_neutral'
#convert to plink					
$Rscript $gds2plink 	--gds_file	$VCF'_neutral.gds'	\
						--out 		$VCF'_neutral'		\
						--pop_file	$pop_file
##clean up
#rm -r $TMP/
