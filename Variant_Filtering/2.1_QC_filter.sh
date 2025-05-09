#!/bin/bash
#SBATCH --array 1-24
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=parallel
#SBATCH --time=1-0:00
#SBATCH --job-name=QC_filtering
#SBATCH -o /nfs/scratch/oostinto/stdout/QC_filtering.%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/QC_filtering.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2
LG=LG${SLURM_ARRAY_TASK_ID}

###load packages
#module load GCC/10.3.0
#module load OpenMPI/4.1.1
#module load R/4.1.0
#module load R/4.0.2
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1

Rscript=/nfs/home/oostinto/bin/R-4.3.1/bin/Rscript
#Rscript paths
AB_script=$SCRATCH/scripts/population_genomics/4_variant_filtering/R/allelelic_imbalance_4.0.R

###resrouces
REF=$SCRATCH/projects/$PROJECT/resources/reference_genomes/nuclear/Chrysophrys_auratus.v.1.0.all.assembly.units.fasta
AB_exclude=$SCRATCH/projects/$PROJECT/resources/sample_info/high_coverage_samples.list

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
DIR=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp
mkdir -p $DIR
TMP=$DIR/$LG'_'$SET

###subsample vcf
bcftools view -Oz -o $TMP'_tmp1.vcf.gz' $VCF'_raw.vcf.gz' -r $LG

##filter genotyeps with DP < 3
vcftools	--gzvcf $TMP'_tmp1.vcf.gz' 	\
			--out 	$TMP'_tmp2'	 		\
			--minDP 3 					\
			--remove-indels  			\
			--recode-INFO-all --recode
mv $TMP'_tmp2.recode.vcf' $TMP'_tmp2.vcf'
bgzip $TMP'_tmp2.vcf'

###basic filter parameters
vcftools 	--gzvcf $TMP'_tmp2.vcf.gz'					\
			--out $TMP'_tmp3' 	  --max-missing 0.95	\
			--min-alleles 2       --max-alleles 2 		\
			--min-meanDP  8       --max-meanDP 25 		\
			--minQ 600 			  --maf 0.01			\
			--recode-INFO-all 	  --recode
mv $TMP'_tmp3.recode.vcf' $TMP'_tmp3.vcf'
bgzip $TMP'_tmp3.vcf' 

###allelic imbalance
#Select output from VCF (genotypes)
vcftools 	--gzvcf $TMP'_tmp3.vcf.gz' 	\
			--out 	$TMP'_tmp4'			\
			--extract-FORMAT-info GT 	
#Select output from VCF (allelic depth)
vcftools 	--gzvcf $TMP'_tmp3.vcf.gz' 	\
			--out 	$TMP'_tmp4'			\
			--extract-FORMAT-info AD
#run binomial test to filter sites with allelic imbalance - could require high mem when many SNPs are to be analysed
$Rscript 	$AB_script 	--GT_file  $TMP'_tmp4.GT.FORMAT' 	\
						--AD_file  $TMP'_tmp4.AD.FORMAT' 	\
						--out_file $TMP'_tmp4_qc' 			\
						--conf.level 0.99			 		\
						--plots TRUE 						\
						--remove $AB_exclude	### make a list of sample names you want to exclude
#filter sites for allelic imbalance - AB
vcftools 	--gzvcf $TMP'_tmp3.vcf.gz' 		\
			--out 	$TMP'_tmp5_qc'			\
			--recode-INFO-all --recode		\
			--exclude-positions $TMP'_tmp4_qc.exclude_pval0.01.list'
mv $TMP'_tmp5_qc.recode.vcf' $TMP'_tmp5_qc.vcf'
bgzip -fi $TMP'_tmp5_qc.vcf'
tabix -fp vcf $TMP'_tmp5_qc.vcf.gz'

#came after first QC filter 
#sed -i 's/\t\.:/\t\.\/\.:/g' $TMP'_tmp3.recode.vcf' # change genotypes from . to ./., dont think this is nessesary anymore with bcftools
