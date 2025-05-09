#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=parallel
#SBATCH --time=1-0:00
#SBATCH --job-name=outlier_selection
#SBATCH -o /nfs/scratch/oostinto/stdout/outlier_selection.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/outlier_selection.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###!!!###
# you will need some prior information by runnign these analyses locally to set certain parameters (e.g. K)
###!!!###

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load bcftools/1.10.1
#module load R/4.0.2
module load plink/1.90

Rscript=/nfs/home/oostinto/bin/R-4.3.1/bin/Rscript


#Rscript paths
R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R
OutFlank=$R/OutFlank_bash.R
pcadaptR=$R/pcadapt_shell.R
vcf2R=$R/vcf2Rinput.R
pops=$R/get_pops.R

###resrouces
pop_file=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
TMP=$SCRATCH/projects/$PROJECT/data/snp/$SET/tmp
OUTFLANK=$SCRATCH/projects/$PROJECT/output/$SET/outlier_analyses/outflank
PCADAPT=$SCRATCH/projects/$PROJECT/output/$SET/outlier_analyses/pcadapt

### create directories
#mkdir $TMP
mkdir -p $PCADAPT
mkdir -p $OUTFLANK
##################################### filter outlier loci  #############################
QVAL=0.05
K=1
$Rscript $pcadaptR --plink   $VCF'_qc'     	\
				   --out     $PCADAPT/$SET 	\
				   --K0 	 5      		\
				   --K 	     $K             \
				   --maf     0.05          	\
				   --q 	     $QVAL          \
				   --slw     50000         	\
				   --minNsnp 5             	\
				   --mode	 full

$Rscript $pops --pop_file $pop_file --set_name $2 --out $TMP/$SET'_pops'
$Rscript $OutFlank 	--vcf_file 		 $VCF'_qc.vcf.gz'  		\
					--gds_file 		 $VCF'_qc.gds'			\
					--pop_file 		 $TMP/$SET'_pops.tsv'	\
					--fst_file		 $OUTFLANK/$SET			\
					--out	  		 $OUTFLANK/$SET		 	\
					--plots_out		 $OUTFLANK/$SET			\
					--qval 			 $QVAL					\
					--filter_method  "FST"					\
					--slw			 50000					\
					--min_snps		 2						\
					--by_chr		 TRUE					\
					--compress		 TRUE					\
					--fst_plot		 TRUE					\
					--qval_plot		 TRUE					\
					--PCA_plot 		 TRUE					\
					--Manhattan_plot TRUE					

#filter vcf file for outliers
if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_independent_outliers.tsv' ]
then
	#extract SNP info from outliers
	tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_independent_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_independent_outlier_LOC.tsv'
	#outputting outlier vcf file containing all snps
	vcftools 	--gzvcf $VCF'_qc.vcf.gz' 						\
				--positions $PCADAPT/$SET'_independent_outlier_LOC.tsv'	\
				--out $VCF'_outliers'							\
				--recode-INFO-all 								\
				--recode
	mv $VCF'_outliers.recode.vcf' $VCF'_outliers.vcf'
	bgzip -i $VCF'_outliers.vcf' 
	
	$Rscript $vcf2R 	--gzvcf $VCF'_outliers.vcf.gz' 	\
					--snprelate_out $VCF'_outliers.'	
	
fi			

