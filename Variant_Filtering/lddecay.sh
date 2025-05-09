#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=LDdecay
#SBATCH -o /nfs/scratch/oostinto/stdout/LDdecay.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/LDdecay.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load htslib/1.9
module load vcftools/0.1.16
module load R/4.0.2
module load plink/1.90

plink=~/bin/plink_v1.90/plink

#Rscript paths
R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R
LDdecay=$R/LDdecay_shell.R


###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
OUT=$SCRATCH/projects/$PROJECT/output/$SET/lddecay

### create directories
mkdir -p $OUT

# calc ld with plink
echo "starting plink"
$plink
$plink --bfile $VCF'_qc' 	--allow-extra-chr 			\
	   --double-id 			--set-missing-var-ids @:# 	\
	   --maf 0.01 			--geno 0.1 					\
	   --mind 0.5 			--chr 1 					\
	   --thin 0.1 -r2 gz 	--ld-window 100 			\
	   --ld-window-kb 1000	--ld-window-r2 0 			\
	   --out $OUT/$SET'_qc_plink'

bgzip -cd $OUT/$SET'_qc_plink.ld.gz' | sed 's/[[:blank:]]/\t/g' | sed 's/^\t//g' | sed 's/\t$//g' > $OUT/$SET'_qc_plink.ld'