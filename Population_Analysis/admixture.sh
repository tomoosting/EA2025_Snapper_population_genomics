#!/bin/bash
#SBATCH --cpus-per-task=14
#SBATCH --mem=10G
#SBATCH --partition=parallel
#SBATCH --time=10-0:00
#SBATCH --job-name=admixture
#SBATCH -o /nfs/scratch/oostinto/stdout/admixture_%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/admixture_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

#run ADMIXTURE on a dataset containing neutral indepedantly segregating SNPs!!
#ADMIXTURE runs on a bed file, first convert vcf to bed
###run input
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3
K=$4

#documentation
#https://github.com/laurabenestan/Admixture

#load modules
module load ADMIXTURE/1.3.0
module load sf/0.9-5-R-4.0.0-Python-3.8.2


#programs
gds2plink=$SCRATCH/scripts/population_genomics/4_variant_filtering/R/gds2plink.R
vcf2R=$SCRATCH/scripts/population_genomics/4_variant_filtering/R//vcf2Rinput.R
evalAdmix=/nfs/home/oostinto/bin/evalAdmix/evalAdmix
export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH

###set paths
SNP=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
OUT=$SCRATCH/projects/$PROJECT/output/$SET/admixture/$FILTER
EXT=$SET'_'$FILTER
mkdir -p $OUT

####resrouces
#pop_file=$SCRATCH/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'
##convert to plink					
#Rscript $vcf2R 		--gzvcf $SNP'_'$FILTER.vcf.gz 	\
#					--snprelate_out $SNP'_'$FILTER	
#
#
#Rscript $gds2plink 	--gds_file	$SNP'_'$FILTER'.gds'	\
#					--out 		$SNP'_'$FILTER			\
#					--pop_file	$pop_file

################################## ADMIXTURE ########################################
for Ki in $( seq $K ) 
do
#admixture
	admixture --cv=10 -B10000 -j10 $SNP'_'$FILTER'.bed' $Ki | tee log${Ki}.out
#evaladmix
	bgzip -c $EXT.$Ki.P > $EXT.$Ki.P.gz
	$evalAdmix  -plink  $SNP'_'$FILTER		\
				-fname	$EXT.$Ki.P.gz 		\
				-qname  $EXT.$Ki.Q			\
				-o 		$EXT.$Ki.evaladmix	\
				-autosomeMax 25				\
				-minMaf 0.05				\
				-P 10
	rm $EXT.$Ki.P.gz
done

#move file to output folder
	mv log*.out $OUT/
	mv $SET* 	$OUT/
	grep -h CV  $OUT/log*.out > $OUT/cross_validation.txt

