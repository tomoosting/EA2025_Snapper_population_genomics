#!/bin/bash
#SBATCH -a 1-24
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=parallel
#SBATCH --time=5-0:00
#SBATCH --job-name=pixy
#SBATCH -o /nfs/scratch/oostinto/stdout/pixy_%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/pixy_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

export PATH=/nfs/home/oostinto/bin/htslib-1.18:$PATH
export PATH=/nfs/home/oostinto/bin/bcftools-1.18:$PATH
#module load bcftools/1.10.1
#module load htslib/1.10.1

###run input
PROJECT=$1
SET=$PROJECT'_'$2
LG=${SLURM_ARRAY_TASK_ID}

POP=$SCRATCH/projects/$PROJECT/resources/sample_info/$SET'_pixy.tsv'
BAM=$SCRATCH/projects/$PROJECT/resources/bam_lists/$SET'_bam.list'
REF=$SCRATCH/projects/$PROJECT/resources/reference_genomes/nuclear/Chrysophrys_auratus.v.1.0.all.assembly.units.fasta
REG=$( head -n $LG $REF.fai | tail -n 1 | cut -f 1 )

###set paths
VCF=$SCRATCH/projects/$PROJECT/data/snp/$SET/$SET
OUT=$SCRATCH/projects/$PROJECT/output/$SET/pixy
TMP=$SCRATCH/projects/$PROJECT/output/$SET/pixy/tmp

mkdir -p $TMP

##genotyping and filter
#echo "performing genotyping and filtering"
#bcftools mpileup 	-f $REF 							\
#					-b $BAM 							\
#					-r $REG								\
#					-a 'INFO/AD,FORMAT/AD,FORMAT/DP'	|
#bcftools call 		-m 									\
#					-Ou 								\
#					-f GQ 								|	 
#bcftools +fill-tags	-Ou 								\
#					-- -t AC,AF,AN,MAF,NS				|
#bcftools filter		-Ou									\
#					-S .								\
#					--exclude 'FMT/DP<3 | FMT/GQ<20'	|
#bcftools view		-Oz									\
#					-M2									\
#					--exclude 'STRLEN(REF)!=1 | STRLEN(ALT) >=2 | QUAL<600 | AVG(FMT/DP)<8 | AVG(FMT/DP)>25 ' 	\
#					-o $TMP/$REG'_'$SET'_raw.allsites.vcf.gz'

#create index
echo "creating index"
tabix -f $TMP/$REG'_'$SET'_raw.allsites.vcf.gz'

# run pixy
echo "running pixy"
source activate /nfs/scratch/oostinto/conda/pixy
pixy	--stats pi dxy fst								\
		--vcf $TMP/$REG'_'$SET'_raw.allsites.vcf.gz'	\
		--populations $POP								\
		--window_size 5000								\
		--n_cores 8										\
		--output_folder $OUT							\
		--output_prefix $REG'_'$SET	
