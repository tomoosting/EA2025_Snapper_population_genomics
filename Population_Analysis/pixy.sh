#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=regular
#SBATCH --time=1-0:00
#SBATCH --job-name=pixy

###run input
PROJECT=$1 #snapper
SET=$PROJECT'_'$2 #norm 
WINDOW=$3 #10000

POP=$SCRATCH/projects/$PROJECT/datasets/$SET/sample_data/$SET'_pixy.tsv'

###set paths
VCF_DIR=$SCRATCH/projects/$PROJECT/datasets/$SET/genotype_data/$SET'_allsites'
OUT_DIR=$SCRATCH/projects/$PROJECT/datasets/$SET/population_analyses/selection/pixy

# run pixy
echo "running pixy"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pixy
pixy	--stats pi dxy fst watterson_theta tajima_d \
		--vcf $VCF_DIR/$SET'_allsites.vcf.gz'       \
		--populations $POP							\
		--window_size $WINDOW						\
		--n_cores 8									\
		--output_folder $OUT_DIR				    \
		--output_prefix $SET'_window'$WINDOW	
conda deactivate

