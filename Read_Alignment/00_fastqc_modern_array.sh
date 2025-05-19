#!/bin/bash
#SBATCH -a 2-84
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-1:00
#SBATCH --partition=quicktest
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=fastqc_array
#SBATCH -o /nfs/scratch/oostinto/stdout/fastqc_array_%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/fastqc_array_%A_%a.err
#SBATCH --mail-user=tom.oosting@vuw.ac.nz
 
module load fastqc/0.11.7

PROJECT=$1
LIB=$2
N=${SLURM_ARRAY_TASK_ID}

lib_file=/nfs/scratch/oostinto/projects/$PROJECT/resources/paleomix_lists/$PROJECT'_'$LIB.txt
dos2unix $lib_file

sample=$( head -n $N $lib_file | tail -n 1 | cut -f 1 )
ext=$( head -n $N $lib_file | tail -n 1 | cut -f 2 )

raw_read_dir=/nfs/scratch/oostinto/projects/$PROJECT/raw_data/illumina
fastqc_out=/nfs/scratch/oostinto/projects/$PROJECT/data/fastqc

fastqc $raw_read_dir/$sample'_'$ext* -t 12 --noextract -o $fastqc_out
