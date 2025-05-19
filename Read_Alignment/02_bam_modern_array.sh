#!/bin/bash
#SBATCH -a 1-382
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --time=1-0:00
#SBATCH --partition=parallel
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --job-name=bam_pipeline_array
#SBATCH -o /nfs/scratch/oostinto/stdout/bam_pipeline_array_%A_%a.out
#SBATCH -e /nfs/scratch/oostinto/stdout/bam_pipeline_array_%A_%a.err
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

module load paleomix/1.2.13.3
module load bwa-kit/0.7.15
module load AdapterRemoval/2.2.3
module load picard/2.18.20
module load java/jdk/1.8.0_121
module load gatk/4.0.8.1
module load picard/2.18.20
module load samtools/1.9
module load R/3.5.1
module load R/CRAN/3.5
module load mapDamage/2.0.9
module load fastqc/0.11.7

PROJECT=$1
LIB=$2
N=${SLURM_ARRAY_TASK_ID}

lib_file=/nfs/scratch/oostinto/projects/$PROJECT/resources/paleomix_lists/$PROJECT'_'$LIB.txt
sample=$( head -n $N $lib_file | tail -n 1 | cut -f 1 )
ext=$( head -n $N $lib_file | tail -n 1 | cut -f 2 )

raw_read_dir=/nfs/scratch/oostinto/projects/$PROJECT/raw_data/illumina
read_dir=/nfs/scratch/oostinto/projects/$PROJECT/data/paleomix/$sample
fastqc_out=/nfs/scratch/oostinto/projects/$PROJECT/data/fastqc

#pre fastqc
echo "perform fastqc on raw data"
fastqc $raw_read_dir/$sample'_'$ext* -t 12 --noextract -o $fastqc_out

#bam pipeline
echo "running PALEOMIX"
paleomix bam_pipeline run $read_dir/*.yaml

#post fastqc 
echo "perform fastqc on trimmed data"
fastqc $read_dir/*_m1.0/reads/*/*/Lane_1/reads.pair*.truncated.gz -t 12 --noextract -o $fastqc_out/trimmed_reads
mv $fastqc_out/trimmed_reads/reads.pair1.truncated_fastqc.html $fastqc_out/trimmed_reads/$sample'_'$ext.pair1.truncated_fastqc.html
mv $fastqc_out/trimmed_reads/reads.pair2.truncated_fastqc.html $fastqc_out/trimmed_reads/$sample'_'$ext.pair2.truncated_fastqc.html
mv $fastqc_out/trimmed_reads/reads.pair1.truncated_fastqc.zip  $fastqc_out/trimmed_reads/$sample'_'$ext.pair1.truncated_fastqc.zip
mv $fastqc_out/trimmed_reads/reads.pair2.truncated_fastqc.zip  $fastqc_out/trimmed_reads/$sample'_'$ext.pair2.truncated_fastqc.zip

#removing soft clipped reads
echo "removing soft clipped reads from bam"
in_bam=$read_dir/$sample'_m1.0.nuclear.realigned.bam'
list_file=$read_dir/$sample'_exclude_clipped_reads_list.txt'
out_bam=$read_dir/$sample'_m1.0.nuclear.realigned.no_clipped_reads.bam'
samtools view $in_bam |awk '$6 ~ /H|S/ {print $1}' |sort -u > $list_file
samtools view -h $in_bam | fgrep -wvf $list_file | samtools view -b - -o $out_bam
samtools index $out_bam

#clean up
echo "deleting unwanted files"
rm -r $read_dir/${sample}_m1.0*/
rm $read_dir/*.coverage
rm $read_dir/*.depths
rm $read_dir/*_exclude_clipped_reads_list.txt
rm $read_dir/*.mtgenome.bai
rm $read_dir/*.mtgenome.bam
rm $read_dir/*.nuclear.bam
rm $read_dir/*.nuclear.bai
rm $read_dir/*.nuclear.realigned.bai
rm $read_dir/*.nuclear.realigned.bam





