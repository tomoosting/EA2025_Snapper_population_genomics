#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=10-0:00
#SBATCH --partition=regular
#SBATCH --mail-type=END
#SBATCH --job-name=selscan
#SBATCH -o /scratch/p268537/stdout/selscan_%j.out
#SBATCH -e /scratch/p268537/stdout/selscan_%j.err

# =============================================================================
# run_xpnsl.sh
# Runs selscan XP-nSL on two populations from VCF input (unphased or phased),
# then normalises scores with norm. Designed for non-human animals without a
# genetic map (physical positions only).
#
# Usage example: 
#   bash run_xpnsl.sh snapper_norm east west
#
# Dependencies:
#   selscan >= 3.0.0   (https://github.com/szpiech/selscan)
#   bcftools           (for VCF subsetting)
#export program path
export PATH=/home2/p268537/bin/selscan/src:$PATH 
#load modules
module load BCFtools/1.22-GCC-14.2.0
module load HTSlib/1.22.1-GCC-14.2.0
# =============================================================================

set -euo pipefail

# =============================================================================
# USER-DEFINED INPUTS — edit these before running
# =============================================================================

###run input
PROJECT=$1 #snapper
DATA_SET=$2 #norm
POP_1=$3 #NE
POP_2=$4 #SW
SNP_SET=$5 #qc

SET=$PROJECT'_'$DATA_SET
SNP=$PROJECT'_'$DATA_SET'_'$SNP_SET
VCF=/scratch/p268537/projects/$PROJECT/datasets/$SET/genotype_data/$SNP/$SNP.vcf.gz

#OUTPUT
OUT_DIR=/scratch/p268537/projects/$PROJECT/datasets/$SET/population_analyses/selection/selscan
OUT_EXT=$SET'_'$POP_1'_'$POP_2
mkdir -p $OUT_DIR

#sample lists for subsetting VCFs
INDS_POP_1=/scratch/p268537/projects/$PROJECT/datasets/$SET'_'$POP_1/sample_data/$SET'_'$POP_1.tsv
INDS_POP_2=/scratch/p268537/projects/$PROJECT/datasets/$SET'_'$POP_2/sample_data/$SET'_'$POP_2.tsv

# Chromosome/scaffold to run — loop over all chromosomes below,
# or set a single value for testing: e.g. CHROMS="chr1"
# List all chromosomes/scaffolds separated by spaces.
CHROMS="LG1 LG2 LG3 LG4 LG5 LG6 LG7 LG8 LG9 LG10 LG11 LG12 LG13 LG14 LG15 LG16 LG17 LG18 LG19 LG20 LG21 LG22 LG23 LG24 LG25"

# =============================================================================
# PARAMETER SETTINGS — justified in detail below the script
# =============================================================================

# --max-gap: maximum allowed physical gap (bp) between consecutive SNPs.
# Pairs of SNPs separated by more than this are treated as non-informative
# for haplotype extension. For non-human genomes without a genetic map,
# this guards against spuriously long iHH integration across assembly gaps,
# centromeres, or poorly covered regions.
# Recommended: 200000 (200 kb). Adjust downward for fragmented assemblies.
MAX_GAP=200000

# --max-extend: maximum physical distance (bp) over which haplotype
# extension is allowed from the focal SNP. Caps integration to a
# biologically meaningful window and prevents long scaffolds with sparse
# SNPs from producing inflated iHH values.
# Recommended: 1000000 (1 Mb). For organisms with short LD blocks
# (e.g. highly outcrossing species), reduce to 500000.
MAX_EXTEND=500000

# --min-snps: minimum number of SNPs that must remain within the
# extension window for a site to be scored. Sites with too few
# flanking SNPs produce unreliable scores.
# Recommended: 10. Reduce to 5 only for very sparse SNP panels.
#MIN_SNPS=10

# --maf: minimum minor allele frequency filter applied internally
# by selscan. Sites below this threshold are skipped. This is a
# safety net — your VCF should already be MAF-filtered, but this
# catches any residual near-monomorphic sites that distort iHH.
# Recommended: 0.05
MAF=0.05

# --norm bins: number of frequency bins used during normalisation.
# norm standardises XP-nSL scores within allele frequency bins.
# 100 bins is standard and appropriate for most dataset sizes.
# Reduce to 50 if your total SNP count after filtering is < 50,000.
NORM_BINS=100

# =============================================================================
# STEP 1 — Run selscan XP-nSL per chromosome
# =============================================================================

for CHR in ${CHROMS}; do
    echo "filtering missing for missing data"
    bcftools filter -e 'F_MISSING > 0' --threads $SLURM_CPUS_PER_TASK -Oz -o $OUT_DIR/${CHR}_${SET}_no_missing.vcf.gz $VCF
    bcftools index --threads $SLURM_CPUS_PER_TASK $OUT_DIR/${CHR}_${SET}_no_missing.vcf.gz

    echo "subsetting VCFs by CHROM for selscan..."
    echo "  Pop 1: ${POP_1}"
    bcftools view -r $CHR -S $INDS_POP_1 --threads $SLURM_CPUS_PER_TASK -Oz -o $OUT_DIR/$CHR'_'$SET'_'$POP_1.vcf.gz $OUT_DIR/${CHR}_${SET}_no_missing.vcf.gz
    echo "  Pop 2: ${POP_2}"
    bcftools view -r $CHR -S $INDS_POP_2 --threads $SLURM_CPUS_PER_TASK -Oz -o $OUT_DIR/$CHR'_'$SET'_'$POP_2.vcf.gz $OUT_DIR/${CHR}_${SET}_no_missing.vcf.gz
    
    echo "removing temporary VCF file"
    rm $OUT_DIR/${CHR}_${SET}_no_missing.vcf.gz
    
    echo "cross population analysis for ${CHR}..."
    selscan                                               \
        --vcf        $OUT_DIR/$CHR'_'$SET'_'$POP_1.vcf.gz \
        --vcf-ref    $OUT_DIR/$CHR'_'$SET'_'$POP_2.vcf.gz \
        --out        $OUT_DIR/$CHR'_'$OUT_EXT             \
        --threads    $SLURM_CPUS_PER_TASK                 \
        --max-gap    $MAX_GAP                             \
        --max-extend-nsl $MAX_EXTEND                      \
        --maf        $MAF                                 \
        --unphased                                        \
        --xpnsl

   echo "cleanup temporary files..."
   rm $OUT_DIR/$CHR'_'$SET'_'$POP_1.vcf.gz*
   rm $OUT_DIR/$CHR'_'$SET'_'$POP_2.vcf.gz*
done

# =============================================================================
# STEP 2 — Normalise with norm
# =============================================================================

WINDOW_SIZE=10000

echo "Normalising XP-nSL scores..."
selscan norm                                \
    --xpnsl                                 \
    --files $OUT_DIR/*_"$OUT_EXT".xpnsl.out \
    --qbins  10                             \
    --winsize $WINDOW_SIZE                  \
    --bp-win

# =============================================================================
# STEP 3 — Concatenate normalised files per-chromosome output files
# =============================================================================

echo "Concatenating chromosome outputs..."
OUT_COMBINED=$OUT_DIR/$OUT_EXT.xpnsl.out.norm
> $OUT_COMBINED
> ${OUT_COMBINED}.10kb.windows

XPNSL_COMBINED=$OUT_DIR/$OUT_EXT.xpnsl.out.norm
HEADER_WRITTEN=0
for CHR in ${CHROMS}; do
    CHROM_FILE=$OUT_DIR/$CHR'_'$OUT_EXT.xpnsl.out.norm
    if [ -f "${CHROM_FILE}" ]; then
        if [ "${HEADER_WRITTEN}" -eq 0 ]; then
            cat "${CHROM_FILE}" >> "${OUT_COMBINED}"
            HEADER_WRITTEN=1
        else
            # Skip header line on subsequent files
            tail -n +2 "${CHROM_FILE}" >> "${OUT_COMBINED}"
        fi
    else
        echo "  WARNING: expected output not found: ${CHROM_FILE}" >&2
    fi
done

XPNSL_NORM_COMBINED=$OUT_DIR/$OUT_EXT.xpnsl.out.norm.10kb.windows
HEADER_WRITTEN=0
for CHR in ${CHROMS}; do
    CHROM_FILE_WINDOWS=$OUT_DIR/$CHR'_'$OUT_EXT.xpnsl.out.norm.5kb.windows
    if [ -f "${CHROM_FILE_WINDOWS}" ]; then
        if [ "${HEADER_WRITTEN}" -eq 0 ]; then
            cat "${CHROM_FILE_WINDOWS}" >> "${OUT_COMBINED}.5kb.windows"
            HEADER_WRITTEN=1
        else
            # Skip header line on subsequent files
            tail -n +2 "${CHROM_FILE_WINDOWS}" >> "${OUT_COMBINED}.5kb.windows"
        fi
    else
        echo "  WARNING: expected output not found: ${CHROM_FILE_WINDOWS}" >&2
    fi
done
