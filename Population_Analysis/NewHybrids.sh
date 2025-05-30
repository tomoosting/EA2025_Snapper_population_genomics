#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --partition=parallel
#SBATCH --time=10-0:00
#SBATCH --job-name=NewHybrids
#SBATCH -o /nfs/scratch/oostinto/stdout/NewHybrids.%j.out
#SBATCH -e /nfs/scratch/oostinto/stdout/NewHybrids.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tom.oosting@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3
EXT=$SET'_'$FILTER

###program
newhybrids=~/bin/newhybrids/newhybrids-no-gui-linux.exe

###set paths
READ_DIR=$SCRATCH/projects/$PROJECT/output/$SET/newhybrids/$EXT

###run parameters
Nburnin=100000
Nsweeps=1000000
Nsteps=5
Ntot=$(( $Nburnin + $Nsweeps - $Nsteps ))

##run newhybrids
cd $READ_DIR
$newhybrids --data-file $EXT.txt --gtyp-cat-file GtypeFreq.txt --burn-in $Nburnin --num-sweeps $Nsweeps --print-traces Pi $Nsteps --no-gui 1> trace.out

#genrate trace file
echo "Sample"        > sample_tmp.log
seq 0 $Nsteps $Ntot >> sample_tmp.log
cat trace.out | awk '/^PI_TRACE:/ {sub(/^PI_TRACE:/, ""); print}' > trace_tmp.log
paste sample_tmp.log trace_tmp.log > trace.log
rm sample_tmp.log trace_tmp.log

#Command Line Switches for Standard Analysis
#-d , --data-file        F               pathname to the data file
#-c , --gtyp-cat-file    F               path to file holding the genotype category probabilities
#-g , --gtyp-freq-probs  S R0 R1 R2      specify genotype frequency category S with probabilities R0 R1 and R2
#     --alle-prior-file  F               path to the file holding information about the priors for allele frequencies
#     --theta-prior      C               Set the type of prior used for the prior on theta
#     --pi-prior         C               Set the type of prior used for the prior on pi
#     --pi-prior-vec     R1...Rk         Parameters for pi prior if choosing Dirichlet or fixed with option pi-prior
#     --burn-in          J               run for J sweeps of burn-in
#     --num-sweeps       J               run for J sweeps AFTER burn-in
#-s , --seeds            J1 J2           seeds for the random number generator
#     --no-gui                           disable the GLUT/OpenGL MCMC visualizer
#Controlling Output
#     --print-traces     S J ...         Tell program to print trace of variable type S every J sweeps