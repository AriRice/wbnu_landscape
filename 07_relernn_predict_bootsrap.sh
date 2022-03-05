#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=relernn2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --partition matador
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --gpus-per-node=2

source activate recomb

# used the Manthey et al. 2021 Certhia mutation rate

SEED="42"
MU="2.506e-9"
DIR="./relernn_output/"
VCF="./nuthatch_relernn.vcf"
GENOME="./wbnu.genome.bed"
GENTIME="2"

# Predict recombination rates
ReLERNN_PREDICT \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED} \
    --unphased
    
# bootstrap and correct
ReLERNN_BSCORRECT \
    --projectDir ${DIR} \
    --seed ${SEED} \
    --nCPU 40 \
    --nSlice 100 \
    --nReps 100
    

    
