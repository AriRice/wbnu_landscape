#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=relernn
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
GENTIME="1"

# Simulate data
ReLERNN_SIMULATE \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --seed ${SEED} \
    --assumedGenTime ${GENTIME} \
    --unphased

# Training
ReLERNN_TRAIN \
    --projectDir ${DIR} \
    --seed ${SEED} \
    --nCPU 40
    
    
