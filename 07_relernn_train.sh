#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=relernn
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --partition matador
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=9G
#SBATCH --gpus-per-node=1

source activate recomb

# used the Manthey et al. 2021 Certhia mutation rate

SEED="42"
DIR="./relernn_output/"

# Simulate data
ReLERNN_TRAIN \
    --projectDir ${DIR} \
    --seed ${SEED} \

    
    
