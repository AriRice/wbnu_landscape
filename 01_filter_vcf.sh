#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-39

# define input files from helper file during genotyping
input_vcf=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/04_nuthatch_recomb

# pull out header and add to filtered vcf file
grep "#" ${workdir}/03_vcf/${input_vcf}.g.vcf > ${workdir}/03_vcf/${input_vcf}.filtered.vcf

# filter our rows that have low quality filters, genotyped sites with quality less than 20, and null alleles (* in col 4)
grep -v "#" ${workdir}/03_vcf/${input_vcf}.g.vcf | grep -v "LowQual" | awk '$6 >= 20 || $6 ~ /^\./' | awk '$5 !~ /*/' >> ${workdir}/03_vcf/${input_vcf}.filtered.vcf

# run vcftools to get variant and invariant sites found in at least 0.5 (14 of 27) individuals
vcftools --vcf ${workdir}/03_vcf/${input_vcf}.filtered.vcf --max-missing 0.5 --minDP 6 --max-meanDP 40 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/04_stats/${input_vcf}

# bgzip and tabix index files that will be subdivided into windows
bgzip ${workdir}/04_stats/${input_vcf}.recode.vcf

tabix -p vcf ${workdir}/04_stats/${input_vcf}.recode.vcf.gz

# run vcftools with SNP output for relernn, minimum 0.72 (20 of 27 individuals) retained
vcftools --vcf ${workdir}/03_vcf/${input_vcf}.filtered.vcf --max-missing 0.72 --minDP 6 --max-meanDP 40 --min-alleles 2 --max-alleles 2 --mac 1 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/05_relernn/${input_vcf}
