#!/bin/bash
#SBATCH --job-name=v2tgatk
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/shared/Project_3_Resources/Group1/OandE/%x.out
#SBATCH --error=/shared/Project_3_Resources/Group1/OandE/%x.err

#activate Conda and change directory
source $HOME/.bash_profile
conda activate /shared/Project_3_Resources/Group1/shared_envs/gatk
cd /shared/Project_3_Resources/Group1/Data

vcf=LAB_NEN_ODN.clean_BI.ann.vcf.gz
ref=C_excelsa_V5.fa

gzvcfs=*.vcf.gz

#run gatk VarianttoTable

for file in $gzvcfs
do
     if [[ $file == LAB.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP --output nfLAB_raw.table
     elif [[ $file == NEN.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP --output nfNEN_raw.table
     elif [[ $file == ODN.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP --output nfODN_raw.table
     fi
done
