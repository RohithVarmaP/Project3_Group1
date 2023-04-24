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

#Prepare reference dictionary (reference index must be prepared in advance using samtools faidx)
gatk CreateSequenceDictionary -R $ref

#filter large VCF into 3 population VCFs (NEN, LAB and ODN)
gatk SelectVariants -V $vcf --select "AF > 0.00" -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -O LAB_AF_filter.$
gatk SelectVariants -V $vcf --select "AF > 0.00" -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -O NEN_AF_filter.vcf
gatk SelectVariants -V $vcf --select "AF > 0.00" -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O ODN_AF_filter.vcf

#set new variables and gzip vcfs
vcfs=*.vcf
conda deactivate

source $HOME/.bash_profile
conda activate /shared/Project_3_Resources/Group1/shared_envs/faststructure

for file in $vcfs
do
     bgzip -c $file > "$file".gz
     tabix -p vcf "$file".gz
     rm $file
done

conda deactivate

gzvcfs=*.vcf.gz

#run gatk VarianttoTable
source $HOME/.bash_profile
conda activate /shared/Project_3_Resources/Group1/shared_envs/gatk

for file in $gzvcfs
do
     if [[ $file == LAB_AF_filter.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP -raw --output LAB_raw.table
     elif [[ $file == NEN_AF_filter.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP -raw --output NEN_raw.table
     elif [[ $file == ODN_AF_filter.vcf.gz ]]; then
          gatk VariantsToTable -V $file -R $ref -F CHROM -F POS -F AC -F AN -F DP -raw --output ODN_raw.table
     fi
done
