#!/bin/bash
#SBATCH --job-name=Scantools_g1
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/shared/Project_3_Resources/Group1/OandE/%x%j.out
#SBATCH --error=/shared/Project_3_Resources/Group1/OandE/%x%j.err

#activate Conda and change directory
source $HOME/.bash_profile
conda activate /shared/Project_3_Resources/Group1/shared_envs/gatk

vcf=/shared/Project_3_Resources/Group1/Data/LAB_NEN_ODN.clean_BI.ann.vcf.gz
ref=/shared/Project_3_Resources/Group1/Data/C_excelsa_V5.fa

#filter large VCF into 3 population VCFs (NEN, LAB and ODN)
gatk SelectVariants -R $ref -V $vcf --select "AF > 0.00" -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -O /shared/Project_3_Resources/Group1/ScanTools/LAB.vcf
gatk SelectVariants -R $ref -V $vcf --select "AF > 0.00" -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -O /shared/Project_3_Resources/Group1/ScanTools/NEN.vcf
gatk SelectVariants -R $ref -V $vcf --select "AF > 0.00" -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O /shared/Project_3_Resources/Group1/ScanTools/ODN.vcf

# variants to table

gatk VariantsToTable -V /shared/Project_3_Resources/Group1/ScanTools/LAB.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O /shared/Project_3_Resources/Group1/ScanTools/LAB.table
gatk VariantsToTable -V /shared/Project_3_Resources/Group1/ScanTools/NEN.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O /shared/Project_3_Resources/Group1/ScanTools/NEN.table
gatk VariantsToTable -V /shared/Project_3_Resources/Group1/ScanTools/ODN.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O /shared/Project_3_Resources/Group1/ScanTools/ODN.table

conda deactivate
conda activate /shared/Project_3_Resources/Group1/shared_envs/scantools_env

# Into recode012.py
# LAB
python3 /shared/Project_3_Resources/Group1/ScanTools/recode012.py -i /shared/Project_3_Resources/Group1/ScanTools/LAB.table -o /shared/Project_3_Resources/Group1/ScanTools/ -pop LAB
# NEN
python3 /shared/Project_3_Resources/Group1/ScanTools/recode012.py -i /shared/Project_3_Resources/Group1/ScanTools/NEN.table -o /shared/Project_3_Resources/Group1/ScanTools/ -pop NEN
# ODN
python3 /shared/Project_3_Resources/Group1/ScanTools/recode012.py -i /shared/Project_3_Resources/Group1/ScanTools/ODN.table -o /shared/Project_3_Resources/Group1/ScanTools/ -pop ODN

# merge for pairwise comparison
# LAB vs NEN
python3 /shared/Project_3_Resources/Group1/ScanTools/Sian_sort_for_ScanTools.py '/shared/Project_3_Resources/Group1/ScanTools/LAB.table.recode.txt /shared/Project_3_Resources/Group1/ScanTools/NEN.table.recode.txt ' /shared/Project_3_Resources/Group1/ScanTools/ LAB_Vs_NEN
# LAB vs ODN
python3 /shared/Project_3_Resources/Group1/ScanTools/Sian_sort_for_ScanTools.py '/shared/Project_3_Resources/Group1/ScanTools/LAB.table.recode.txt /shared/Project_3_Resources/Group1/ScanTools/NEN.table.recode.txt ' /shared/Project_3_Resources/Group1/ScanTools/ LAB_Vs_ODN

# Into bmp.py
# LAB vs NEN
python3 /shared/Project_3_Resources/Group1/ScanTools/sian_bpm.py -i /shared/Project_3_Resources/Group1/ScanTools/LAB_Vs_NEN.concat.txt -o /shared/Project_3_Resources/Group1/ScanTools/ -prefix LAB_Vs_NEN_15 -ws 1000 -ms 15 -np 2
# LAB vs ODN
python3 /shared/Project_3_Resources/Group1/ScanTools/sian_bpm.py -i /shared/Project_3_Resources/Group1/ScanTools/LAB_Vs_ODN.concat.txt -o /shared/Project_3_Resources/Group1/ScanTools/ -prefix LAB_Vs_ODN_15 -ws 1000 -ms 15 -np 2
