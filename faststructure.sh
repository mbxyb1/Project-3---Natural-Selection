#!/bin/bash
#SBATCH --job-name=faststructure
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/shared/Project_3_Resources/Group4/OandE/%x.out
#SBATCH --error=/shared/Project_3_Resources/Group4/OandE/%x.err

#activate Conda and change directory
source $HOME/.bash_profile
conda activate /shared/Project_3_Resources/Group4/shared_envs/faststructure
cd /shared/Project_3_Resources/Group4/faststructure

bgzip -c reheadered_4dg_dips.clean_BI.ann.vcf > reheadered_4dg_dips.clean_BI.ann.vcf.gz
tabix -p vcf reheadered_4dg_dips.clean_BI.ann.vcf.gz

#run Plink to convert vcf to bed format
plink --allow-extra-chr --vcf reheadered_4dg_dips.clean_BI.ann.vcf.gz --make-bed --out reheadered

for i in {1..15}
do
     #Run faststructure using K values (number of populations) of 1-15 with 10-fold cross validation
     structure.py -K $i --input=reheadered --cv=10 --output=faststructure_outputK"$i" --format=bed

     #Generate plots of the faststructure output to help better infer the best K value to display population structure
     chooseK.py --input=faststructure_outputK"$i" > chooseK_outK"$i".txt
done

