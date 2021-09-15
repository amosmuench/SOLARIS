#!/bin/bash
#SBATCH --job-name=Snake
#SBATCH --mail-user=amos.muench@fu-berlin.de  
#SBATCH --mail-type=end
#SBATCH --ntasks=1                           
#SBATCH --mem-per-cpu=4096                  
#SBATCH --time=08:00:00                         
#SBATCH --qos=standard 
#SBATCH --nodes=1



module load Miniconda3/4.7.10
source activate /home/amosmuench/snakemake
source activate /home/amosmuench/snakemake/envs/snakemake_pulp/

cd /home/amosmuench/SOLARIS/
snakemake -j 4 --use-conda 