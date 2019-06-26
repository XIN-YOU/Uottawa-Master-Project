#!/bin/bash
#SBATCH --job-name=pricetable
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xinyou.work@gmail.com
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-02:00:00  #0 days 2 hours
#SBATCH --mem=4GB
# commands for your job go here
module load r

R CMD BATCH /root/PROJECT/ab13.R

