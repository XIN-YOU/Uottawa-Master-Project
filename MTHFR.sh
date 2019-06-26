#!/bin/bash
#SBATCH --job-name=MTHFD1L
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xyou016@uottawa.ca
#SBATCH --output=STD.out
#SBATCH --error=STD.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:00:00  #0 days 2 hours
#SBATCH --mem=15GB
# commands for your job go here
module load r

R CMD BATCH /global/home/hpc4429/XINYOU/gene1000/Rcode/MTHFR.R