#!/bin/bash


#SBATCH --job-name=structure_mb_2.5
#SBATCH -A beagle
#SBATCH -t 140:00:00
#SBATCH -N 1-1
#SBATCH -n 10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<email>
#SBATCH -o structure.o
#SBATCH -e structure.e

cd $SLURM_SUBMIT_DIR

module load anaconda
source activate myenv
structure -m mainparams -e ../extraparams -o 2.5