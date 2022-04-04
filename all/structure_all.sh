#!/bin/bash


#SBATCH --job-name=structure_all_2.1
#SBATCH -A beagle
#SBATCH -t 336:00:00
#SBATCH -N 1-1
#SBATCH -n 1
#SBATCH --mail-type=ALL
#SBATCH --<email>
#SBATCH -o structure.o
#SBATCH -e structure.e

cd $SLURM_SUBMIT_DIR

module load anaconda
source activate myenv
structure -m mainparams -e extraparams.popalpha -o

