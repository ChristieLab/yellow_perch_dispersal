#!/bin/bash

#SBATCH -A darwin
#SBATCH --job-name=populations
#SBATCH -o populaitions.out
#SBATCH -e populations.err
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH -exclusive
#SBATCH --mail-user= <email>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd $SLURM_SUBMIT_DIR

module load stacks/2.5

#single pop preliminary run for vcf output. no pop delim, as dataset is too large. 
populations -P ./gstacks 
	--vcf \
	-r 0.30 \ 
	-O populations \