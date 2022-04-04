#!/bin/bash

#SBATCH -A beagle
#SBATCH --job-name=gstacksv2.5
#SBATCH -o gstacks.out
#SBATCH -e gstacks.err
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 300:00:00
#SBATCH -exclusive
#SBATCH --mail-user= <email>
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

cd $SLURM_SUBMIT_DIR

module load stacks/2.5

gstacks -I ./refalignedbam -M ./popmap.txt -t 1 -O ./gstacks \
        --rm-unpaired-reads \
        --details \

done