#!/bin/bash

#PBS -N process_radtags123
#PBS -q darwin
#PBS -l nodes=1:ppn=10,naccesspolicy=singleuser
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M [email]

cd $PBS_O_WORKDIR


ml purge
ml bioinfo stacks

#process for processing 1 plate (10 tot) of raw illumina reads

process_radtags -1 ./AA19_123/AA19_123_CKDL200146451-1a_HTC5MDSXX_L1_1.fq.gz -2 ./AA19_123/AA19_123_CKDL200146451-1a_H$
                 -b ./AA19_123/19-123.txt -e sbfI \
                                -i gzfastq -c -q -r \
                                --filter_illumina \
                                --bestrad \
                                 -t 140 \
                                -o ./processed_data/P1 \