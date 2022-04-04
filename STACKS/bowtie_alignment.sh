#!/bin/bash -l

#PBS -N <jobname>
#PBS -q darwin
#PBS -l nodes=1:ppn=10,naccesspolicy=singleuser
#PBS -l walltime=300:00:00
#PBS -m abe
#PBS -M <email>

#module version based on most up to date as of 25 Feb 2020
module --force purge
module load bioinfo
module load bowtie2/2.3.5.1
module load samtools/1.9

projhome=/scratch/snyder/c/cschraid/ypseq

proc=$projhome/working/processed_data/all ## clone_filtered, demultiplexed reads
al_out=$projhome/working/processed_data/refaligned ## output directory


# Ref genome index (using bowtie2-build)
bowtie_db=$projhome/ncbi/ncbi-genomes-2020-02-25/pflavescens

samples=$(ls $proc/*.1.fq.gz) ## make a list of samples to loop through

#rm $al_out/bowtie.log

echo 'Aligning samps to' $bowtie_db > $al_out/bowtie.log

for i in $samples; do
        samp=$(echo $i | cut -d / -f 8 | cut -f 1 -d .) #check this output before running full loop to make sure it is delimiting your data correctly
        echo 'Aligning' $samp

        echo $samp >> $al_out/bowtie.log

        ## conservatively, use --very-sensitive 
        bowtie2 --very-sensitive -p 32 -x $bowtie_db \
        -1 $proc/${samp}.1.fq.gz -2 $proc/${samp}.2.fq.gz \
        -S $al_out/${samp}.sam 2>> $al_out/bowtie.log
        samtools view -b $al_out/${samp}.sam | samtools sort --threads 32 \
        -o $al_out/${samp}.bam

done