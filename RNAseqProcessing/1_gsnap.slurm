#!/bin/bash

#SBATCH -t 168:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 150   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "salt_gsnap"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-salt_gsnap.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# load modules
module load gmap-gsnap/2019-05-12-zjqshxf
module load parallel/20170322-nct4iah

# build genome
gmap_build -d D5 -D db Dgenome2_13.fasta

# put the snp files in the format it likes better 
cut -f1,2 D13.snp4.0 | sed 's/\t/_/g' | sed 's/Chr/>Chr/g' | paste - D13.snp4.0 | sed 's/\t/ /1' | sed 's/\t/:/1' | sed 's/\t/ /1' | sed 's/\t//g' > snp4.0
cut -f1,2 D13.snp4.1 | sed 's/\t/_/g' | sed 's/Chr/>Chr/g' | paste - D13.snp4.1 | sed 's/\t/ /1' | sed 's/\t/:/1' | sed 's/\t/ /1' | sed 's/\t//g' > snp4.1
cut -f1,2 D13.snp4.2 | sed 's/\t/_/g' | sed 's/Chr/>Chr/g' | paste - D13.snp4.2 | sed 's/\t/ /1' | sed 's/\t/:/1' | sed 's/\t/ /1' | sed 's/\t//g' > snp4.2
cut -f1,2 D13.snp4.4 | sed 's/\t/_/g' | sed 's/Chr/>Chr/g' | paste - D13.snp4.4 | sed 's/\t/ /1' | sed 's/\t/:/1' | sed 's/\t/ /1' | sed 's/\t//g' > snp4.4

# build SNP indices
cat snp4.0 | iit_store -o snp4.0
cat snp4.1 | iit_store -o snp4.1
cat snp4.2 | iit_store -o snp4.2
cat snp4.4 | iit_store -o snp4.4

# cp *.iit ./db/D5/D5.maps

snpindex -D db -d D5 -v snp4.0 snp4.0.iit
snpindex -D db -d D5 -v snp4.1 snp4.1.iit
snpindex -D db -d D5 -v snp4.2 snp4.2.iit
snpindex -D db -d D5 -v snp4.4 snp4.4.iit
 
ls [AD][25]*gz | sed 's/_[12]\..*//g'| sort | uniq > diploid.list
ls [AD][25]*web*gz | sed 's/_[12]\..*//g'| sort | uniq > diploid.web.list
ls AD4*gz | sed 's/_[12]\..*//g'| sort | uniq > AD4.list
ls PS*gz | sed 's/_[12]\..*//g'| sort | uniq > AD2.list
ls PS*web*gz | sed 's/_[12]\..*//g'| sort | uniq > AD2.web.list
ls TM*gz | sed 's/_[12]\..*//g'| sort | uniq > AD1.list


# run gsnap
parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.0 -A sam {}_1.fq.gz {}_2.fq.gz > {}.sam 2>> {}.gsnap.log' :::: diploid.list
parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.0 -A sam {}_1.web.fq.gz {}_2.web.fq.gz > {}.web.sam 2>> {}.gsnap.web.log' :::: diploid.web.list

parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.1 -A sam {}_1.fq.gz {}_2.fq.gz > {}.sam 2>> {}.gsnap.log' :::: AD1.list

parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.2 -A sam {}_1.web.fq.gz {}_2.web.fq.gz > {}.web.sam 2>> {}.gsnap.web.log' :::: AD2.web.list
parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.2 -A sam {}_1.fq.gz {}_2.fq.gz > {}.sam 2>> {}.gsnap.log' :::: AD2.list

parallel -j3 'gsnap --gunzip -n 1 -N 1 --merge-distant-samechr -Q -t 50 --use-sarray=0 -D db -d D5 -v snp4.4 -A sam {}_1.fq.gz {}_2.fq.gz > {}.sam 2>> {}.gsnap.log' :::: AD4.list

