#!/bin/bash

# Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

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
module load samtools/1.9-k6deoga
module load parallel/20170322-36gxsog

# https://github.com/IGBB/polycat; version 2.0.0b
ml subread/1.6.0-ak6vxhs

ls *.sam | parallel "samtools view -bu -@ 2 {} | samtools sort -m 4G -o {.}.sort.bam" 

ls [AD][25]*.sort.bam | parallel "./polycat D13.snp4.0 {}"
ls TM1*.sort.bam | parallel "./polycat D13.snp4.1 {}"
ls PS7*.sort.bam | parallel "./polycat D13.snp4.2 {}"
ls AD4*.sort.bam | parallel "./polycat D13.snp4.4 {}"

featureCounts -a D5.primaryOnly.gtf -o salt.counts -g 'gene_id' --primary  -p  -C  -T 15 --ignoreDup *.bam
