### Contents

The genome used (filename: Dgenome2_13.fasta) is the the D5_JGI genome (https://www.cottongen.org/species/Gossypium_raimondii/jgi_genome_221) with only chromosomes.

The files D13.snp4.x are the original SNP indices from doi: 10.1534/g3.112.005298    
D13.snp4.0: SNP index for the diploids, _G. arboreum_ and _G. raimondii_     
D13.snp4.1: SNP index for _G. hirsutum_     
D13.snp4.2: SNP index for _G. barbadense_     
D13.snp4.4: SNP index for _G. mustelinum_     

1_gsnap.slurm: preparation of gsnap files and read mapping      
2_sam2bam.slurm: transform sam2bam, polycat to partition, and featureCounts to count reads

salt.counts.gz: count file output by subread     
salt.counts.summary.gz: summary file output by subread

*.list files: lists of reads categorized as diploid, _G. hirsutum_, _G. barbadense_, or _G. mustelinum_

D5.primaryOnly.gtf: GTF file of annotations derived from https://www.cottongen.org/species/Gossypium_raimondii/jgi_genome_221. Only primary transcripts are considered.
