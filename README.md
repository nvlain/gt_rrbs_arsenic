Commands and scripts used in Laine et al. 2021 Does arsenic pollution influence DNA methylation patterns in a wild bird population? An experimental approach 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8277128/

Workflow:  
QC: Trimgalore  
Align and call: Bismark  
Analysis: Methylkit 

1. 01_QC_trimgalore.sh 

trim_galore --fastqc --quality 20 --paired --rrbs --basename sample1 /mnt/nfs/bioinfdata/data_other/genomes/NIOO/AnE/Parus_major/raw_data/reduced_library_sequencing/RRBS/Arsenic_RRBS_SuviRuuskanen/1_Pmaj_1_S1_L001_R1_001.fastq.gz /mnt/nfs/bioinfdata/data_other/genomes/NIOO/AnE/Parus_major/raw_data/reduced_library_sequencing/RRBS/Arsenic_RRBS_SuviRuuskanen/1_Pmaj_1_S1_L001_R2_001.fastq.gz;
t 

Before and after FastQC with multiQC. 

2. 02_alignment_bismark.sh 

bismark --genome /mnt/nfs/bioinfdata/home/NIOO/veronikal/projects/methylation_arsenic/genome -1 /mnt/nfs/bioinfdata/home/NIOO/veronikal/projects/methylation_arsenic/veronika_protocol/01_QC/sample20_val_1.fq.gz -2 /mnt/nfs/bioinfdata/home/NIOO/veronikal/projects/methylation_arsenic/veronika_protocol/01_QC/sample20_val_2.fq.gz; 

3. 03_methcalling_bismark.sh 

Before this script I did couple meth calls without the --ignore and checked the mbias to see how many bases needs to be trimmed. I followed the Bismark manual recommendations in this. 

4. 04_analysis_methylkit_combinedCs.R 
