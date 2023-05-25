Commands and scripts used in Laine et al. 2021 Does arsenic pollution influence DNA methylation patterns in a wild bird population? An experimental approach 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8277128/

Workflow:  
QC: Trimgalore  
Align and call: Bismark  
Analysis: Methylkit 

1. 01_QC_trimgalore.sh 

trim_galore --fastqc --quality 20 --paired --rrbs --basename sample1 1_Pmaj_1_S1_L001_R1_001.fastq.gz 1_Pmaj_1_S1_L001_R2_001.fastq.gz 

Before and after FastQC with multiQC. 

2. 02_alignment_bismark.sh 

bismark --genome /projects/methylation_arsenic/genome -1 sample20_val_1.fq.gz -2 sample20_val_2.fq.gz 

3. 03_methcalling_bismark.sh 

Before this script I did couple meth calls without the --ignore and checked the mbias to see how many bases needs to be trimmed. I followed the Bismark manual recommendations in this. 

4. 04_analysis_methylkit_combinedCs.R 
