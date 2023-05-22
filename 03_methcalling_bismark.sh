##bash
for f in ./02_alignment/*.bam; do
  bismark_methylation_extractor -p --no_overlap --report --bedGraph --scaffolds --cytosine_report --ignore 3 --ignore_r2 2 --ignore_3prime_r2 1 --genome_folder /mnt/nfs/bioinfdata/home/NIOO/veronikal/projects/methylation_arsenic/genome "$f" 
done