#!/bin/bash

# Example using RSEM to quantify the V'DJer results
# Requires samtools to be in path.
# Run this after running demo.bash
# Quantification results will appear in: rsem_results.isoforms.results

# Set the following to your RSEM installation directory
RSEM=/path/to/rsem-1.2.21

if [ -s "vdj_contigs.fa" ]; then
  samtools view -1 -bS  vdjer.sam > vdjer.bam

  # Sort by name, while keeping read pairs together for input into RSEM
  # Note: you may need to increase the memory passed to sort for larger datasets.
  samtools view -H vdjer.bam > vdjer.namesort.sam
#  samtools view vdjer.bam | tr '\t' '~' | paste - - | sort -S 20G | awk '{print $1 "\n" $2}' | tr '~' '\t' >> vdjer.namesort.sam
  samtools view vdjer.bam | tr '\t' '~' | paste - - | sort -S 3G | awk '{print $1 "\n" $2}' | tr '~' '\t' >> vdjer.namesort.sam
  samtools view -bS vdjer.namesort.sam > vdjer.namesort.bam

  # RSEM prepare ref
  $RSEM/rsem-prepare-reference  vdj_contigs.fa rsem_vdj

  # RSEM calc expression
  $RSEM/rsem-calculate-expression -p 8 --paired-end --bam vdjer.namesort.bam rsem_vdj rsem_results
else
  touch rsem_results.isoforms.results
fi