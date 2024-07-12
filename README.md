# Recruitment plot
Script to map short read data to reference isolate, single-cell genomes, or metagenome-assembled genomes (MAGs) and output recruitment plot to assess microbial population relative abundance and/or structure.
This script uses [LAST search](https://genome.cshlp.org/content/21/3/487.long). It takes into account coordinates if multiple contigs are in the reference. 
Requires SeqIO from Biopython, defaultdict, numpy, pandas, seaborn, and matplotlib

## Usage: python LAST_recruitment_plot.py fastq_reads reference.fna output_prefix "#000000"

fastq_reads: the FASTQ read file to map
reference.fna: the MAG/contig file to map the reads against
output_prefix: prefix for all output files
"#000000": HEX color code (in quotation marks) for reads in plot 
