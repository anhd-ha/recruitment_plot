import os, sys, re, subprocess, shlex
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import pandas as pd
import seaborn
from matplotlib import pyplot as plt

### Sample command: python LAST_recruitment_plot_v4.py fastq_reads reference.fna output_prefix "#000000"

inputfile = sys.argv[1] # the FASTQ file to map
reference = sys.argv[2] # the contig file to map against
output = sys.argv[3]    # prefix for all output files
plot_color= sys.argv[4] # HEX color of reads in plot (in quotation marks)

# CUTOFFS: Change these if you want to make the search more or less stringent
bit_cutoff = 50
alnlength_cutoff = 50
percid_cutoff = 70

# make LAST db
refdb = reference +".lastdb"
cmd = "lastdb "+ refdb +" "+ reference
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stdout.txt", "w"), stdout=open("stderr.txt", "w"))

# run LAST search
cmd = "lastal -P 32 "+ refdb +" "+ inputfile +" -f BlastTab -Q 1"
#cmd = "lastal "+ refdb +" "+ inputfile +" -f BlastTab -Q 1"
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdin=open("stderr.txt", "w"), stdout=open(output+".lastout", "w"))

# get dictionary of contig lengths
contig2length = {}
for record in SeqIO.parse(reference, "fasta"):
	contig2length[record.id] = len(record.seq)

# get added values that should be added on to alignment coordinates
contig_bounds = open(output+".contig_boundaries.txt", "w")
sorted_dict = sorted(contig2length, key=contig2length.get, reverse=True)
added_values = defaultdict(float)
adj = 0
for i in sorted_dict:
	#print(i, contig2length[i])
	added_values[i] = adj
	adj = adj + contig2length[i]
	contig_bounds.write(i +"\t"+ str(added_values[i]) +"\n")
contig_bounds.close()

# parse through lastout file and generate final coordinates and output file
handle = open(output+".lastout", "r")
bit_dict = defaultdict(float)
hit2coord = {}
hit2percid = {}
for i in handle.readlines():
	if i.startswith("#"):
		pass
	else:
		line = i.rstrip()
		tabs = line.split("\t")
		query = tabs[0]
		contig = tabs[1]
		percid = float(tabs[2])
		score = float(tabs[11])
		aln_length = float(tabs[3])

		if score > bit_cutoff and percid > percid_cutoff and  aln_length > alnlength_cutoff and score > bit_dict[query]:
			bit_dict[query] = score
			start = float(tabs[8])
			end = float(tabs[9])
			mean = np.mean([start, end])
			final_coord = mean + added_values[contig]
			hit2coord[query] = final_coord
			hit2percid[query] = percid

coord_out = open(output+".finalcoords.txt", "w")
coord_out.write("readid\tcoord\tpercid\n")
for j in hit2coord:
	coord_out.write(j +"\t"+ str(hit2coord[j]) +"\t"+ str(hit2percid[j]) +"\n")

coord_out.close()


# function for seaborn plot  
def plot_recruitment_plot(coord_output, contig_bound, plot_out, plot_title):
    seaborn.set(color_codes=True)
    seaborn.set_style("ticks")
    fig = plt.figure()
    ax = plt.subplot(111)
 
    reads = pd.read_csv(coord_output, delimiter="\t")
    ax.scatter(reads.coord, reads.percid, s=0.2, color=str(plot_color))

    # draw contig boundaries 
    with open(contig_bound) as bounds:
#        next(bounds) # skip position 0
        for bound in bounds:
            bound = bound.split("\t")
            contig_bound = float(bound[1])
            ax.axvline(x=contig_bound, color="#545453")

    plt.ylabel('% Identity')
    plt.xlabel('Genome Position (bp)')
    plt.title(plot_title)
    fig.savefig(plot_out, format="png", dpi=400)
 
    plt.close()

# Generate recruitment plot
coord = str(output)+".finalcoords.txt"
contig = str(output)+".contig_boundaries.txt"
out = str(output)+".png"
title = str(output)
plot_recruitment_plot(coord, contig, out, title)


# cleanup
os.remove(output+".lastout")
os.remove("stderr.txt")
os.remove("stdout.txt")
for j in os.listdir("."):
	if j.startswith(reference +".lastdb"):
		os.remove(j)













