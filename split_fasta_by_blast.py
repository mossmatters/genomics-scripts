#!/usr/bin/env python

#A general script for splitting a FASTA file based on (tabular) BLAST results.

import sys
from Bio import SeqIO

if len(sys.argv) < 3:
	print "Usage: python split_fasta_by_blast.py blastfile fastafile"
	print "Blast must have been done with the same fasta file"
	print "Will output a 'hit' and a 'nothit' file to the current directory"

blast_fn = sys.argv[2]
fasta_fn = sys.argv[1]

basename = blast_fn.split(".")[0]
targetdb = blast_fn.split('.')[1]

hits_fn = "{}.{}.hit.fasta".format(basename,targetdb)
nothit_fn = "{}.{}.nothit.fasta".format(basename,targetdb)

hits_file = open(hits_fn,'w')
nothits_file = open(nothit_fn,'w')

blast_results = {line.split()[0]:line.split()[1:] for line in open(blast_fn).readlines()}
hit_ids = set(blast_results.keys())

hit_seqs = 0
nothits_seqs =0 
for prot in SeqIO.parse(fasta_fn,'fasta'):
	if prot.id in hit_ids:
		#Filter for hits that are at least 50% of the length of the query protein
		if int(blast_results[prot.id][2]) > len(prot.seq) * 0:#.5:
			hit_seqs += 1
			SeqIO.write(prot,hits_file,'fasta')
		else:
			nothits_seqs += 1
			SeqIO.write(prot,nothits_file,'fasta')
	else:
		nothits_seqs += 1
		SeqIO.write(prot,nothits_file,'fasta')

print "{}\nProtein Hits: {}\nNot Hit: {}".format(fasta_fn,hit_seqs,nothits_seqs)


hits_file.close()
nothits_file.close()