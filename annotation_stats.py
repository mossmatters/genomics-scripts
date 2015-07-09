#!/usr/bin/env python
# Read in two ublast results in blast tabular format:
# 	Transcriptome assembly to Protein database (Blastx-like)
# 	Protein database to transcriptome assembly (tBlastN-like)
#
#  The script also needs the location of the original protein file, 
#		to get the Subject length for Ortholog Hit Ratios
 	
# Output statistics:
# 	Number of unique best hits in assembly
# 	Number of proteins with hits in assembly
# 	Number of Reciprocal Best Hits
# 	Ortholog Hit Ratio


from __future__ import division
import sys,argparse
from Bio import SeqIO

helptext = '''
Usage: python annotation_stats.py -b blastx.results -t tblastn.results -p reference_proteins.fasta
'''

def input_files_exist(file_list):
	import os
	for f in file_list:
		if os.path.isfile(f):
			pass
		else:
			sys.stderr.write("File {} not found!".format(f))
			return False
	return True

def read_subset(subset_filename):
	return set([x.rstrip() for x in open(subset_filename).readlines()])

def get_protein_dict(protein_filename):
	return SeqIO.to_dict(SeqIO.parse(protein_filename,'fasta'))

def blastx_parse(blastx_filename,protein_dict,subset_fn=None):
	if subset_fn:
		subset = read_subset(subset_fn)
	blastx_dict = {}
	OHR = []
	
	blastx_file = open(blastx_filename)
	with open(blastx_filename) as blastx_file:
		for line in blastx_file:
			hit = line.split('\t')
	
	#NOTE: Change the next two lines if the delimiters are different for your files!
			#contig = hit[0]#.split(' ')[0]
			#peptide = hit[1] #.split('|')[0].upper()
			
			contig,peptide = hit[0:2]
	#OHR = Alignment Length / Subject Length	
			OHR.append(float(hit[3])/len(protein_dict[peptide].seq))
			blastx_dict[contig] = peptide
	return blastx_dict,OHR

def tblastn_parse(tblastn_filename,subset_fn = None):
	if subset_fn:
		subset = read_subset(subset_fn)
	tblastn_dict = {}
	with open(tblastn_filename) as tblastn_file:
		for line in tblastn_file:
			hit = line.rstrip().split('\t')
	
	#NOTE: Change the next two lines if the delimiters are different for your files!
			peptide = hit[0].rstrip() #.split("|")[0].upper()
			contig = hit[1] #.split("|")[0].upper()
	#print contig,peptide
			tblastn_dict[peptide] = contig
	return tblastn_dict

def find_reciprocal_best_hits(blastx_dict,tblastn_dict):
	RBH = 0
	for i in blastx_dict.keys():
		try:
			j = tblastn_dict[blastx_dict[i]]
			if i == j:
				RBH += 1
		except KeyError:
			continue
	return RBH
	
def find_collapse_factor(tblastn_dict):
	cf_dict={}
	for prot in tblastn_dict:
		contig = tblastn_dict[prot]
		if contig in cf_dict:
			cf_dict[contig] +=1
		else:
			cf_dict[contig] = 1
	avg_cf = sum(cf_dict.values())/len(cf_dict)
	return avg_cf
	
def print_stats(blastx_dict,OHR,tblastn_dict,RBH,avg_cf):
	print "Annotation statistics"
	print "---------------------"	
	print "%i contigs had a best hit in the protein database!" % len(blastx_dict)
	print "There were %i unique peptides found." % len(set(blastx_dict.values()))
	print "The average Ortholog Hit Ratio was %.2f" % (sum(OHR)/len(OHR))
	print "\nReverse Search Statistics"
	print "---------------------------"
	print "%i peptides had a best hit in the assembly!" %len(tblastn_dict)
	print "There were %i unique contigs with hits." % len(set(tblastn_dict.values()))	
	print "There were %i reciprocal best hits!" % RBH	
	print "The average collapse factor was %.2f" % avg_cf
	

def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-b','--blastx_fn', help="BLASTx tabular format file, standard fields",default=None)
	parser.add_argument('-t','--tblastn_fn', help="tBLASTn tabular format file, standard fields",default=None)
	parser.add_argument('-p','--protein_fn',help="FASTA file of reference proteins",default=None)
	#parser.add_argument('--sub',help="Optional list of transcript IDs to use for calculations.",default=None)
	
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()
	
	if args.blastx_fn and args.tblastn_fn and args.protein_fn:
		if input_files_exist([args.blastx_fn,args.tblastn_fn,args.protein_fn]):
			blastx_filename = args.blastx_fn
			tblastn_filename = args.tblastn_fn
			protein_filename = args.protein_fn
		else:
			sys.exit(1)
	else:
		parser.print_help()
		sys.exit(1)
	
	protein_dict = get_protein_dict(protein_filename)
	blastx_dict,OHR = blastx_parse(blastx_filename,protein_dict)	
	tblastn_dict = tblastn_parse(tblastn_filename)	

	RBH = find_reciprocal_best_hits(blastx_dict,tblastn_dict)
	avg_cf = find_collapse_factor(tblastn_dict)
	
	print_stats(blastx_dict,OHR,tblastn_dict,RBH,avg_cf)

if __name__ == '__main__':main()

#Read the protein file into a Dict for easy hashing.



#Read the blastx-like results

#Read tblastn-like results

#print blastx_dict.keys()[0:10]
#print tblastn_dict.keys()[0:10]

#Find Reciprocal Best Hits

#Find Collapse Factor

	
