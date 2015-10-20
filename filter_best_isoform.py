#!/usr/bin/env python

import sys,os,argparse
from Bio import SeqIO

helptext = """This script will take the output of Trinity's RSEM scripts and return
the best isoform from the original Trinity.fasta file for each gene. The user may
select the definition of 'best' as being TPM, FPKM, or IsoPct. You may also optionally
specify a cutoff for each isoform (e.g. FKPM > 2). Results go to stdout.
"""

def get_column_from_filter_mode(filtering_mode):
	filter_modes = {'TPM':5,
					'FPKM':6,
					'IsoPct':7,
					'length':2
					}
	return filter_modes[filtering_mode]
					
def parse_rsem(rsem_output_filename,filtering_mode,cutoff):
	"""Return a list of the best isoforms from each geneID"""
	rsem_dict = {}
	filter_column = get_column_from_filter_mode(filtering_mode)
	rsem_output_file = open(rsem_output_filename)
	header = rsem_output_file.readline()
	
	for line in rsem_output_file.readlines():
#		print rsem_dict
		line = line.split()
		geneID = line[1]
		filter_value = float(line[filter_column])
		try:
			if filter_value > rsem_dict[geneID][1]:
				if filter_value > cutoff:
					rsem_dict[geneID] = (line[0],filter_value)
		except KeyError:
			if filter_value > cutoff:
				rsem_dict[geneID] = (line[0],filter_value)
	rsem_output_file.close()
	return [rsem_dict[x][0] for x in rsem_dict]



def filter_isoforms(rsem_output_filename,trinity_output_filename,filtering_mode,cutoff):
	"""Prints the best isoform from each geneID in a Trinity.fasta file, using the RSEM data"""
	isoform_list = parse_rsem(rsem_output_filename,filtering_mode,cutoff)
	trinity_transcripts = SeqIO.to_dict(SeqIO.parse(trinity_output_filename,'fasta'))
	numseq=0
	for isoform in isoform_list:
		numseq += SeqIO.write(trinity_transcripts[isoform],sys.stdout,'fasta')
	return numseq	
	
	

def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("--rsem_output","-r",help="RSEM isoform output file.",default=None)
	parser.add_argument("--trinity_output","-t",help="Trinity output file",default=None)
	parser.add_argument("--filtering_mode","-f",help="Filtering mode, will choose the highest value in either 'length','FPKM', 'TPM', or 'IsoPct'",default='FPKM')
	parser.add_argument("--cutoff",type=float,help="Optional filtering cutoff using the same filteirng mode.",default=None)
	
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)
	args = parser.parse_args()
	
	if not os.path.isfile(args.rsem_output):
		sys.stderr.write("RSEM output file not found at: {}".format(args.rsem_output))
		sys.exit(1)
	if not os.path.isfile(args.trinity_output):
		sys.stderr.write("Trinity fasta file not found at: {}".format(args.trinity_output))
		sys.exit(1)
	trinity_transcripts = filter_isoforms(args.rsem_output,args.trinity_output,args.filtering_mode,args.cutoff)		
	sys.stderr.write('{} sequences written\n'.format(trinity_transcripts))
if __name__ == '__main__':main()