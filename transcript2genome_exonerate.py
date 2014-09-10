#!/usr/bin/env python

import argparse,sys,os,subprocess,shutil,errno
from Bio import SeqIO

helptext='''
This script is intended to use the results of a BLAST search of a transcriptome 
against a draft genome to conduct more specific exonerate alignments. Each transcript 
is paired with its genome scaffold hit using the BioPython Indexing feature.
'''

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def parse_blast(blast_fn):
	'''Given a tabular BlAST output, return a dictionary of transcript:genome_scaffold hits.'''
	blast_dict = {}
	with open(blast_fn) as blast_file:
		for line in blast_file:
			line = line.split()
			try:
				blast_dict[line[0]].append(line[1])
			except KeyError:
				blast_dict[line[0]] = [line[1]]
	
	return blast_dict

def fasta_index(fasta_fn):
	return SeqIO.index(fasta_fn,'fasta')

def sort_blast(blast_dict):
	"""Given the blast dictionary, create a new dictionary mapping all transcripts to each scaffold"""
	scaffold_dict = {}
	for transcript in blast_dict:
		for scaffold in blast_dict[transcript]:
			if scaffold in scaffold_dict:
				if transcript not in scaffold_dict[scaffold]:
					scaffold_dict[scaffold].append(transcript)
			else:
				scaffold_dict[scaffold] = [transcript]
	return scaffold_dict

def fasta_generator(scaffold_dict,genome_idx,transcriptome_idx):
	"""Given the scaffold-transcript dictionary, generate scaffold/transcript fasta file pairs for parallel to work on."""
	scaffold_num = 1
	for scaffold in scaffold_dict:
			scaffold_fn = "exonerate_tmp/scaffold_{}.fasta".format(scaffold_num)
			transcript_fn = "exonerate_tmp/transcripts_{}.fasta".format(scaffold_num)
			transcript_seqs = [transcriptome_idx[t] for t in scaffold_dict[scaffold]]
			
			SeqIO.write(genome_idx[scaffold],scaffold_fn,'fasta')
			SeqIO.write(transcript_seqs,transcript_fn,'fasta')
			scaffold_num += 1
	return scaffold_num - 1

def parallel_exonerate(num_scaffolds,exonerate_format,exonerate_bestn,num_cpu = 4):
	#completedscaffolds = 0
	if exonerate_format == "human":
		exonerate_call = "exonerate -m cdna2genome -t exonerate_tmp/scaffold_{{}}.fasta -q exonerate_tmp/transcripts_{{}}.fasta --showvulgar no -V 0 -n {}".format(exonerate_bestn)
	elif exonerate_format == "gff":
		exonerate_call = "exonerate -m cdna2genome -t exonerate_tmp/scaffold_{{}}.fasta -q exonerate_tmp/transcripts_{{}}.fasta --showvulgar no --showalignment no --showtargetgff -V 0 -n {}".format(exonerate_bestn)
		
	parallel_call = "parallel -j {} --eta {} ::: {{1..{}}}".format(num_cpu,exonerate_call,num_scaffolds)
	sys.stderr.write(parallel_call)	
	subprocess.call(parallel_call,shell=True)
	#	completedScaffolds += 1
	#	progress_bar(float(completed_scaffolds)/num_scaffolds)



def pairwise_exonerate(genome_idx,transcriptome_idx,blast_dict,exonerate_format,exonerate_bestn):
	numTranscripts = float(len(blast_dict))
	completedTranscripts = 0
	for transcript in blast_dict:
		for scaffold in blast_dict[transcript]:
			SeqIO.write(transcriptome_idx[transcript],"temp_transcript.fasta","fasta")
			SeqIO.write(genome_idx[scaffold],"temp_scaffold.fasta","fasta")
		
			if exonerate_format == "human":
				exonerate_call = "exonerate -m cdna2genome -t temp_scaffold.fasta -q temp_transcript.fasta --showvulgar no -V 0 -n {}".format(exonerate_bestn)
			elif exonerate_format == "gff":
				exonerate_call = "exonerate -m cdna2genome -t temp_scaffold.fasta -q temp_transcript.fasta --showvulgar no --showalignment no --showtargetgff -V 0 -n {}".format(exonerate_bestn)
	
			subprocess.call(exonerate_call,shell=True)
		completedTranscripts += 1
		progress_bar(float(completedTranscripts)/numTranscripts)
		
def progress_bar(progress,barLength=10):
	"""Print to stderr the progress of a function"""
	status = ''
	if isinstance(progress,int):
		progress = float(progress)
	
	if progress < 0:
		progress = 0
		status = "Halt... \r\n"
	if progress >= 1:
		progress = 1
		status = "Done... \r\n"	
	
	block = int(round(barLength*progress))	
	text = "\rPercent: [{0}] {1:.2f}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
	sys.stderr.write(text)
	sys.stderr.flush()

def main():
	parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-b","--blast",dest="blast_fn",help="Tabular BLAST output, preferably using max_target_seqs=1",default=None)
	parser.add_argument("-t","--transcriptome",dest="transcriptome_fn", help = "Transcriptome FASTA file.",default=None)
	parser.add_argument("-g","--genome",dest="genome_fn",help = "Genome Scaffold FASTA file",default=None)
	
	parser.add_argument('--format',dest="exonerate_format",help="Exonerate output format type: [human, vulgar, gff, cigar]",choices=["human","vulgar","gff","cigar"],default='human')
	parser.add_argument('--bestn',dest="exonerate_bestn",help="bestn flag for Exonerate",default=1,type=int)
	parser.add_argument('--cpu',dest='num_cpu',help='Number of CPUs to run GNU parallel.',default=4,type=int)
	parser.add_argument('--no-cleanup',dest='cleanup',help="Don't delete temporary files generated for Exonerate.",action='store_false')
	
	parser.set_defaults(cleanup=True)
	
	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)
	args = parser.parse_args()
	
	input_files = [args.blast_fn,args.transcriptome_fn,args.genome_fn]
	for f in input_files:
		if not f:
			sys.stderr.write("Please specify all three input file names")
			parser.print_help(sys.stderr)
			sys.exit(1)
		if not os.path.isfile(f):
			sys.stderr.write("Cannot find file {}".format(f))
			sys.exit(1)
	
	genome_idx = fasta_index(args.genome_fn)
	sys.stderr.write("{} scaffolds loaded\n".format(len(genome_idx)))
	transcriptome_idx = fasta_index(args.transcriptome_fn)
	sys.stderr.write("{} transcripts loaded\n".format(len(transcriptome_idx)))
	blast_dict = parse_blast(args.blast_fn)
	sys.stderr.write("{} blast results loaded\n".format(len(blast_dict)))
	
	mkdir_p('exonerate_tmp')
	scaffold_dict = sort_blast(blast_dict)
	sys.stderr.write("{} scaffolds with hits, writing temporary files\n".format(len(scaffold_dict)))
	num_scaffolds = fasta_generator(scaffold_dict,genome_idx,transcriptome_idx)
	parallel_exonerate(num_scaffolds,args.exonerate_format,args.exonerate_bestn,num_cpu=args.num_cpu)
	
	if args.cleanup:
		shutil.rmtree('exonerate_tmp/')
	
	
	#pairwise_exonerate(genome_idx,transcriptome_idx,blast_dict,args.exonerate_format,args.exonerate_bestn)
			

if __name__ == '__main__':main()