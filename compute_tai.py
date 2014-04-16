#!/usr/bin/python2.7

"""
NAME:        compute_tai.py
AUTHOR:      Pascal Timshel, pascal.timshel@gmail.com 
DESCRIPTION: Program for calculating tAI (tRNA adaptation index, dos Reis et al. 2004).

Improvements:
	- Passe references to objects or remake to class tAI	
	- Regex for ID in fasta header
	- Redefine print_tai_pos to use csv.write module
"""


from sys import argv
from argparse import ArgumentParser
import os, datetime, getpass

from math import exp, log
from collections import defaultdict

from Bio.Seq import Seq # for reverse complement etc.
from Bio import SeqIO


def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """

	parser = ArgumentParser(description='Script for calculating tAI (tRNA adaptation index, dos Reis et al. 2004)')

	parser.add_argument("-f", type=str, dest="fasta_file", 
                        default="", help="Fasta input file", required=True)
	parser.add_argument("-gcn", type=str, dest="gcn_tab", 
                        default="", help="tRNA gene copy number (GCN) file", required=True)
	parser.add_argument("-o", type=str, dest="output_file", 
                        default="tai.res", help="Output file for tAI calculation [default tai.res]")

	parser.add_argument("-p", dest="position_calc", action='store_true',
                         help="Calculate position specific tAI and print to file tai.pos") # default=False is implied by action='store_true'
	
	args = parser.parse_args()

	# PRINT RUNNING DESCRIPTION 
	now = datetime.datetime.now()
	print '# ' + ' '.join(argv)
	print '# ' + now.strftime("%a %b %d %Y %H:%M")
	print '# USER: ' + getpass.getuser()
	print '# CWD: ' + os.getcwd()
	if os.name != "nt":
		print '# HOST: ' + ' '.join([ os.uname()[0] , os.uname()[4] , 
        	                         os.uname()[2] , os.uname()[1] ])

	print '# COMMAND LINE PARAMETERS SET TO:'
	for arg in dir(args):
		if arg[:1]!='_':
			print '# \t' + "{:<30}".format(arg) +\
				  "{:<30}".format(getattr(args, arg))

	return args

def print_args(args):
	print "Input fasta file is:", args.fasta_file
	print "tRNA gene copy number table:", args.gcn_tab
	print "Output file for tAI table is:", args.output_file
	print "Calculate position specific tAI:", str(args.position_calc)


def read_gcn_tab(filename):
	""" Reads GCN file skipping first line (header).
	The file must be structured with the following three tab seperated columns: 
		| anticodon | codon | GCN | """
	tRNAtab = {}
	with open(filename, 'r') as f:
		next(f) # skipping header of table
		for line in f:
			line = line.rstrip('\n')
			cols = line.split('\t')
			anticdn = cols[0] # anticodon
			gcn = int(cols[2]) # gene copy number
			tRNAtab[anticdn] = gcn
			#print "anticdn:gcn ==", anticdn+":"+gcn
	return(tRNAtab)

def calc_tAIweights(tRNAtab):
	"""
	Function to calculate the weights described in dos Reis et al. 2004.
	The s_ij wobble parameters for the tAI calculation
	Special order of the s_X:Y wobble parameters
	1-4: A:T, G:C, T:A, C:G
	5-8: G:T, A:C, A:A, T:G
	9:           L:A
	"""
	sij = (1.00,  1.00,  1.00,   1.00, 
	        0.59,  0.72,  0.0001, 0.32, 0.11) # weights are optimized by dos Reis et al. 2004.
	wi = {}
	anticodons = []
	codons = []
	GCN = []
	W = []
	
	## the anticodons built by position in the order we need.
	acpos1 = "AGTC" * 16
	acpos2 =  "A" * 16 + "G" * 16 + "T" * 16 + "C" * 16
	acpos3 = ("A" *  4 + "G" *  4 + "T" *  4 + "C" *  4) * 4
	
	for i in range(len(acpos1)):
		ac = acpos1[i] + acpos2[i] + acpos3[i]
		anticodons.append(ac)
		ac_seq = Seq(ac) # bioPython
		cdn = ac_seq.reverse_complement()
		codons.append(str(cdn))
		GCN.append(tRNAtab.get(ac, 0))
		# KEEP for DEBUGGIN
		#if i == 0: print "ac\tcodons\tGCN"
		#print "%s\t%s\t%s" % (ac, codons[i], GCN[i])
	
	## TODO: make the below code work instead for string replication method
	#bases = ['T', 'C', 'A', 'G']
	#codons1 = [a+b+c for a in bases for b in bases for c in bases]
	#print "#######"	
	#print "\n".join(codons1)
		
	## Calculate 'W'
	for i in xrange(0,len(codons),4):
		#print sij[0], GCN[i], sij[4], GCN[i+1]
		#print type(sij[0]), type(GCN[i]), type(sij[4]), type(GCN[i+1])
		Wi0 = sij[0] * GCN[i] + sij[4] * GCN[i+1] # INN -> NNT, NNC, NNA
		Wi1 = sij[1] * GCN[i+1] + sij[5] * GCN[i]   # GNN -> NNT, NNC
		Wi2 = sij[2] * GCN[i+2] + sij[6] * GCN[i]   # TNN -> NNA, NNG
		Wi3 = sij[3] * GCN[i+3] + sij[7] * GCN[i+2] # CNN -> NNG
		W.extend([Wi0, Wi1, Wi2, Wi3])
	
	## CORRECTED (re-included 14,10 Ile, Val)
	## 11, 34, 35, 50 should be removed before the max is taken	
	i_remove = [11,34,35,50]
	i_used = range(0,64) # gives 0...63
	i_used = [x for x in i_used if x not in i_remove] # filtering list
	#OR i_used = range(0,11) + range(12,34) + range(36,50) + range(51,64)
	
	## get the max(W_i) and then the mean(w_i)
	## definition of the w_mean should not include '0's (!!!)
	W_used  = [W[i] for i in i_used]
	W_max = max(W_used)
	non0wi = 0;
	
	#Formula (2) in dos Reis 2004
	for i in i_used:
		cdn = codons[i]
		if W[i] != 0:
			non0wi += 1
			wi[cdn] = W[i] / float(W_max)
	
	# Calculating geometric mean of all w_i with W_i != 0
	# Log transforming first for numerical stability
	w_gmean = exp( sum([log(a) for a in wi.values()]) / non0wi )
	
	#Formula (2) in dos Reis 2004#
	## wi for TAA, TAG, TGA are 0
	for i in i_used:
		cdn = codons[i]
		#if wi[cdn] == 0: wi[cdn] = w_gmean
		if wi.get(cdn,0) == 0: wi[cdn] = w_gmean
		cdnSeq = Seq(cdn) # TODO delete?
		ac = cdnSeq.reverse_complement() # TODO delete?
	return(wi)

def count_codons(seq):
	""" Calculating codon counts """
	codon_counts = defaultdict(int) # important to use default dict
	cdn_list = [seq[i:i+3] for i in range(0,len(seq), 3)]
	for cdn in cdn_list:
		codon_counts[cdn] += 1
	return codon_counts


def tAI(cdn_count, wi):
	"""
	tRNA adaptation index (dos Reis et al. 2004)
	arguments:
	cdn_count: dict of codon counts for an mRNA
	wi: dict of w_i values based on tRNA GCN
	"""
	tai = 0
	n = 0
	for cdn in cdn_count.keys():
		if cdn in wi: # NO KEY/VALUE Error
			tai += log(wi[cdn]) * cdn_count[cdn]
			n += cdn_count[cdn]
	return ( exp(tai/n) )

def tAIbyCodonIndex(seq, wi):
	"""
	Position specific tAI: returns weights for each codon position
	"""
	tAI_by_position = []
	for i in xrange(0, len(seq)-2, 3): # i <= len(seq)-3
		cdn = seq[i:i+3]
		tAI_by_position.append(wi.get(cdn,0))
	return tAI_by_position


def parse_seq(filename, wi, pos_calc):
	"""
	Parsing fasta sequence and generates tAI dict for each entry
	"""
	stats = defaultdict(dict)
	
	for record in SeqIO.parse(filename, "fasta"):
		id = record.id
		mRNAseq = record.seq.transcribe()
		if len(mRNAseq) % 3 != 0:
			e_msg = "Length of mRNAseq is not a multiplum of 3. Sequence is NOT a coding seqeunce."
			raise Exception("Error in fasta entry ID =" + id + "\n" + e_msg + "\n")
		
		# Counting codons
		#stats[id]["codons"] = count_codons(str(record.seq))
		cdn_count = count_codons(str(record.seq))
		
		stats["tai"][id] = tAI(cdn_count, wi)
		
		if pos_calc:
			stats["tai_pos"][id] = tAIbyCodonIndex(str(record.seq), wi)
	
	return(stats)
		
def print_tAI_table(filename, tai_dict):
	with open(filename, "w") as f:
		f.write("ID\ttAI\n")
		for (id, tai) in sorted(tai_dict.items()):
			f.write("%s\t%.4f\n" % (id, tai))

def print_tAI_pos_table(filename, tai_pos_dict):
	"""
	This printing function is suited for importing data into a dataframe in R
	"""
	list_len = []
	for (id, tai_pos_list) in tai_pos_dict.items():
		list_len.append(len(tai_pos_dict[id]))
	max_list_len = max(list_len)
	with open(filename, "w") as f:
		#f.write("ID\t" + "\t" * (max_list_len-1) + "\n") # CHECK THIS
		f.write("ID\t" + "\t".join(str(x) for x in range(0,max_list_len)) + "\n") # CHECK THIS
		for (id, tai_pos_list) in sorted(tai_pos_dict.items()):
			f.write("%s\t" % id)
			#list_format = ["%.4f" % val for val in tai_pos_list]
			str_format  = ("%.4f\t" * (len(tai_pos_list)-1)) + "%.4f\n"
			f.write(str_format % tuple(tai_pos_list))

	
def main(argv):	
	args = ParseArguments()

	# Prints commandline arguments
	print_args(args)

	# Read GCN file
	tRNAtab = read_gcn_tab(args.gcn_tab)
	
	# Calculate weights from GCN table
	wi = calc_tAIweights(tRNAtab)
	
	# Read fasta file and calculate tai values
	tai_dict = parse_seq(args.fasta_file, wi, args.position_calc)
	
	# Print results	
	print_tAI_table(args.output_file, tai_dict["tai"])
	
	if args.position_calc:
		print_tAI_pos_table("tai.pos", tai_dict["tai_pos"])
	

if __name__ == '__main__':
	main(argv)
