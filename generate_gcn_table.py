#!/usr/bin/python2.7

"""
NAME:        generate_gcn_table.py
AUTHOR:      Pascal Timshel, pascal.timshel@gmail.com 
DESCRIPTION: Program for generating GCN table
"""

from sys import argv
from argparse import ArgumentParser
import os, datetime, getpass

from Bio import SeqIO


def ParseArguments():
	""" Handles program parameters and returns an argument class containing 
	all parameters """

	parser = ArgumentParser(description='Description of program')

	parser.add_argument("-f", type=str, dest="fasta_file", 
                        default="", help="Fasta input file", required=True)
	parser.add_argument("-o", type=str, dest="output_file", 
                        default="", help="Output file")
	parser.add_argument("-p", type=float, dest="first_parameter", 
                        default=0.5, help="First Parameter (decimal)")
	parser.add_argument("-n", type=int, dest="second_parameter", 
                        default=0, help="Second Parameter (integer)")
	
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

def my_function(args):
	print "Input file is:", args.fasta_file
	print "Output file is:", args.output_file
	print "Paramter 1 is:", args.first_parameter
	print "Paramter 2 is:", args.second_parameter



def main(argv):	
    args = ParseArguments()

    # PUT FUNCTIONS AND CLASSES ABOVE AND USE THEM HERE
    my_function(args)

if __name__ == '__main__':
   main(argv)
