#import necessary modules for HCMV project
import os
import glob 
import argparse
from Bio import SeqIO
from Bio import Entrez


#Define a function for each tool/step in wrapper 

#Step 1 

#Convert SRA files from NCBI to fastq paired-end reads. Construct path based on SRR numbers
mainoutputfile = open("miniProject.log", "w+") #produce output file requested called miniProject.log for whole wrapper
path = os.getcwd() #returns current working directory of all code processes


def get_ncbi_data(SRR):
    wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
    rename = 'mv ' + SRR + '.1 ' + SRR
    os.system(wget)
    os.system(rename)


def run_fastq(SRR):
    fastqfiles = 'fastq-dump -I --split-files ' + SRR + '.1'
    os.system(fastqfiles)
    print(fastqfiles)









  

#include argparse to parse input data from the command line
parser = argparse.ArgumentParser(description='input SRR files and split paired reads.')
parser.add_argument('SRR', metavar='N', type=str, nargs='+', help='SRRs you want to assemble and BLAST')
parser.add_argument('--input', action='store_true', help='run with input data')
args = parser.parse_args()




if not args.input:
    for i in args.SRR:
        if i not in files:
            get_ncbi_data(i)
            run_fastq(i)
