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




#want to parse
parser = argparse.ArgumentParser()
parser.add_argument("--SRR", nargs="+", required=True, help = "SRR input to assemble and BLAST")
parser.add_argument("--fastq", action="store_true", help = "SRR input data in paired end fastq format")
p = parser.parse_args()





SRR = []
#get SRRs from argument
for SRR in p.SRR:
  SRR.append(SRR)




  


get_ncbi_data(SRR)
run_fastq(SRR)

