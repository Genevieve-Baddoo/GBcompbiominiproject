#import necessary modules for HCMV project
import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW

#Define a function for each tool/step in wrapper 

#Step 1 

#Convert SRA files from NCBI to fastq paired-end reads. Construct path based on SRR numbers
#defined SSRs as a list
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
#print(SRR)

def fastq(SRR):
  wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
  fastq = 'fastq-dump -I --split-files '+ str(SRR) #what you would use in terminal command line
  rename = 'mv '+ SRR + '.1 ' + SRR
  os.system(wget) #os.system runs the command line
  os.system(fastq)
  os.system(rename)

for i in SRR:
  fastq(i)


#Step 2

#define function to create fasta file from CDS of NCBI Accession EF999921
#made a path to py file for 2nd step
def fasta_file():
  getfasta = "/usr/bin/python3 miniproject2.py" #path to seperate py file in hcmv repo
  os.system(getfasta)





















fasta_file()
