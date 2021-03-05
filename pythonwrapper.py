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




#Step 3

#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code
def use_Kallisto(SRR):
    kallisto_idx = 'time kallisto index -i hcmv_index.idx EF999921_CDS.fasta'
    os.system(kallisto_idx)
    kallisto_quant = 'time kallisto quant -i hcmv_index.idx -o kallisto_output/' + str(SRR) +' -b 30 -t 4 '+ str(SRR) + '_1.fastq '+ str(SRR)+ '_2.fastq'
















fasta_file()
use_Kallisto(SRR)
