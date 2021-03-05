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
SRRs = ['SRR5660030', 'SRR5660033','SRR5660044','SRR5660045']

def fastq(SRRs):
    for SRR in SRRs:
        SRR_link = 'https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/'+ SRR + '/'+ SRR + '.1'
        wget = 'wget' + ' ' + SRR_link
        fastq_dump = 'fastq-dump -I --split-files' + ' ' + SRR + '.1'
        os.system(wget)
        os.system(fastq_dump



#Step 2

#define function to create fasta file from CDS of NCBI Accession EF999921
#made a path to py file for 2nd step
def fasta_file():
  getfasta = "/usr/bin/python3 miniproject2.py" #path to seperate py file in hcmv repo
  os.system(getfasta)




#Step 3

#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code


def kallisto():
  run_idx = "kallisto index -i hcmv_index.idx EF999921_CDS.fasta"
  os.system(run_idx)
  for SRR in SRRs:
    path = os.cwd()
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + SRR + ' -b 30 -t 4 ' + SRR + '.1_1.fastq ' + SRR + '.1_2.fastq'
    os.system(run_kallisto)





















fasta_file()
kallisto()
