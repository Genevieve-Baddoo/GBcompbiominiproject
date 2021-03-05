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
        os.system(fastq_dump)


def fasta():
    getfasta = '/usr/bin/python3 miniproject2.py'
    os.system(getfasta)











#run functions

fastq(SRRs)
fasta()
