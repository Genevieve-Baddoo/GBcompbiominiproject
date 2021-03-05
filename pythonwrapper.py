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

#Step 3

#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code
def kallisto_idx():
  run_idx = "kallisto index -i hcmv_index.idx EF999921_CDS.fasta"
  os.system(run_idx)



path = os.getcwd()

def kallisto_quant(SRRs):
  for SRR in SRRs:
    quant = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + SRR + ' -b 30 -t 4 ' + SRR + '.1_1.fastq ' + SRR + '.1_2.fastq'
    os.system(quant)


#define function to generate the sleuth input
def sleuth_input(SRRs):
    # output file that goes in R
    output = open('sleuth_infile.txt', 'w')
    # initial line in file
    output.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')
    # based on SRR number, write condition and path to output file
    path = os.getcwd()
    for SRR in SRRs:
        path1 = path + '/' + 'results_' + SRR
        print(path1)
        if int(SRR[
               3:]) % 2 == 0:  
            output.write(str(SRR) + '\t' + '2dpi' + '\t' + path1 + '\n')
        else:
            output.write(str(SRR) + '\t' + '6dpi' + '\t' + path1 + '\n')
    output.close()









#run functions
fastq(SRRs)
fasta()
kallisto_idx()
kallisto_quant(SRRs)
sleuth_input(SRRs)

