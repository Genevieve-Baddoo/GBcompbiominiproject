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
def use_kallisto(SRR):
    kallisto_idx = 'time kallisto index -i HCMV_index.idx EF999921_CDS.fasta'
    os.system(kallisto_idx)
    path = os.getcwd()
    run_kallisto = 'time kallisto quant -i HCMV_index.idx -o' + path + '/results_' + str(SRR) + ' -b 30 -t 4 ' + str(SRR) + '.1_1.fastq ' + str(SRR) + '.1_2.fastq'
    os.system(run_kallisto)



def sleuth_input(SRR):
    # output file that goes in R
    output = open('sleuth_infile.txt', 'w')
    # initial line in file
    output.write('sample' + '\t' + 'condition' + '\t' + 'path' + '\n')
    # based on SRR number, write condition and path to output file
    path = os.getcwd()
    for i in SRR:
        path1 = path + '/' + 'results_' + i
        print(path1)
        if int(i[
               3:]) % 2 == 0:  
            output.write(str(i) + '\t' + '2dpi' + '\t' + path1 + '\n')
        else:
            output.write(str(i) + '\t' + '6dpi' + '\t' + path1 + '\n')
    output.close()

def sleuth():
    sleuth_command = 'Rscript sleuth.R'
    os.system(sleuth_command)
    log_file = open("miniProject.log", "w")
    read_sleuth = open('sleuth_outfile.txt', 'r').readlines()
    for i in read_sleuth:
        log_file.write(i + '\n')































fasta_file()
use_kallisto(SRR)
sleuth_input(SRR)
sleuth()
