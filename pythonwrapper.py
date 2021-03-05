#import necessary modules for HCMV project
import os
import os.path
from os import path
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



kallisto_output = "/homes/gbaddoo/hcmvminiproject"
os.mkdir(kallisto_output)

if path.exists('kallisto_output'):
    print ("File exist")
else:
    os.mkdir(kallisto_output)


#Step 3

#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code
def use_Kallisto(SRR):
    kallisto_idx = 'time kallisto index -i hcmv_index.idx EF999921_CDS.fasta'
    os.system(kallisto_idx)
    kallisto_quant = 'time kallisto quant -i hcmv_index.idx -o kallisto_output/' + str(SRR) +' -b 30 -t 4 '+ str(SRR) + '_1.fastq '+ str(SRR)+ '_2.fastq'




#define function to generate the sleuth input
def sleuth_input(SRR):
    #input file for sleuth package (read this txt file in read.table for R script)
    outfile = open('sleuth_infile.txt','w') 
    #initial line in file (column names)
    outfile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to outnput file
    for i in SRR:
        path = 'kallisto_output/' + i
        if SRR.index(i)%2==0:
          outfile.write(str(i) + '\t' + '2dpi' + '\t' + path + '\t' + '\n')
        else:
          outfile.write(str(i) + '\t' + '6dpi' + '\t' + path + '\t' + '\n')
    outfile.close()


#define function to run sleuth in R
#want to read sleuth output and add to miniProject.log file
def sleuth():
    sleuth_command = 'Rscript sleuth.R'
    os.system(sleuth_command)
    log_file = open("miniProject.log", "w")
    sleuth_read = open('sleuth_infile.txt', 'r').readlines()
    for i in sleuth_read:
      log_file.write(i + '\n')








fasta_file()
use_Kallisto(SRR)
sleuth_input(SRR)
sleuth()
