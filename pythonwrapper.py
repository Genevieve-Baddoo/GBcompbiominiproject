#import necessary modules for HCMV project
import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW



#Define a function for each tool/step in wrapper. 
#Add test data (small # of input reads) once done

#Convert SRA files from NCBI to fastq paired-end reads. Construct path based on SRR numbers
#create list with acession numbers listed
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

def fastq(SRR):
  wget = 'wget https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-12/SRR' + SRR + '/' + SRR + '.1'
  fastq = 'fastq-dump -I --split-files '+ str(SRR) #what you would use in terminal command line
  rename = 'mv '+ SRR + '.1 ' + SRR
  os.system(wget) #os.system runs the command line
  os.system(fastq)
  os.system(rename)

for i in SRR:
  fastq(i)


#define function to create fasta file from CDS of NCBI Accession EF999921
#made a path to py file for 2nd step
def fasta_file():
  getfasta = "/usr/bin/python3 miniproject2.py" #path to seperate py file in hcmv repo
  os.system(getfasta)



#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code
def use_Kallisto(SRR):
    kallisto_idx = 'time kallisto index -i hcmv_index.idx CDS_EF999921.fasta'
    os.system(kallisto_idx)
    kallisto_quant = 'time kallisto quant -i hcmv_index.idx -o ./' + str(SRR) +' -b 30 -t 4 '+ str(SRR) + '_1.fastq '+ str(SRR)+ '_2.fastq'
    os.system(kallisto_quant)


#define function to create file for input into R sleuth
def sleuthInput(SRR):
    #input file for sleuth package
    outfile = open('sleuth_infile.txt','w')
    #initial line in file (column names)
    outfile.write('sample'+ '\t' + 'condition' + '\t' + 'path' + '\n')
    #based on SRR number, write condition and path to outnput file
    for i in SRR:
        path = "/" + i
        if int(i[2:])%2==0:
            outfile.write(str(i)+'\t'+'2dpi'+path+'\t')
        else:
            outfile.write(str(i)+'\t'+'6dpi'+path+'\t')
    outfile.close()


#define function to run sleuth in R
#want to read sleuth output and adds to miniProject.log file
def sleuth():
    sleuth_command = 'Rscript sleuth.R'
    os.system(sleuth_command)
    sleu_out = "Rsleuthout.txt"
    readsleuth = open(sleu_out).readlines()
    for i in read_sleuth:
      outputfile.write(i + 'n')
    
#defined outputfile = open("miniProject.log", "w") in miniproject2.py


