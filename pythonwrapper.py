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


#Step 1

#Convert SRA files from NCBI to fastq paired-end reads. Construct path based on SRR numbers
#create list with acession numbers listed
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

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



#kallisto_output = "/homes/gbaddoo/hcmvminiproject"
#os.mkdir(kallisto_output)

#if path.exists('kallisto_output'):
    #print ("File exist")
#else:
    #os.mkdir(kallisto_output)



#Step 3

#define function to take SRR # and create/run Kallisto on the command line
#use CDS fasta file created in the previous py code
def use_Kallisto(SRR):
    kallisto_idx = 'time kallisto index -i hcmv_index.idx EF999921_CDS.fasta'
    os.system(kallisto_idx)
    kallisto_quant = 'time kallisto quant -i hcmv_index.idx -o kallisto_output/' + str(SRR) +' -b 30 -t 4 '+ str(SRR) + '_1.fastq '+ str(SRR)+ '_2.fastq'
    os.system(kallisto_quant)


#define function to create file for input into R sleuth
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


#Step 4

#define function to run bowtie2 on the command line

def bowtie2(SRR):
  #create HCMV bowtie2 index from EF999921 fasta file
  bowtie2_build = "bowtie2-build CDS_EF999921.fasta HCMV"
  os.system(bowtie2_build)
  for i in SRR:
    run_bowtie2 = "bowtie2 --no-unal --al-conc " + SRR + " --quiet -x HCMV -1 "+ SRR +".1.fastq -2 " + SRR + ".2.fastq -S " + SRR + "map.sam"
    os.system(run_bowtie2)

#define function to write reads before and after bowtie2 filtering

def write_reads(SRR):
  log_file = open("miniProject.log", "a") #add output to log file
  for i in SRR:
    name1 = "mv " + SRR + ".1 " + SRR + ".01.fastq"
    os.system(name1) #rename bowtie2 output files to end in fastq
    name2 = "mv " + SRR + ".2 " + SRR + ".02.fastq"
    os.system(name2) #rename bowtie2 output files to end in fastq
    fastq_before = open(SRR + ".1.fastq")
    fastq_after = open(SRR + ".01.fastq")
    count_before = 0 #number of reads in each transcriptome before and after Bowtie2 mapping
    count_after = 0
    for line in fastq_before: #count reads before and after
      count_before = count_before+1
    for line in fastq_after:
      count_after = count_after+1
    count_before = int(count_before/4)
    count_after = int(count_after/4)
    if i == SRR[0]:
      log_file.write("Donor 1 (2dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif i == SRR[1]:
      log_file.write("Donor 1 (6dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif i == SRR[2]:
      log_file.write("Donor 3 (2dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif i == SRR[3]:
      log_file.write("Donor 3 (6dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
  log_file.close() #close miniProject.log file



#Step 5

#define function to run SPAdes

def spades(SRR):
  run_spades = "spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 " + SRR[0] + ".01.fastq --pe1-2 " + SRR[0] + ".02.fastq --pe2-1 " + SRR[1] + ".01.fastq --pe2-2 " + SRR[1] + ".02.fastq --pe3-1 " + SRR[2] + ".01.fastq --pe3-2 " + SRR[2] + ".02.fastq --pe4-1 " + SRR[3] + ".01.fastq --pe4-2 " + SRR[3] + ".02.fastq -o spades" 
  os.system(run_spades) #run on command line
  log_file = open("miniProject.log", "a") #add output to log file
  log_file.write(run_spades + "\n") #write to log file
  log_file.close() #close file


#Step 6
#define function to count the number of contigs
def num_contigs():
  num = "python3 miniproject6.py"
  os.system(num)

#Step 7
#define function to calculate how many base pairs are in the contigs
def bp_contigs():
  length = "python3 miniproject7.py"
  os.system(length)











#Lastly, run functions by calling them
fasta_file()
use_Kallisto(SRR)
sleuth_input(SRR)
#sleuth()
bowtie2(SRR)
#write_reads(SRR)
spades(SRR)
num_contigs()
bp_contigs()
