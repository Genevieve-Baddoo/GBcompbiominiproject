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
        wget = 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/' + SRR + '/' + SRR + '.1'
        rename = 'mv ' + SRR + '.1 ' + SRR
        fastq_dump = 'fastq-dump -I --split-files ' + SRR + '.1'
        os.system(wget)
        os.system(rename)
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


def sleuth():
    sleuth_command = 'Rscript sleuth.R'
    os.system(sleuth_command)
    log_file = open("miniProject.log", "w")
    read_sleuth = open('sleuth_outfile.txt', 'r').readlines()
    for SRR in read_sleuth:
        log_file.write(SRR + '\n')



#function to run bowtie2
def bowtie2(SRRs):
  #build bowtie2 index from EF999921 fasta file
  run_build = "bowtie2-build EF999921_CDS.fasta HCMV"
  os.system(run_build)
  for SRR in SRRs:
    run_bowtie2 = "bowtie2 --no-unal --al-conc " + SRR + " --quiet -x HCMV -1 "+ SRR +".1_1.fastq -2 " + SRR + ".1_2.fastq -S " + SRR + "map.sam"
    os.system(run_bowtie2)


#define function to count reads before and after filtering with bowtie2
def num_reads(SRRs):
  log_file = open("miniProject.log", "a")
  for SRR in SRRs:
    name1 = "mv " + SRR + ".1_1." + SRR + ".1.fastq"
    os.system(name1) #rename bowtie2 output files to end in fastq
    name2 = "mv " + SRR + ".1_2." + SRR + ".2.fastq"
    os.system(name2) #rename bowtie2 output files to end in fastq
    fastq_before = open(SRR + ".1_1.fastq")
    fastq_after = open(SRR + ".1.fastq")
    count_before = 0
    count_after = 0
    for line in fastq_before:
      count_before = count_before+1
    for line in fastq_after:
      count_after = count_after+1
    count_before = int(count_before/4)
    count_after = int(count_after/4)
    if SRR == SRR[0]:
      log_file.write("Donor 1 (2dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif SRR == SRR[1]:
      log_file.write("Donor 1 (6dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif SRR == SRR[2]:
      log_file.write("Donor 3 (2dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")
    elif SRR == SRR[3]:
      log_file.write("Donor 3 (6dpi) had " + str(count_before) + " read pairs before Bowtie2 filtering and " + str(count_after) + " read pairs after.\n")    
  log_file.close()

def run_spades(SRRs):
    spades = 'spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 ' + SRRs[0] + '.1.fastq --pe1-2 '+ SRRs[0] + '.2.fastq --pe2-1 ' + SRRs[1] + '.1.fastq --pe2-2 ' + SRRs[1] + '.2.fastq --pe3-1 ' + SRRs[2] + '.1.fastq --pe3-2 ' + SRRs[2] +'.2.fastq --pe4-1 ' + SRRs[3] + '.1.fastq --pe4-2 ' + SRRs[3] + '.2.fastq -o spades'
    os.system(spades)
    miniProject_log = open("miniProject.log", "a")
    miniProject_log.write(spades + '\n')
    miniProject_log.close()























#run functions
fastq(SRRs)
fasta()
kallisto_idx()
kallisto_quant(SRRs)
sleuth_input(SRRs)
sleuth()
bowtie2(SRRs)
num_reads(SRRs)
run_spades(SRRs)
