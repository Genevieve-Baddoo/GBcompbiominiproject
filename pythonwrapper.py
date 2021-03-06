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


#Step 4

#function to run bowtie2
def bowtie2(SRRs):
  #build bowtie2 index from EF999921 fasta file
  run_build = "bowtie2-build EF999921_CDS.fasta HCMV"
  os.system(run_build)
  for SRR in SRRs:
    run_bowtie2 = "bowtie2 --no-unal --al-conc " + SRR + " --quiet -x HCMV -1 "+ SRR +".1_1.fastq -2 " + SRR + ".1_2.fastq -S " + SRR + "map.sam"
    os.system(run_bowtie2)

#find the number of reads in each transcriptome before and after the Bowtie2 mapping
def num_reads(SRRs):
    for i in range(0,len(SRRs)):
        if i == 0:
            donor = "Donor 1 (2dpi)"
        elif i == 1:
            donor = "Donor 1 (6dpi)"
        elif i == 2:
            donor = "Donor 3 (2dpi)"
        elif i ==3:
            donor = "Donor 3 (6dpi)"

        beforebow_1 = open(str(SRRs[i])+'.1_1.fastq').readlines()
        beforebow_2 = open(str(SRRs[i])+'.1_2.fastq').readlines()
        before_reads = (len(beforebow_1) + len(beforebow_2))/8

        afterbow_1 = open('Bowtie2_'+ str(SRRs[i])+'.1.fastq').readlines()
        afterbow_2 = open('Bowtie2_'+ str(SRRs[i])+'.2.fastq').readlines()
        after_reads = (len(afterbow_1) + len(afterbow_2))/8


        log_file = open("miniProject.log", "a")
        log_file.write(donor + "had " + str(before_reads) + " read pairs before Bowtie2 filtering and " + str(after_reads) + " read pairs after." + "\n")
        log_file.close()


#Step 5
#define function to run spades


def run_spades(SRRs):
    spades = 'spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 ' + SRRs[0] + '.1.fastq --pe1-2 '+ SRRs[0] + '.2.fastq --pe2-1 ' + SRRs[1] + '.1.fastq --pe2-2 ' + SRRs[1] + '.2.fastq --pe3-1 ' + SRRs[2] + '.1.fastq --pe3-2 ' + SRRs[2] +'.2.fastq --pe4-1 ' + SRRs[3] + '.1.fastq --pe4-2 ' + SRRs[3] + '.2.fastq -o spades'
    os.system(spades)
    miniProject_log = open("miniProject.log", "a")
    miniProject_log.write(spades + '\n')
    miniProject_log.close()


#Step 6

#function to conunt and subset contigs longer than 1000 bp  
def contig_subset():
  subset = "python3 miniproject6.py" #python code in repo
  os.system(subset)


#Step 7

#function to calculate assembly length  
def length_assembly():
  calc = "python3 miniproject7.py" #python code in repo
  os.system(calc)


#Step 8
  
#function to run blast using the largest contig  
#make sure to blast against a local database
def blast_longestcontigs():
    blast_input = '/sequence.fasta'
    makeblastdb_command = 'makeblastdb -in ' + path + blast_input +' -out ' + path + '/betaherpesvirinae -title betaherpesvirinae -dbtype nucl'
    print(path + input_blast)
    os.system(makeblastdb_command)

    blastn_command = 'blastn -query ' + path + '/longest_contig.fasta -db '+ path + 'betaherpesvirinae -max_target_seqs 10 -out ' + path + '/blast_results.txt -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"'
    os.system(blastn_command)
    log_file.write('sacc' + '\t' + 'pident' + '\t' + 'length' + '\t' + 'qstart' + '\t' + 'qend' + '\t' + 'sstart' + '\t' + 'send' + '\t' + 'bitscore' + '\t' + 'eval' + '\t' + 'stitle' + '\n') #write to log
    read_blast_results = open(path + 'blast_results.txt', 'r').readlines()
    for i in read_blast_results:
        log_file.write(str(i))
        log_file.write('\n')





















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
contig_subset()
length_assembly()
blast_longestcontigs()
