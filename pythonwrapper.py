#import necessary modules for HCMV project
import os
import argparse
from Bio import SeqIO
from Bio import Entrez
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW

'''

1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi).
retrieve the following transcriptomes from two patient donors from SRA and convert to paired-end fastq files. Use wget.
Donor 1 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896360
Donor 1 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896363
Donor 3 (2dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896374
Donor 3 (6dpi): https://www.ncbi.nlm.nih.gov/sra/SRX2896375

'''


#can define a function for each tool. NEED TO COMMENT OUT SECTIONS WHEN DONE
#will need to add test data (small # of input reads

#convert SRA files from NCBI to fastq paired-end reads. Construct path based on SRR numbers

#create list with acession numbers listed
SRR = ['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']

def fastq(SRR):
  wget = 'wget https://sra-download.st-va.ncbi.nlm.nih.gov/sos1/sra-pub-run-12/SRR' + SRR + '/' + SRR + '.1'
  fastq = 'fastq-dump -I --split-files '+ str(SRR)
  rename = 'mv '+ SRR + '.1 ' + SRR
  os.system(wget) #os.system runs the command line
  os.system(fastq)
  os.system(rename)

for i in SRR:
  fastq(i)

'''

2. We will quantify TPM in each sample using kallisto, but first we need to build a transcriptome index for HCMV (NCBI accession EF999921)
Use Biopython to retrieve and generate the appropriate input and then build the index with kallisto.
Need to extract the CDS features from the GenBank format. Write the following to your log file (replace # with the number of coding sequences in the HCMV genome):
The HCMV genome (EF999921) has # CDS.

wrote code on a seperate py file!

'''

#define function to create fasta file from CDS of NCBI Accession EF999921
def transcriptome_index():
  getfasta = "/usr/bin/python3 miniproject2.py" 
  os.system(getfasta)

