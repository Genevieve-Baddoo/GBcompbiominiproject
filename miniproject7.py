#import biopython 
from Bio.Seq import Seq
from Bio import SeqIO

#open files needed
contigs = open("contigs_1000.fasta")
log_file = open("miniProject.log", "a")

num_bps = 0 #initialize count of all base pairs

#pull out sequences out of file
for contig in SeqIO.parse(contigs, "fasta"):
  #retrieve length and add to num_ bps
  length = len(contig.seq)
  num_bps = num_bps + length

#update log file
log_file.write("There are " + str(num_bps) + " bp in the assembly.\n")
log_file.close()
