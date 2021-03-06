#import biopython 
from Bio.Seq import Seq
from Bio import SeqIO

#open files
contigs = open("spades/contigs.fasta") #created after SPAdes run 
log_file = open("miniProject.log", "a")
output = open("contigs_1000.fasta", "w")
count = 0 #initialize count

#pull sequences out of contig
for contig in SeqIO.parse(contig, "fasta"):
  if len(contig.seq) > 1000:
    #count contigs that meet threshold
    count = count + 1 
    #output contigs that meet threshold
    output.write(">" + str(contig.id) + "\n")
    output.write(str(contig.seq) + "\n")

#update log with the number of contigs
log_file.write("There are " + str(count) + " contigs > 1000 bp in the assembly.\n")
log_file.close() 
output.close()

