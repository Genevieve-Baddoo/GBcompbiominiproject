
HCMV (Human cytomegalovirus) Mini Project.
This respisotory contains a python pipeline that takes SRA input reads from the HCMV virus to produce an assembly which is blasted and analyzed. 
It also contains a test dataset which is a subset of the input fastq files from the fastq function. The SRR IDs for the test dataset are: SRR5660030 SRR5660033 SRR5660044 SRR5660045


Software Tools Required
=========

* Python3
* Biopython
* Kallisto
* R
* SPAdes
* Bowtie2


Main Python Script
==================

pythonwrapper.py
----------------

To run this repo, please clone this github to your current working directory:

`https://github.com/Genevieve-Baddoo/hcmvminiproject.git`

Change directory:
`cd hcmvminiproject`




**To run this script:**







**Project Steps**


Step 1. We would like to compare HCMV transcriptomes 2- and 6-days post-infection (dpi). First, retrieve the following transcriptomes from two patient donors from SRA and convert to paired-end fastq files. You can use wget (by constructing the path based on the SRR numbers for each of these samples).
* [Donor 1 (2dpi)](https://www.ncbi.nlm.nih.gov/sra/SRX2896360)
*
*
*



