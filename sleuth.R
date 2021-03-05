#load sleuth package 
library(sleuth)
#load the dplyr package for data frame filtering 
library(dplyr)
#load data.table package to read in table data  
library(data.table)

#read in the table describing samples and kallisto output 
stab<-read.table("sleuth_infile.txt",header=TRUE,stringsAsFactors=FALSE,sep='\t')

#initialize sleuth object
so<-sleuth_prep(stab)


#fit a model comparing two conditions 
so<-sleuth_fit(so,~condition,'full')

#fit the reduced model to compare in the likelihood ratio test
so<-sleuth_fit(so,~1,'reduced')

#perform the likelihood ratio test for differential expression between conditions (2dpi & 6dpi)
so<-sleuth_lrt(so,'reduced','full')

#extract the test results from the sleuth object 
sleuth_table<-sleuth_results (so,'reduced:full','lrt',show_all=FALSE)

#filter most significant results (FDR/qval <0.05) and sort by pval 
sleuth_significant<-dplyr::filter(sleuth_table,qval<=0.05)%>%dplyr::arrange(pval)
sleuth_output<-sleuth_significant %>% select(target_id,test_stat,pval,qval)

#write FDR <0.05 transcripts to file 
#sleuth_output<-na.omit(sleuth_output)
write.table(sleuth_output,file='sleuth_outfile.txt', quote=FALSE,row.names=FALSE)
