
rm(list = ls())
#setwd(home/ama42bf/Dokumente/Data_files/STOP_data/STOP_665)
load('~/Dokumente/Data_files/STOP_data/STOP_665/Pre_GWAS_file.Rda')
f17<-Pre_GWAS_file
f17$com_g1g2 = paste(f17$Gene1, f17$Gene2, sep="_")
f17$com_g2g1 = paste(f17$Gene2, f17$Gene1, sep="_")
f18<-f17[,c(11,3,4,8,10)]
f19<-f17[,c(12,3,5,9,10)]

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape)
library(splitstackshape)
#####Gene1
tim<-Sys.time()
for (i in 1:nrow(f18)){
gwas_acc_table<-as.data.frame(matrix(nrow=f18[i,3],ncol=2))

colnames(gwas_acc_table)<-c('KO_ids', 'com_ids')

gwas_acc_table$KO_ids<-unlist(strsplit(f18[i,4],split=','))
r<-unlist(strsplit(f18[i,5],split=','))
gwas_acc_table$com_ids<-rep(NA)

  for (j in 1:nrow(gwas_acc_table)) {
   if(gwas_acc_table[j,1] %in% r){
   gwas_acc_table[j,2]<-paste(1)
  } else {gwas_acc_table[j,2]<-paste(0)}
}
if(sum(as.numeric(gwas_acc_table[,2])) == f18[i,2]){print(TRUE)
  
  print('Success----file being saved!')

  assign  (paste("Y",f18[i,1],sep="_"),gwas_acc_table)

  save(gwas_acc_table,file = paste(f18[i,1],'.Rda'))

  }else {print('Error!!! file not saved.')}
cat(i)
}
Sys.time()-tim

#####Gene2##########
tim<-Sys.time()
for (i in 1:nrow(f19)){
  gwas_acc_table<-as.data.frame(matrix(nrow=f19[i,3],ncol=2))
  
  colnames(gwas_acc_table)<-c('KO_ids', 'com_ids')
  
  gwas_acc_table$KO_ids<-unlist(strsplit(f19[i,4],split=','))
  r<-unlist(strsplit(f19[i,5],split=','))
  gwas_acc_table$com_ids<-rep(NA)
  
  for (j in 1:nrow(gwas_acc_table)) {
    if(gwas_acc_table[j,1] %in% r){
      gwas_acc_table[j,2]<-paste(1)
    } else {gwas_acc_table[j,2]<-paste(0)}
  }
  if(sum(as.numeric(gwas_acc_table[,2])) == f18[i,2]){print(TRUE)
    print('Success----file being saved!')
    assign  (paste("Y",f19[i,1],sep="_"),gwas_acc_table)
    save(gwas_acc_table,file = paste(f19[i,1],'.Rda'))
  }else {print('Error!!! file not saved.')}
  cat(i)
}
Sys.time()-tim
