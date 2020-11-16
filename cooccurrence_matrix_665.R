
rm(list = ls())

##### Regulation expression

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape)
library(splitstackshape)


load("/home/ama42bf/Dokumente/Data_files/727_data/727_Expression_data.Rda") ### d<-Expression data matrix

allele_count_filter<-subset(ara_stop_main_datatable, AC < 2)
merg_data<-allele_count_filter[,c(1,2,11)]

######### Filter STOP GAIN genes only

gene_data<-allele_count_filter[,1:2]
gene_data[,3]<-allele_count_filter$IDs
colnames(gene_data)[3]<-c('accession_IDs')
gene_list<-unique(gene_data$Gene)
d_<-t(d)
e<-as.matrix(subset(d_, rownames(d_) %in% gene_list))


######## GET GENE DATA and change column name accordingly


#A<-matrix(unlist(strsplit(colnames(e),split='X')),ncol=2,byrow=T)[,1]
#colnames(e)<-A
#x_<-as.data.frame(t(x))
#colnames(x_)<-(x)[,1]
#x_<-x_[-c(1),]
#rownames(x)<-x[,1]
#x$gene_id<-NULL
sel_genes<-(rownames(e))
gene_fil<-filter(gene_data, Gene %in% sel_genes)
#D<-as.data.frame(lapply(x_, as.numeric))
d<-as.matrix(t(e))

##### make subsets of data 


a<-separate_rows(gene_fil, accession_IDs, sep = ",")
acc<-colnames(e)
gene_data_fil<-filter(a, accession_IDs %in% acc)

###### SNP_wise

SW<-ddply(gene_data_fil, "SNP", summarize, accession_IDs = paste(accession_IDs, collapse = ","))
stop_<-allele_count_filter[,1:2]
merg_data<-merge(SW, stop_, by ='SNP')

Tim <- Sys.time()
for (i in 1:nrow(merg_data)){
  merg_data[i,4]<-mean(subset(d,rownames(d)%in% unlist(strsplit(as.character(merg_data[i,2]),split=',')),which( colnames(d)%in% merg_data[i,3] )))
  cat(i)
}
Sys.time()-Tim

Tim <- Sys.time()
for (i in 1:nrow(merg_data)){
  merg_data[i,5]<-mean(subset(d,!rownames(d)%in% unlist(strsplit(as.character(merg_data[i,2]),split=',')),which( colnames(d)%in% merg_data[i,3] )))
  cat(i)
}
Sys.time()-Tim


Tim <- Sys.time()
b<-na.omit(merg_data)
for( i in 1:nrow(b)) {
  if(length(unique(unlist(strsplit(as.character(b[i,2]),split=','))))>1){
    a1<-subset(d,rownames(d)%in% unlist(strsplit(as.character(b[i,2]),split=',')))[,colnames(d)%in% b[i,3] ]
    b1<-subset(d,!rownames(d)%in% unlist(strsplit(as.character(b[i,2]),split=',')))[,colnames(d)%in% b[i,3] ]
    
    b[i,6]<-t.test(a1,b1,var.equal = TRUE,paired=FALSE)$p.value
    
  } else {next}
  cat(i)
}
Sys.time()-Tim
colnames(b)[4:6]<-c('mean_KO', 'mean_WT','p-value')


####### filter data

b_<-na.omit(b)

Stop_ara_T1<-subset(b_,(b_$`p-value` <  5*10^(-2)))

##################################################################

genefile<-Stop_ara_T1
gene_data<-as.data.frame(matrix(nrow=nrow(genefile),ncol=2)) 
colnames(gene_data)=c('Gene','IDs') 

gene_data$Gene<-genefile$Gene
gene_data$IDs<-genefile$accession_IDs

###calculate accession centric table
tab<-ddply(gene_data, "Gene", summarize, IDs = paste(IDs, collapse = ","))
df3 <- cSplit(tab, "IDs", sep = ",", direction = "long")
pre_idcount<-unique(df3)
ID_Count<-data.frame ( table ( pre_idcount$IDs ) )
colnames(ID_Count)<-c('IDs','no_of_stop')
e<-ddply(pre_idcount,.(IDs),summarise,  Gene = paste(unique(Gene),collapse = ','))
df4<-ddply(pre_idcount,.(Gene),summarise,  IDs = paste(unique(IDs),collapse = ','))
accession_centric<-merge(ID_Count,e,by="IDs")

### calculate genewise table
gene_count<-data.frame ( table ( pre_idcount$Gene ) )
genewise<-gene_data
colnames(genewise)<-c('Gene','IDs')
gw<-cSplit(genewise, "IDs", sep = ",", direction = "long")
genewise<-ddply(gw,.(Gene) ,summarise,  IDs = paste(unique(IDs),collapse = ','))


### co-occurence 
d_<-genewise
colnames(d_)<-c('Gene', 'Accession')
c<-d_$Gene
Tim <- Sys.time()
A <- matrix(NA, nrow = length(unique(genefile$Gene)),ncol = length(unique(genefile$Gene)))  
colnames(A)<-c
for (s in 1:length(unique(genefile$Gene))) {
  for (k in 1:length(unique(genefile$Gene))){
    A[s,k] <- length(which(strsplit(as.character(d_$Accession[s]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
  }
  cat(s)
}
Sys.time()-Tim

### network table
rownames(A)<-colnames(A)
ind <- which( upper.tri(A,diag=F) , arr.ind = TRUE )
Q<- data.frame( col = dimnames(A)[[2]][ind[,2]] ,
                row = dimnames(A)[[1]][ind[,1]] ,
                val = A[ind] )
Q<-Q[,c(2,3,1)]
f<-Q
f[,4]<-f$val
f[,2]<-f$col
f[,3]<-f$V4
f[,4]<-NULL
colnames(f)<-c('Gene1', 'Gene2', 'connection')
colnames(gene_count)[1]<-'Gene1'
f1<-merge(f, gene_count, by= 'Gene1')
colnames(gene_count)[1]<-'Gene2'
f2<-merge(f1, gene_count, by= 'Gene2')
colnames(f2)[4:5]<-c('KO_gene1', 'KO_gene2')

###testing co-occurrence
which(f2$connection > f2$KO_gene2)
which(f2$connection > f2$KO_gene1)

### calculate significance level (p-value) of all connections
#connection<-f2$connection
#gene1_count<-f2$KO_gene1
#gene2_count<-f2$KO_gene2

#f2$p_value_over<-phyper(connection,gene1_count,665-gene1_count,gene2_count)
#f2$p_value_under<-phyper(connection,gene1_count,665-gene1_count,gene2_count,lower.tail=FALSE)

#sig_over<-subset(f2,f2$p_value_over <  (5*10^(-2)/nrow(f2)))
#sig_under<-subset(f2,f2$p_value_under <  (5*10^(-2)/nrow(f2)))

####calculate distance between genes

#load("~/Dokumente/Data_files/tair10.rda")

genelist2<-unique(f2$Gene2)
sub_ara11<-subset(ara11, ara11$Gene %in% genelist2)

colnames(sub_ara11)[5]<-'Gene2'
f3<-merge(f2, sub_ara11, by = 'Gene2')

colnames(f3)<-c("Gene2", "Gene1" , "connection","KO_gene2" ,"KO_gene1","Chr_2" ,  "Start_2", "Stop_2", "orientation_2")

colnames(sub_ara11)[5]<-'Gene1'
f4<-merge(f3, sub_ara11, by = 'Gene1')

library(matrixStats)

colnames(f4)[14]<-'distance'
f4$distance <- (rowMaxs(as.matrix(f4[c(7,11)]))) - (rowMins(as.matrix(f4[c(8,12)])))

f5<-filter(f4, (f4$distance > 100000) & (f4$Chr == f4$Chr_2))
f6<-filter(f4, f4$Chr != f4$Chr_2)

f7<-rbind(f5,f6)

#length(which(f7$p_value_over <  (5*10^(-2)/nrow(f7))))
#length(which(f7$p_value_under <  (5*10^(-2)/nrow(f7))))

#sig_over_100kb<-subset(f7, f7$p_value_over <  (5*10^(-2)/nrow(f7)))
#sig_under_100kb<-subset(f7, f7$p_value_under <  (5*10^(-2)/nrow(f7)))

Tim<-Sys.time()
thres<-0.05/nrow(f7)
for (i in 1:nrow(f7)){
f7[i,15]<- qhyper(0.05,f7[i,4],664-f7[i,4],f7[i,5])
f7[i,16]<- qhyper(0.05,f7[i,4],664-f7[i,4],f7[i,5],lower.tail=FALSE)
f7[i,17]<- qhyper(thres,f7[i,4],664-f7[i,4],f7[i,5])
f7[i,18]<- qhyper(thres,f7[i,4],664-f7[i,4],f7[i,5],lower.tail=FALSE)
f7[i,19]<- dhyper(f7[i,3],f7[i,4],664-f7[i,4],f7[i,5])
cat(i)
}
Sys.time()-Tim

f8<-subset(f7, V21<0.05/nrow(f7))
colnames(f8)[c(17:21)]<-c('low_lim_0.05', 'up_lim_0.05', 'low_lim_thres', 'up_lim_thres', 'Pvalue')
f8[,22]<-NA
colnames(f8)[22]<-'sig'
d1<-subset(f8, f8[,3]<=f8[,19])
d1$sig<-'U'
d2<-subset(f8, f8[,3]>=f8[,20])
d2$sig <-'O'
f9<-rbind(d1,d2)
f10<-f9[,c(1:5,17:22)]
f11<- cbind(f10, sample(22222:99999, nrow(f10), replace = FALSE))
f12<-f11[,c(12,4:5)]
f13<-f12[apply(f12[,-1], MARGIN = 1, function(x)all(x>=100)),] ### extract genes knocked out in more than or equal to 100 accessions
colnames(f13)[1]<-'sam'
colnames(f11)[12]<-'sam'
f13$KO_gene1<-NULL
f13$KO_gene2<-NULL
f14<- merge(f13, f11 ,by = 'sam')


#### add accessions for both genes

colnames(d_)<-c('Gene1','G1_accession_IDs')
f15<-merge(f14,d_, by = 'Gene1')
d_2<-d_
colnames(d_2)<-c('Gene2','G2_accession_IDs')
f16<-merge(f15,d_2, by = 'Gene2')

#### find common accessions from KO_Gene1 and KO_Gene2 for GWAS

tim<-Sys.time()
for (i in 1:nrow(f16)) {
  k<-unlist(strsplit(f16[i,13],split = ','))
  j<-unlist(strsplit(f16[i,14], split = ','))
  l<-paste(intersect(k,j),collapse=',')
  f16[i,15]<-l
  cat(i)
}
Sys.time()-tim
S<-f16
for (i in 1:nrow(S)) {
       S[i,9]<-length(unlist(strsplit(S[i,5],split = ',')))
       }


#save(f7, file = paste('syn_665_T1_down.Rda'))
save.image(file = 'syn_665_T1_down.RData')

