
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape)
library(devtools)
library(splitstackshape)

#################################### SAMPLE DATA GENERATION T1 ###################
load("~/Dokumente/Data_files/STOP_data/STOP_665/Stop_data_for_cooccur.Rda") ### named S
r<-S$count

load("~/Dokumente/Data_files/SYN_665/SYN_data_for_cooccur.Rda") ### named S_
colnames(S_)[3]<-'count'
load("~/Dokumente/Data_files/SYN_665_afterGT_analysis_updated.Rda") ### named syn_data4
############## ******* SAMPLE DATASET******** ############

simul_results<-as.data.frame(matrix(nrow=100,ncol=2))
colnames(simul_results)<-c('Under_rep','Over_rep')
for(i in 1:100 ){
  
no_of_rows <- 2545
no_of_unique_gene <- 1671
temp <- S_

while(n_distinct(temp$Gene) != no_of_unique_gene) {
  gene <- sample(unique(S_$Gene),no_of_unique_gene)
  temp <- S_[S_$count %in% r & S_$Gene %in% gene, ]
}
part1  <- temp %>% group_by(Gene) %>% sample_n(floor(no_of_rows/no_of_unique_gene))
part2 <- temp %>% anti_join(part1) %>% sample_n(no_of_rows - nrow(part1)) 
genefile <- bind_rows(part1, part2)
snp_list<-genefile$SNP
genefile<-subset(syn_data4, SNP %in% snp_list)

####################### ********MAIN FUNCTION********** ###########

gene_data<-as.data.frame(matrix(nrow=nrow(genefile),ncol=2)) 
colnames(gene_data)=c('Gene','IDs') 

gene_data$Gene<-genefile$Gene
gene_data$IDs<-genefile$IDs

### calculate accession centric table #####
tab<-ddply(gene_data, "Gene", summarize, IDs = paste(IDs, collapse = ","))
df3 <- cSplit(tab, "IDs", sep = ",", direction = "long")
pre_idcount<-unique(df3)
ID_Count<-data.frame ( table ( pre_idcount$IDs ) )
colnames(ID_Count)<-c('IDs','no_of_stop')
e<-ddply(pre_idcount,.(IDs),summarise,  Gene = paste(unique(Gene),collapse = ','))
df4<-ddply(pre_idcount,.(Gene),summarise,  IDs = paste(unique(IDs),collapse = ','))
accession_centric<-merge(ID_Count,e,by="IDs")
View(accession_centric)

### calculate genewise table ####
gene_count<-data.frame ( table ( pre_idcount$Gene ) )
genewise<-genefile[,c(3,2)]
colnames(genewise)<-c('Gene','IDs')
gw<-cSplit(genewise, "IDs", sep = ",", direction = "long")
gene_id_count<-data.frame ( table ( gw$Gene ) )
genewise<-ddply(gw,.(Gene) ,summarise,  IDs = paste(unique(IDs),collapse = ','))
View(genewise)
### co-occurence #####

d_<-genewise
colnames(d_)<-c('Gene', 'Accession')
c<-d_$Gene
Time <- Sys.time()
A <- matrix(NA, nrow = length(unique(genefile$Gene)),ncol = length(unique(genefile$Gene)))  
colnames(A)<-c
for (s in 1:length(unique(genefile$Gene))) {
  for (k in 1:length(unique(genefile$Gene))){
    A[s,k] <- length(which(strsplit(as.character(d_$Accession[s]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
  }
  cat(s)
}
Sys.time()-Time

### network table #####
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

### testing co-occurrence
which(f2$connection > f2$KO_gene2)
which(f2$connection > f2$KO_gene1)

####calculate distance between genes

load("~/Dokumente/Data_files/tair10.rda")

genelist2<-unique(f2$Gene2)
sub_tairT1down<-subset(tair10, tair10$Gene %in% genelist2)

colnames(sub_tairT1down)[5]<-'Gene2'
f3<-merge(f2, sub_tairT1down, by = 'Gene2')

colnames(f3)<-c("Gene2", "Gene1" , "connection","KO_gene1" ,"KO_gene2" ,"Chr_2" ,  "Start_2", "Stop_2", "orientation")

colnames(sub_tairT1down)[5]<-'Gene1'
f4<-merge(f3, sub_tairT1down, by = 'Gene1')

library(matrixStats)
f4$distance<-rep(NA)
f4$distance <- (rowMaxs(as.matrix(f4[c(7,11)]))) - (rowMins(as.matrix(f4[c(8,12)])))

f5<-filter(f4, (f4$distance > 100000) & (f4$Chr == f4$Chr_2))
f6<-filter(f4, f4$Chr != f4$Chr_2)

f7<-rbind(f5,f6)

Time<-Sys.time()
thres<-0.05/nrow(f7)
for (i in 1:nrow(f7)){
  f7[i,15]<- qhyper(0.05,f7[i,4],665-f7[i,4],f7[i,5])
  f7[i,16]<- qhyper(0.05,f7[i,4],665-f7[i,4],f7[i,5],lower.tail=FALSE)
  f7[i,17]<- qhyper(thres,f7[i,4],665-f7[i,4],f7[i,5])
  f7[i,18]<- qhyper(thres,f7[i,4],665-f7[i,4],f7[i,5],lower.tail=FALSE)
  f7[i,19]<- dhyper(f7[i,3],f7[i,4],665-f7[i,4],f7[i,5])
  cat(i,'_')
}
Sys.time()-Time
###### Over under representation #####
f7[,20]<-NA
colnames(f7)[20]<-'sig'
d1<-subset(f7, f7[,3]<f7[,17])
d1$sig<-'U'
d2<-subset(f7, f7[,3]>f7[,18])
d2$sig <-'O'
f8<-rbind(d1,d2)
colnames(f8)[c(15:20)]<-c('lower_limit_0.05', 'upper_limit_0.05', 'lower_limit_thres', 'upper_limit_thres', 'Pvalue','sig')

simul_results[i,1]<-length(which(f8$sig == 'U'))
simul_results[i,2]<-length(which(f8$sig == 'O'))       
}

save(simul_results, file = paste('/home/ama42bf/Dokumente/Data_files/syn_665_T1_down_simul.Rda'))
#save.image(file = 'syn_665_T1_down.RData')


#########################################################******************#########################################################################