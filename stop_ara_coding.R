
#########################################################################
#### extract the knockout accession IDs from big vcf file through the following code:
#
## extract head first from the large snpeff vcf file
#
#   head -n 12 1001genomes_snp-short-indel_only_ACGTN_v3.1.vcf.snpeff > newfile1.txt
#
## to delete unused columns
#
#   cut -f 3,4,5,6,7,8,9 --complement stop_ID.txt > syn_ids.txt
#
## to merge both header and filtered vcf file
#   cat syn_id.txt syn_ids.txt > c.txt 
#   
## save both files differently then merge

#   paste -d"\t\t" pos.txt ids.txt > c.txt
#   
#   awk '{print $1"- "$2}' test_del.txt > pos.txt
#      awk -F'\t' '
#      NR==1 {for(i=1;i<=NF;i++) h[i]=$i; next}
#      {
#           sep="";
#           for(x=1;x<=NF;x++) {
#           if($x ~ /1\|1/) {
#      printf "%s%s", sep, h[x], 
#      sep=", ";
#             }
#         }
#         print ""
#   }' c.txt > ids.txt
##########################################################################


rm(list = ls())

load("/home/ama42bf/Source_data_files/Ann_gene.rda")
load("/home/ama42bf/Source_data_files/ara11.rda")

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape)
library(splitstackshape)
library(stringr)
setwd("~/Dokumente/Data_files/STOP_data/stop_ara")


#********+++++++++++++++++++++++++ NEW START HERE !!!!!! +++++++++++++++++++++++++++++
#########################################################################

#extract gene names from vcf file
stop<-read.delim('1001G_ara11_stop.vcf', sep = '', header= FALSE)

anno_Stop<-stop[,c(1,2,4,5,8)]
colnames(anno_Stop)<-c('Chr','POS','Ref','Alt','snpeff')
anno_Stop<-unite(anno_Stop, SNP, c(Chr, POS), sep = '-', remove=FALSE)
ara_stop<-read.delim('stop_final.txt', sep = '', header = FALSE)
colnames(ara_stop)<-c('SNP', 'IDs')
ann<-as.data.frame(anno_Stop[,c(1,6)])
test<-as.data.frame(str_split_fixed(ann$snpeff,'\\|', 6)) ### split column where there is a pipe
ara_stop$Gene<-test$V5
stop_genelist<-ara_stop$Gene
fil1 <- filter(ara11, Gene %in% stop_genelist)
d_<-merge(anno_Stop, ara_stop, by = 'SNP')
d<-merge(fil1, d_, by = 'Gene')
colnames(d)[7]<-'Chr'
d$Chr.x<-NULL

#### calculate the position of SNP and its reference location from start and start position

for (i in 1:nrow(d)) {
  if (d[i,4] == '-'){
    d[i,12]<-d[i,3] - d[i,7]
  }else if (d[i,4] == '+'){
    d[i,12]<-d[i,7] - d[i,2]
  }
}
colnames(d)[12]<-'pos_on_genome'
d$total<-d$Stop-d$Start
d$rel_pos<-d$pos_on_genome/d$total

#d$ref_start_syn = d$POS -d$Start
#d$ref_stop_syn = d$POS -d$Stop
#d$total<-d$Stop-d$Start
#d$rel_start<-d$ref_start_syn/d$total
#d$rel_stop<-d$ref_stop_syn/d$total
#d[,5]<-NULL
#stop_data<-d[,c(2,5)]

###############################################################

### removing all snps with more than 1 alt allele
count_alt_alle<-d[,c(5,9)]
test2 <- cSplit(count_alt_alle, "Alt", sep = ",", direction = "long")
Allele_Count<-as.data.frame(table(test2$SNP ) )
colnames(Allele_Count)<-c('SNP', 'AC')
df<-merge(d, Allele_Count, by = 'SNP')

ara_stop_main_datatable<-df
save(ara_stop_main_datatable, file = paste('ara_stop_main_datatable.Rda'))

##############################################################
### filter snps with more than 1 alt allele
df<-ara_stop_main_datatable
allele_count_filter<-subset(df, AC < 2)

### extracting first 20% snps on the genome 

filter_20per_on_genome<-subset(allele_count_filter, rel_pos <= 0.20)
for (i in (1:4344)) {
  filter_20per_on_genome[i,16]<- length(unlist(strsplit(as.character(filter_20per_on_genome[i,11]),split=',')))
}
filter_20per_on_genome_ex_singleKO<-subset(filter_20per_on_genome, acc_length > 1)
save(filter_20per_on_genome, file = paste('filter_20per_on_genome.Rda'))
save(filter_20per_on_genome_ex_singleKO, file = paste('filter_20per_on_genome_ex_singleKO.Rda'))

###############################################################


#####  calculate accession centric ( which IDs have been knocked out in how many genes)

genefile<-filter_20per_on_genome_ex_singleKO
gene_data<-as.data.frame(matrix(nrow=nrow(genefile),ncol=2)) 
colnames(gene_data)=c('Gene','IDs') 

gene_data$Gene<-genefile$Gene
gene_data$IDs<-genefile$IDs

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
sub_ara11<-subset(ara11, Gene %in% genelist2)

colnames(sub_ara11)[5]<-'Gene2'
f3<-merge(f2, sub_ara11, by = 'Gene2')

colnames(f3)<-c("Gene2", "Gene1" , "connection","KO_gene2" ,"KO_gene1","Chr_2" ,  "Start_2", "Stop_2", "orientation_2")

genelist3<-unique(f3$Gene1)
sub_ara11_a<-subset(ara11, Gene %in% genelist3)

colnames(sub_ara11_a)[5]<-'Gene1'
f4<-merge(f3, sub_ara11_a, by = 'Gene1')

library(matrixStats)
f4$distance<-NA
#colnames(f4)[14]<-'distance'
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
for (i in 307246:1437844){
  f7[i,15]<- qhyper(0.05,f7[i,4],1135-f7[i,4],f7[i,5])
  f7[i,16]<- qhyper(0.05,f7[i,4],1135-f7[i,4],f7[i,5],lower.tail=FALSE)
  f7[i,17]<- qhyper(thres,f7[i,4],1135-f7[i,4],f7[i,5])
  f7[i,18]<- qhyper(thres,f7[i,4],1135-f7[i,4],f7[i,5],lower.tail=FALSE)
  f7[i,19]<- dhyper(f7[i,3],f7[i,4],1135-f7[i,4],f7[i,5])
  cat(i)
}
Sys.time()-Tim

f8<-subset(f7, V19<0.05/nrow(f7))
colnames(f8)[c(15:19)]<-c('low_lim_0.05', 'up_lim_0.05', 'low_lim_thres', 'up_lim_thres', 'Pvalue')
f8[,20]<-NA
colnames(f8)[20]<-'sig'
d1<-subset(f8, f8[,3]<=f8[,17])
d1$sig<-'U'
d2<-subset(f8, f8[,3]>=f8[,18])
d2$sig <-'O'
f9<-rbind(d1,d2)
stop_20per_ara_cooccurence_analysis_complete_table<-f7
save(stop_20per_ara_cooccurence_analysis_complete_table, file = paste('stop_20per_ara_cooccurence_analysis_complete_table.Rda'))
stop_20per_ara_significant_cooccur_table_with_U_and_O<-f9
save(stop_20per_ara_significant_cooccur_table_with_U_and_O, file = paste('stop_20per_ara_significant_cooccur_table_with_U_and_O.Rda'))
#############################

gene_data_20<-as.data.frame(matrix(nrow=2395,ncol=2))
colnames(gene_data_20)=c('Gene','IDs')
gene_data_20$Gene<-filter_20per_on_genome_ex_singleKO$Gene
gene_data_20$IDs<-filter_20per_on_genome_ex_singleKO$IDs
tab_20<-ddply(gene_data_20, "Gene", summarize, IDs = paste(IDs, collapse = ","))
df3_20<- cSplit(tab_20, "IDs", sep = ",", direction = "long")
pre_idcount_20<-unique(df3_20)
ID_Count_20<-data.frame ( table ( pre_idcount_20$IDs ) )
colnames(ID_Count_20)<-c('IDs','no_of_stop')
e_20<-ddply(pre_idcount_20,.(IDs),summarise,  Gene = paste(unique(Gene),collapse = ','))
df4_20<-ddply(pre_idcount_20,.(Gene),summarise,  IDs = paste(unique(IDs),collapse = ','))
accession_centric_stop_ara11_20<-merge(ID_Count_20,e_20,by="IDs")

save(accession_centric_stop_ara11_20, file = paste('accession_centric_stop_ara11_20.Rda'))
gene_count_stop_ara11_20<-data.frame ( table ( pre_idcount_20$Gene ) )
save(gene_count_stop_ara11_20, file = paste('gene_idcount_stop_ara11_20.Rda'))

genewise_20<-gene_data_20
gw<-cSplit(genewise_20, "IDs", sep = ",", direction = "long")
genewise<-ddply(gw,.(Gene) ,summarise,  IDs = paste(unique(IDs),collapse = ',')) ### doing this step to make the table genewise
d_<-genewise
colnames(d_)<-c('Gene', 'Accession')
c<-d_$Gene

Tim <- Sys.time()
A <- matrix(NA, nrow = 1698,ncol = 1698)
colnames(A)<-c
for (s in 1:3514) {
  for (k in 1:3514){
    A[s,k] <- length(which(strsplit(as.character(d_$Accession[s]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
  }
  cat(s,'\n')
}
Sys.time()-Tim
rownames(A)<-colnames(A)
ind <- which( upper.tri(A,diag=T) , arr.ind = TRUE )
Q<- data.frame( col = dimnames(A)[[2]][ind[,2]] ,
                row = dimnames(A)[[1]][ind[,1]] ,
                val = A[ind] )

## Q will be the resulting matrix with 3 columns 
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
colnames(f2)[4:5]<-c('KO_gene2', 'KO_gene1')
which(f2$connection > f1$KO_gene2)
which(f2$connection > f1$KO_gene1)

connection<-f2$connection
gene1_count<-f2$KO_gene1
gene2_count<-f2$KO_gene2
f2$p_value_over<-phyper(connection,gene1_count,665-gene1_count,gene2_count)
f2$p_value_under<-phyper(connection,gene1_count,665-gene1_count,gene2_count,lower.tail=FALSE)

##### Regulation expression

load("/home/ama42bf/Dokumente/Data_files/727_data/727_Expression_data.Rda") ### d<-Expression data matrix
allele_count_filter<-subset(df, AC < 2)
merg_data<-ara_stop_main_datatable[,c(1,2,11)]
Tim <- Sys.time()
for (i in 1:nrow(merg_data)){
  merg_data[i,4]<-mean(subset(d,rownames(d)%in% unlist(strsplit(merg_data[i,3],split=',')),which( colnames(d)%in% merg_data[i,2] )))
  cat(i)
}
Sys.time()-Tim

Tim <- Sys.time()
for (i in 1:nrow(merg_data)){
  merg_data[i,5]<-mean(subset(d,!rownames(d)%in% unlist(strsplit(merg_data[i,3],split=',')),which( colnames(d)%in% merg_data[i,2] )))
  cat(i)
}
Sys.time()-Tim


require(svMisc)
Tim <- Sys.time()
b<-na.omit(merg_data)
for( i in 1:nrow(b)) {
  if(length(unique(unlist(strsplit(b[i,3],split=','))))>1){
    a1<-subset(d,rownames(d)%in% unlist(strsplit(b[i,3],split=',')))[,colnames(d)%in% b[i,2] ]
    b1<-subset(d,!rownames(d)%in% unlist(strsplit(b[i,3],split=',')))[,colnames(d)%in% b[i,2] ]
    
    b[i,6]<-t.test(a1,b1,var.equal = TRUE,paired=FALSE)$p.value
    
  } else {next}
  cat(i)
}
Sys.time()-Tim
colnames(b)[4:6]<-c('mean_KO', 'mean_WT','p-value')


####### CALCULATE Distance

name_of_dataframe$distance <- (rowMaxs(as.matrix(name_of_dataframe[c(Startgene1,Startgene2)]))) - (rowMins(as.matrix(name_of_dataframe[c(Stopgene1,Stopgene2)])))

####### PLOT data

b_<-na.omit(b)

Syn_reg_T1<-subset(b_,(b_$`p-value` <  5*10^(-2)))

Syn_reg_T2<-subset(b_,(b_$`p-value` <  5*10^(-2)/nrow(b_)))

##################################################################

### GO SLim (who are these SNPs?)

##################################################################

load('/home/ama42bf/Dokumente/Data_files/727_data/GO_Molecular_FUNCTIONS.Rda')

#T2_Syn[,8]<- round(T2_syn$mean_KO) - round(T2_syn$mean_WT)
T2_Syn[,8]<-abs(b_$mean_WT)>abs(b_$mean_KO)

colnames(T2_Syn[8])<-'direction'
T2_up<-subset(T2_Syn, direction == 'TRUE')
T2_down<-subset(T2_Syn, direction == 'FALSE')

unique_gene_list<-unique(T2_up$Gene)

d<-filter(t1, locus_name %in% unique_gene_list)
graph_table<-as.data.frame(matrix(nrow=15, ncol=2))
rownames(graph_table)<-c('hydrolase_activity', "kinase_activity" , "transferase_activity" , "other_enzyme_activity" , "transcription_factor_activity", "DNA_or_RNA_binding", "nucleic_acid_binding" , "nucleotide_binding",  "protein_binding" , "receptor_binding_or_activity", "other_binding", "structural_molecule_activity", "transporter_activity", "unknown_molecular_functions", "other_molecular_functions")

graph_table[1,2]<-count(filter(d, GO_slim_term == 'hydrolase activity'))
graph_table[2,2]<-count(filter(d, GO_slim_term == 'kinase activity'))
graph_table[3,2]<-count(filter(d, GO_slim_term == "transferase activity"))
graph_table[4,2]<-count(filter(d, GO_slim_term == "other enzyme activity"))
graph_table[5,2]<-count(filter(d, GO_slim_term == "transcription factor activity"))
graph_table[6,2]<-count(filter(d, GO_slim_term == "DNA or RNA binding"))
graph_table[7,2]<-count(filter(d, GO_slim_term == "nucleic acid binding"))
graph_table[8,2]<-count(filter(d, GO_slim_term == "nucleotide binding"))
graph_table[9,2]<-count(filter(d, GO_slim_term == "protein binding"))
graph_table[10,2]<-count(filter(d, GO_slim_term == "receptor binding or activity"))
graph_table[11,2]<-count(filter(d, GO_slim_term == "other binding"))                  
graph_table[12,2]<-count(filter(d, GO_slim_term == "structural molecule activity"))
graph_table[13,2]<-count(filter(d, GO_slim_term == "transporter activity"))
graph_table[14,2]<-count(filter(d, GO_slim_term == "unknown molecular functions"))
graph_table[15,2]<-count(filter(d, GO_slim_term == "other molecular functions"))


write.csv(graph_table,file = paste('GOslim_665_T1syn.csv'))

unique_gene_list<-unique(T2_down$Gene)

d<-filter(t1, locus_name %in% unique_gene_list)
graph_table<-as.data.frame(matrix(nrow=15, ncol=2))
rownames(graph_table)<-c('hydrolase_activity', "kinase_activity" , "transferase_activity" , "other_enzyme_activity" , "transcription_factor_activity", "DNA_or_RNA_binding", "nucleic_acid_binding" , "nucleotide_binding",  "protein_binding" , "receptor_binding_or_activity", "other_binding", "structural_molecule_activity", "transporter_activity", "unknown_molecular_functions", "other_molecular_functions")

graph_table[1,2]<-count(filter(d, GO_slim_term == 'hydrolase activity'))
graph_table[2,2]<-count(filter(d, GO_slim_term == 'kinase activity'))
graph_table[3,2]<-count(filter(d, GO_slim_term == "transferase activity"))
graph_table[4,2]<-count(filter(d, GO_slim_term == "other enzyme activity"))
graph_table[5,2]<-count(filter(d, GO_slim_term == "transcription factor activity"))
graph_table[6,2]<-count(filter(d, GO_slim_term == "DNA or RNA binding"))
graph_table[7,2]<-count(filter(d, GO_slim_term == "nucleic acid binding"))
graph_table[8,2]<-count(filter(d, GO_slim_term == "nucleotide binding"))
graph_table[9,2]<-count(filter(d, GO_slim_term == "protein binding"))
graph_table[10,2]<-count(filter(d, GO_slim_term == "receptor binding or activity"))
graph_table[11,2]<-count(filter(d, GO_slim_term == "other binding"))                  
graph_table[12,2]<-count(filter(d, GO_slim_term == "structural molecule activity"))
graph_table[13,2]<-count(filter(d, GO_slim_term == "transporter activity"))
graph_table[14,2]<-count(filter(d, GO_slim_term == "unknown molecular functions"))
graph_table[15,2]<-count(filter(d, GO_slim_term == "other molecular functions"))


write.csv(graph_table,file = paste('GOslim_665_T2syn.csv'))



################### END ####################
















