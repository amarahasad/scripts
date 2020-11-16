####

library(plyr)
library(dplyr)
library(tidyr)
library(reshape)
library(ggplot2)

load('/home/ama42bf/Dokumente/Data_files/727_data/GO_Molecular_FUNCTIONS.Rda')

######## Up down Regulation

T2_727<-read.csv('727_T2.csv')
T2_727$dir<- round(T2_727$mean_KO) - round(T2_727$mean_WT)
T2_up<-subset(T2_727, dir < 0)
T2_down<-subset(T2_727, dir > 0)


genes_file<-T2_syn_down   #### change file name accordingly

unique_gene_list<-unique(genes_file$Gene)

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

write.csv(graph_table,file = paste('GOslim_665_T2syn1.csv'))

write.csv(T2_syn_down,file = paste('T2_syn_665_down.csv'))
write.csv(T2_syn_up,file = paste('T2_syn_665_up.csv'))
write.csv(T1_syn_down,file = paste('T1_syn_665_down.csv'))
write.csv(T1_syn_up,file = paste('T1_syn_665_up.csv'))



GO_FILE1<-read.csv('/home/ama42bf/Dokumente/Data_files/727_data/727_full/SNPwise/GO_727_T2_down.csv', header = TRUE)
GO_File2<- read.csv('/home/ama42bf/Dokumente/Data_files/727_data/727_full/SNPwise/GO_727_T2_up.csv', header = TRUE)

pvalue<-read.csv('GO_Sim_727_T2.csv')
genome_data<-read.csv('GO_Sim_727_T1.csv')


t<-as.data.frame(matrix(nrow = 15, ncol = 5))

t[,1]<-GO_FILE1[,1]

t[,2]<-genome_data$GENOME
t[,3]<-GO_FILE1$V2
t[,4]<-GO_File2$V2
t[,5]<-pvalue$p.value

colnames(t)<-c('X', 'GENOME', 'T2_down', 'T2_up', 'pval')

t1<-t
t1$T2_down = t1$T2_down / sum(t1$T2_down)* 100
t1$GENOME = t1$GENOME / sum(t1$GENOME)* 100
t1$T2_up = t1$T2_up / sum(t1$T2_up)* 100

t2<-t1[,c(1:4)]
t3<-melt(data = t2, id.vars = c('X'))


t2$col1<-t1$pval
t2$col1[which(t1$pval >= 0.05)]<-NA
t2$col1[which(t1$pval <= 0.05&t1$pval>6.738544^-06)]<-'**'
t2$col1[which(t1$pval <6.738544^-06)]<-'***' 




t3$col<-NA
t3$col[which(t3$variable == 'GENOME'&t3$value>0.000)]<- 'REF'
t3[16:30,4]<-t2$col1 
t3[31:45,4]<-t2$col1 


g <- ggplot(t3, aes(x= X, y=value)) + geom_bar(aes(fill = variable), stat = "identity", position = "dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_text(aes(label = col, group = variable), vjust = -0.5, position = position_dodge(width = 1), size = 2 ) + theme(axis.text.x = element_text(face="bold", color="#993333", size=9, angle=45),axis.text.y = element_text(face="bold", color="#993333",  size=10, angle=45)) + theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) + scale_y_continuous(name="Value in Percentage" , breaks=seq(0,50,5))+ labs(title = "Significance level of 727 Accessions in GO_SLIM Terms", face="bold", color='BLACK', size=12) + theme(axis.title.x=element_blank())
cbbPalette <- c( "#0072B2","#E69F00", "#009E73")
g1<-g +  scale_fill_manual(values=cbbPalette)

g1




pdf(file="GO_727_regulatory.pdf", height=8, width=11)
g1
dev.off()






data_sum<-as.data.frame(matrix(nrow = 2, ncol = ))
















