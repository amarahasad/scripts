


library(reshape)
library(ggplot2)



GO_FILE1<-read.csv('/home/ama42bf/Dokumente/Data_files/727_data/727_full/SNPwise/GO_Sim_727_T1.csv', header = TRUE)
GO_File2<- read.csv('/home/ama42bf/Dokumente/Data_files/727_data/727_full/SNPwise/GO_Sim_727_T2.csv', header = TRUE)

t<-GO_FILE1[,c(2:7)]
t[,2]<-GO_FILE1$GENOME
t[,3]<-GO_FILE1$V2
t[,4]<-GO_File2$V2
t[,5]<-GO_FILE1$p.value
t[,6]<-GO_File2$p.value
colnames(t)<-c('X', 'GENOME', 'T1', 'T2', 'T1_pval', 'T2_pval')


t1<-t
t1$T1 = t1$T1 / sum(t1$T1)* 100
t1$GENOME = t1$GENOME / sum(t1$GENOME)* 100
t1$T2 = t1$T2 / sum(t1$T2)* 100

t2<-t1[,c(1:4)]
t3<-melt(data = t2, id.vars = c('X'))


t2$col1<-t1$T1_pval
t2$col2<-t1$T2_pval  
t2$col1[which(t1$T1_pval >= 0.05)]<-NA
t2$col1[which(t1$T1_pval <= 0.05&t1$T1_pval>0.0001)]<-'**'
t2$col1[which(t1$T1_pval <0.0001)]<-'***' 
t2$col2[which(t1$T2_pval >= 0.05)]<-NA
t2$col2[which(t1$T2_pval <= 0.05&t1$T2_pval>0.0001)]<-'**'
t2$col2[which(t1$T2_pval <0.0001)]<-'***' 




t3$col<-NA
t3$col[which(t3$variable == 'GENOME'&t3$value>0.000)]<- 'REF'
t3[16:30,4]<-t2$col1 
t3[31:45,4]<-t2$col2 


g <- ggplot(t3, aes(x= X, y=value)) + geom_bar(aes(fill = variable), stat = "identity", position = "dodge") + theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_text(aes(label = col, group = variable), vjust = -0.5, position = position_dodge(width = 1), size = 2 ) + theme(axis.text.x = element_text(face="bold", color="#993333", size=9, angle=45),axis.text.y = element_text(face="bold", color="#993333",  size=10, angle=45)) + theme( axis.line = element_line(colour = "darkblue", size = 1, linetype = "solid")) + scale_y_continuous(name="Value in Percentage" , breaks=seq(0,50,5))+ labs(title = "Significance level of 727 Accessions in GO_SLIM Terms", face="bold", color='BLACK', size=12) + theme(axis.title.x=element_blank())
cbbPalette <- c( "#0072B2","#E69F00", "#009E73")
g1<-g +  scale_fill_manual(values=cbbPalette)


pdf(file="GO_727_10.pdf", height=8, width=11)
g1
dev.off()


