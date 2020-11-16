

library(plyr)
library(dplyr)
library(tidyr)
library(stringr)


### graph WITH BINS


results<-as.data.frame(matrix(nrow=300,ncol=5))
i=1:500*100
colnames(results)<-paste('bin',i, sep = '')
rownames(results)<-c('NON_significant', 'Significant', 'Highly_significant')

bin<-0.02
for(i in 1:50){
  tst<-which(df$rel_start>(i-1)*bin&df$rel_start<=(i*bin))
  tst1<-df[tst,]
  tst1$col<-tst1$p.value
 
  tst1$col[which(tst1$col>0.05)]<-1
  tst1$col[which(tst1$col<=0.05&tst1$col>=0.05/nrow(b))]<-2
  tst1$col[which(tst1$col>0.05/nrow(b))]<-3

  results[1,i]<-table(tst1$col== '1')["TRUE"]
  results[2,i]<-table(tst1$col== '2')["TRUE"]
  results[3,i]<-table(tst1$col== '3')["TRUE"]
}
write.csv(results, file=paste('exp_gra.csv', sep = ''))

####################################################################################################

##### SNP_WISE

####################################################################################################

############## Graph with 100 bps distance

load("~/Dokumente/Data_files/stop_gain_data.rda")
df<-read.csv('rel_pos_snps0.05.csv')
df<-b_
sb_sg<-stop_gain_data[,c(1,2,4,5)]
test<-merge(df, sb_sg, by = 'SNP')


cal_dist<-separate(test, SNP, c("A", "B"), sep = '\\- ')
cal_dist$pos_start<-abs((cal_dist$Start) - as.numeric(cal_dist$B))
cal_dist$pos_stop<-abs((cal_dist$Stop) - as.numeric(cal_dist$B))


### distance with start

results<-as.data.frame(matrix(nrow=3,ncol=1 + as.numeric(str_match(max(cal_dist$pos_start), "\\d\\d\\d"))))
i=1:1 + as.numeric(str_match(max(cal_dist$pos_start), "\\d\\d\\d"))
colnames(results)<-paste('bin',i, sep = '')
rownames(results)<-c('NON_significant', 'Significant', 'Highly_significant')

bin=100
for(i in 1:1 + as.numeric(str_match(max(cal_dist$pos_start), "\\d\\d\\d"))){
  tst<-which(cal_dist$pos_start>(i-1)*bin&cal_dist$pos_start<=(i*bin))
  tst1<-cal_dist[tst,]
  tst1$col<-tst1$`p-value`

  tst1$col[which(tst1$col>0.05)]<-1
  tst1$col[which(tst1$col<=0.05&tst1$col>0.05/nrow(b))]<-2
  tst1$col[which(tst1$col<0.05/nrow(b))]<-3
  
  results[1,i]<-table(tst1$col== '1')["TRUE"]
  results[2,i]<-table(tst1$col== '2')["TRUE"]
  results[3,i]<-table(tst1$col== '3')["TRUE"]
}
results[is.na(results)] <- 0

write.csv(results, file=paste('727_MT_50%strt.csv', sep = ''))

### distance with stop

results<-as.data.frame(matrix(nrow=3,ncol=1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d"))))
i=1:1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d"))
colnames(results)<-paste('bin',i, sep = '')
rownames(results)<-c('NON_significant', 'Significant', 'Highly_significant')


bin=100
for(i in 1:1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d"))){
  tst<-which(cal_dist$pos_stop>(i-1)*bin&cal_dist$pos_stop<=(i*bin))
  tst1<-cal_dist[tst,]
  tst1$col<-tst1$`p-value`
  
  tst1$col[which(tst1$col>0.05)]<-1
  tst1$col[which(tst1$col<=0.05&tst1$col>0.05/nrow(b))]<-2
  tst1$col[which(tst1$col<0.05/nrow(b))]<-3
  
  results[1,i]<-table(tst1$col== '1')["TRUE"]
  results[2,i]<-table(tst1$col== '2')["TRUE"]
  results[3,i]<-table(tst1$col== '3')["TRUE"]
}
results[is.na(results)] <- 0

write.csv(results, file=paste('727_MT_50%stp.csv', sep = ''))



####################################################################################################

##### GENE_WISE

####################################################################################################

############## Graph with 100 bps distance

load("~/Dokumente/Data_files/stop_gain_data.rda")
df<-read.csv('rel_pos_snps0.05.csv')



#sb_sg<-stop_gain_data[,c(1,2,4,5)]
#test<-merge(df, sb_sg, by = 'Gene')

### distance with start

results<-as.data.frame(matrix(nrow=3,ncol= 2 ))
#i=1 : (1 + as.numeric(str_match(max(cal_dist$pos_start), "\\d\\d\\d")))
#colnames(results)<-paste('bin',i, sep = '')
rownames(results)<-c('NON_significant', 'Significant', 'Highly_significant')

#bin=100
#for(i in 1: (1 + as.numeric(str_match(max(cal_dist$pos_start), "\\d\\d\\d")))){
#  tst<-which(cal_dist$pos_start>(i-1)*bin&cal_dist$pos_start<=(i*bin))
#  tst1<-test[tst,]
df$col<-df$`p-value`
  
  df$col[which(df$col>0.05)]<-1
  df$col[which(df$col<=0.05&df$col>0.05/nrow(b))]<-2
  df$col[which(df$col<0.05/nrow(b))]<-3
  
  NON_significant <-filter(df, df$col== '1')
  Significant <-filter(df, df$col== '2')
  Highly_significant <-filter(df, df$col== '3')



write.csv(results, file=paste('GW_727_50%strt.csv', sep = ''))

### distance with stop

results<-as.data.frame(matrix(nrow=3,ncol=(1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d")))))
i=1:(1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d")))
colnames(results)<-paste('bin',i, sep = '')
rownames(results)<-c('NON_significant', 'Significant', 'Highly_significant')


bin=100
for(i in 1:(1 + as.numeric(str_match(max(cal_dist$pos_stop), "\\d\\d\\d")))){
  tst<-which(cal_dist$pos_stop>(i-1)*bin&cal_dist$pos_stop<=(i*bin))
  tst1<-cal_dist[tst,]
  tst1$col<-tst1$`p-value`
  
  tst1$col[which(tst1$col>0.05)]<-1
  tst1$col[which(tst1$col<=0.05&tst1$col>0.05/nrow(b))]<-2
  tst1$col[which(tst1$col<0.05/nrow(b))]<-3
  
  results[1,i]<-table(tst1$col== '1')["TRUE"]
  results[2,i]<-table(tst1$col== '2')["TRUE"]
  results[3,i]<-table(tst1$col== '3')["TRUE"]
}
results[is.na(results)] <- 0

write.csv(results, file=paste('727_MT_50%stp.csv', sep = ''))



####################################################################################################














































##### Creating bar plots

### sample_data
dat <- read.table(text = "    ONE TWO THREE
+ 1   23  234 324
+ 2   34  534 12
+ 3   56  324 124
+ 4   34  234 124
+ 5   123 534 654",sep = "",header = TRUE)


datm <- melt(cbind(dat, ind = rownames(dat)), id.vars = c('ind'))

library(scales)
ggplot(datm,aes(x = variable, y = value,fill = ind)) + 
  +     geom_bar(position = "fill",stat = "identity") +
  +     # or:
  +     # geom_bar(position = position_fill(), stat = "identity") 
  +     scale_y_continuous(labels = percent_format())




barplot(as.matrix(Stress_vs_cont))



library(dplyr)
library(ggplot2)

### foo is your melt data
ana <- mutate(group_by(foo, x), percent = value / sum(value) * 100) %>%
  ungroup()

### Plot once
bob <- ggplot(data = ana, aes(x = x, y = percent, fill = variable)) +
  geom_bar(stat = "identity") +
  labs(y = "Percentage (%)")

### Get ggplot data
caroline <- ggplot_build(bob)$data[[1]]

### Create values for text positions
caroline$position = caroline$ymax + 1

### round up numbers and convert to character. Get unique values
foo <- unique(as.character(round(ana$percent, digits = 2)))

### Create a column for text
caroline$label <- paste(foo,"%", sep = "")

### Plot again
bob + annotate(x = caroline$x, y = caroline$position,
               label = caroline$label, geom = "text", size=3) 
































plot(1:652,df[,4])
df$col<-df$p.value
tst1$col<- -log10(tst1$p.value)
df$col<- round(df$col  )
plot(1:652,df[,4],col=df$col)
df$col[which(df$col<5)]<-1
df$col[which(df$col>10)]<-3
df$col[which(df$col>=5&df$col<=10)]<-2
colo<-c('black','blue','red')
plot(1:652,df[,4],col=colo[df$col])
plot(1:652,df[,4],col=colo[df$col],pch=16)

