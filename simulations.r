

load('all_chr_net_data.Rda')
load('/home/ama42bf/Dokumente/Data_files/727_data/GO_Molecular_FUNCTIONS.Rda')
library(plyr)
library(dplyr)
library(tidyr)
library(ggrepel)


GO_SLIM_term<-function(GO_SLIM<-read.delim(file = 'ATH_GO_GOSLIM.txt', header = FALSE)){
  
  subset_GO_SLIM<-GO_SLIM
  subset_GO_SLIM<-subset_GO_SLIM[,1:6]
  subset_GO_SLIM[,2:5]<-NULL
  subset_GO_SLIM<-subset_GO_SLIM[1:277537,1:2]
  
  mol_fun<-c('GO:0016787','GO:0016301','GO:0016740','GO:0003824','GO:0003700','GO:0003677','GO:0003723','GO:0003676','GO:0000166','GO:0005515','GO:0005102','GO:0004872','GO:0005488','GO:0005198','GO:0005215','GO:0005554','GO:0003674')
  
  t1<-filter(subset_GO_SLIM, GO_ID %in% mol_fun)
  
  data_subset<-subset(all_chr_net_data,(all_chr_net_data$p_value_over < 10^-9))
  
  suba<-data.frame ( table ( data_subset$gene1 ) )
  subb<-data.frame ( table ( data_subset$gene2 ) )
  s_a<-suba[!(suba$Freq == 0),]
  s_b<-subb[!(subb$Freq == 0),]
  a<-rbind(s_a, s_b)
  unique_gene_list<-unique(a$Var1)
  slim_data<-filter(t1, locus_name %in% unique_gene_list)
  
  
  sample_data <- replicate(1000, sample(nrow(all_chr_net_data), "10", replace=FALSE,  prob = NULL ))
  
  simul_results<-as.data.frame(matrix(nrow=15,ncol=1000))
  
  for(i in 1:1000 ){
    d<-filter(t1, locus_name %in% under50_sam[,i])
    simul_results[1,i]<-count(filter(d, GO_ID == 'GO:0016787'))
    simul_results[2,i]<-count(filter(d, GO_ID == 'GO:0016301'))
    simul_results[3,i]<-count(filter(d, GO_ID == 'GO:0016740'))
    simul_results[4,i]<-count(filter(d, GO_ID == 'GO:0003824'))
    simul_results[5,i]<-count(filter(d, GO_ID == 'GO:0003700'))
    simul_results[6,i]<-count(filter(d, GO_ID == 'GO:0003677')) + count(filter(d, GO_ID == 'GO:0003723'))
    simul_results[7,i]<-count(filter(d, GO_ID == 'GO:0003676'))
    simul_results[8,i]<-count(filter(d, GO_ID == 'GO:0000166'))
    simul_results[9,i]<-count(filter(d, GO_ID == 'GO:0005515'))
    simul_results[10,i]<-count(filter(d, GO_ID == 'GO:0005102')) + count(filter(d, GO_ID == 'GO:0004872'))
    simul_results[11,i]<-count(filter(d, GO_ID == 'GO:0005488'))                  
    simul_results[12,i]<-count(filter(d, GO_ID == 'GO:0005198'))
    simul_results[13,i]<-count(filter(d, GO_ID == 'GO:0005215'))
    simul_results[14,i]<-count(filter(d, GO_ID == 'GO:0005554'))
    simul_results[15,i]<-count(filter(d, GO_ID == 'GO:0003674'))
  }
  
  simul_under_50_<-as.data.frame(t(simul_under_50))
  colnames(simul_under_50_)<-colnames(under_gra)
  
  
  
  for(i in 1:15){
    simul_results[i,5]<-mean(simul_under_50_[,i])
  }
  
  
  ### calculate pvalue under
  for ( i in 1:15 ) {
    simul_results_[i,11]<-length(which(simul_over_9[i,1:1000]<=simul_results_[i,1]))/1000
    simul_results_[i,12]<-length(which(simul_over_9[i,1:1000]>=simul_results_[i,1]))/1000
  }
  
  ### calculate pvalue over
  for ( i in 1:15 ) {
    under_gra_[z,8]<-length(which(simul[z,2:1001]>=under_gra_[z,2]))/1000}
  
  
  
  
  over_50_slim_data<-filter(t1, locus_name %in% over50_list)
  
  simul_results[1,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0016787'))
  simul_results[2,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0016301'))
  simul_results[3,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0016740'))
  simul_results[4,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0003824'))
  simul_results[5,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0003700'))
  simul_results[6,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0003677')) + count(filter(over_50_slim_data, GO_ID == 'GO:0003723'))
  simul_results[7,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0003676'))
  simul_results[8,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0000166'))
  simul_results[9,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005515'))
  simul_results[10,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005102')) + count(filter(over_50_slim_data, GO_ID == 'GO:0004872'))
  simul_results[11,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005488'))                  
  simul_results[12,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005198'))
  simul_results[13,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005215'))
  simul_results[14,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0005554'))
  simul_results[15,10]<-count(filter(over_50_slim_data, GO_ID == 'GO:0003674'))
  
  
  
  dat2 = data.frame(count=simul_results[,5], category=rownames(simul_results))
  dat2$percentage = dat2$count / sum(dat2$count)* 100
  dat2$fraction = dat2$count / sum(dat2$count)
  dat2 = dat2[order(dat2$fraction), ]
  dat2$ymax = cumsum(dat2$fraction)
  dat2$ymin = c(0, head(dat2$ymax, n=-1))
  
  donut = ggplot(dat2, aes(fill = category, ymax = ymax, ymin = ymin, xmax = 100, xmin = 80)) +
    geom_rect(colour = "black") +
    coord_polar(theta = "y") + 
    xlim(c(0, 100)) +
    geom_label_repel(aes(label = paste(round(percentage,2),"%"), x = 100, y = (ymin + ymax)/2),inherit.aes = F, show.legend = F, size = 5)+
    theme(legend.title = element_text(colour = "black", size = 16, face = "bold"), 
          legend.text = element_text(colour = "black", size = 15), 
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    annotate("text", x = 0, y = 0, size = 6, label = "sim_under-50")
  donut
  
  ## to save plot as pdf
  pdf(file='sim_under-50.pdf')
  donut
  dev.off()
  
  
 
p2 = ggplot(dat2, aes(fill=category, ymax=ymax, ymin=ymin, xmax=4, xmin=3)) +
    geom_rect(colour="grey30") +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    theme_bw() +
    theme(panel.grid=element_blank()) +
    theme(axis.text=element_blank()) +
    theme(axis.ticks=element_blank()) +
    labs(title="test1")
  
  
}
  
  





