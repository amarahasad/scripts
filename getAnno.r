### extract information from vcf files 
### download VCF for your gene of interest 
### rename it and it must be named like AT1G10000_snpeff.vcf  !
### Set home directory
### install.packages(c("tidyr", "devtools"))

library(GenomicFeatures)
library(VariantAnnotation)
library(DSR)
library(tidyr)
library(plyr)


getAnno<-function(vcffile='AT5G61590_snpeff.vcf', safe=TRUE){
out.name<-paste('Snp_output',unlist(strsplit(vcffile,split='_snp'))[[1]],sep='_')

vcf<-readVcf(vcffile,'TAIR10')     ### str(vcf) ------>  shows the content of the file 
EFF<-info(vcf)$EFF  ## gets predicted effects 

########## crate a matrix with nrow equals number of markers (vcf@fixed@nrows) ncol equals 6 ######

SNPEFF_MATRIX<-as.data.frame(matrix(nrow=vcf@fixed@nrows,ncol=7))
colnames(SNPEFF_MATRIX)=c('Varaint','Chr','Pos','Ref','Alt','AC','AN')  ## AN is number of Alleles sequenced, ## AC is Allele count ## Alt alternative ALLele 

SNPEFF_MATRIX[,1]<-rownames(info(vcf))
SNPEFF_MATRIX[,7]<-info(vcf)$AN
SNPEFF_MATRIX$Pos<-vcf@rowRanges@ranges@start
SNPEFF_MATRIX$Chr<-rep(vcf@rowRanges@seqnames@values)
SNPEFF_MATRIX$Ref<-as.character(fixed(vcf)[,1])
SNPEFF_MATRIX$Alt<-sapply(alt(vcf), paste, collapse = "|")   ### is prob slow see https://support.bioconductor.org/p/66874/
SNPEFF_MATRIX$AC<-sapply(info(vcf)$AC, paste, collapse = "|")

SNPEFF_EFFECT<-as.data.frame(matrix(nrow=vcf@fixed@nrows,ncol=1))
colnames(SNPEFF_EFFECT)=c('snpeff')
SNPEFF_EFFECT$snpeff<-sapply(info(vcf)$EFF, paste, collapse = ",")  # too much information

##TESTS FOR SNPEFF EFFECT INFORMATION
Split_data<-strsplit(SNPEFF_EFFECT$snpeff, split= ',', fixed = TRUE) ## splits the types of variants
SNPEFF_SPLITTED<-as.data.frame(unlist(Split_data)) ##  unlist the data into multiple rows and converts to data frame
SNPEFF_SPLITTED<-data.frame(cbind(id=rep(NA),snpeff=SNPEFF_SPLITTED)) ## adds a column in data frame with null values
Variant_repeat<-sapply(Split_data,length) ##  creates a number list to the total length of each variant type
SNP_name_list<-rep(SNPEFF_MATRIX[,1],Variant_repeat)  ###  takes the snp info from the column of other matrix and replicates it to the number list
SNPEFF_SPLITTED[,1]<-SNP_name_list  ### paste the list of duplicated names to the 1st column of data frame

SNPEFF_SPLITTED<-separate(SNPEFF_SPLITTED, unlist.Split_data. , into = c('Effect', 'Effect_Impact'),sep='\\(', convert = TRUE)
SNPEFF_SPLITTED<-separate(SNPEFF_SPLITTED, Effect_Impact , into = c('Effect_Impact','Functional_Class','Codon_Change','Amino_Acid_Change','Amino_Acid_length','Gene_Name','Transcript_BioType','Gene_Coding','Transcript_ID','Exon_Rank','Genotype_Number'),sep='\\|', convert = TRUE)

a<-ddply(SNPEFF_SPLITTED, "id", summarize, Effect = paste(Effect, collapse = ", "))
b<-ddply(SNPEFF_SPLITTED, "id", summarize, Effect_Impact = paste(Effect_Impact, collapse = ", "))
c<-ddply(SNPEFF_SPLITTED, "id", summarize, Functional_Class = paste(Functional_Class, collapse = ", "))
d<-ddply(SNPEFF_SPLITTED, "id", summarize, Codon_Change = paste(Codon_Change, collapse = ", "))
e<-ddply(SNPEFF_SPLITTED, "id", summarize, Amino_Acid_Changee = paste(Amino_Acid_Change, collapse = ", "))
f<-ddply(SNPEFF_SPLITTED, "id", summarize, Amino_Acid_length = paste(Amino_Acid_length, collapse = ", "))
g<-ddply(SNPEFF_SPLITTED, "id", summarize, Gene_Name = paste(Gene_Name, collapse = ", "))
h<-ddply(SNPEFF_SPLITTED, "id", summarize, Transcript_BioType = paste(Transcript_BioType, collapse = ", "))
i<-ddply(SNPEFF_SPLITTED, "id", summarize, Gene_Coding = paste(Gene_Coding, collapse = ", "))
j<-ddply(SNPEFF_SPLITTED, "id", summarize, Transcript_ID = paste(Transcript_ID, collapse = ", "))
k<-ddply(SNPEFF_SPLITTED, "id", summarize, Exon_Rank = paste(Exon_Rank, collapse = ", "))
SNPEFF_SPLITTED$Genotype_Number<-gsub( "\\)", "", as.character(SNPEFF_SPLITTED$Genotype_Number))
l<-ddply(SNPEFF_SPLITTED, "id", summarize, Genotype_Number = paste(Genotype_Number, collapse = ", "))
combine_df<-cbind(a[,2],b[,2],c[,2],d[,2],e[,2],f[,2],g[,2],h[,2],i[,2],j[,2],k[,2],l[,2])
colnames(combine_df)=c('Effect','Effect_Impact','Functional_Class','Codon_Change','Amino_Acid_Change','Amino_Acid_length','Gene_Name',
                       'Transcript_BioType','Gene_Coding','Transcript_ID','Exon_Rank','Genotype_Number')
SNPEFF_VARIANT_DATA<- cbind(SNPEFF_MATRIX,combine_df)

####Calculate Summary Statistics of particular gene_ID

Summary_Stats<-as.data.frame(matrix(nrow=1,ncol=8))
colnames(Summary_Stats)=c('No_of_SNPS','Coding_Region','Synonymous_SNPS','Non_Syn_SNPS','NO_of_STOP_GAINED','NO_of_STOP_LOST','NO_of_START_GAINED','NO_of_START_LOST')
Summary_Stats$No_of_SNPS<-length(Split_data)
Summary_Stats$Coding_Region<-length(grep('CODING',SNPEFF_VARIANT_DATA$Gene_Coding))
Summary_Stats$Synonymous_SNPS<-length(grep('SILENT',SNPEFF_VARIANT_DATA$Functional_Class))
Summary_Stats$Non_Syn_SNPS<-length(grep('MISSENSE',SNPEFF_VARIANT_DATA$Functional_Class))
Summary_Stats$NO_of_STOP_GAINED<-length(grep('stop_gained',SNPEFF_VARIANT_DATA$Effect))
Summary_Stats$NO_of_STOP_LOST<-length(grep('stop_lost',SNPEFF_VARIANT_DATA$Effect))
Summary_Stats$NO_of_START_GAINED<-length(grep('_start_codon_gain_',SNPEFF_VARIANT_DATA$Effect))
Summary_Stats$NO_of_START_LOST<-length(grep('_start_codon_lost_',SNPEFF_VARIANT_DATA$Effect))

if(safe==TRUE) { save(SNPEFF_VARIANT_DATA,file=paste(out.name,'.rda',sep=''))
  write.csv(SNPEFF_VARIANT_DATA,file=paste(out.name,'.csv',sep=''))
  write.csv(Summary_Stats,file=paste(out.name,'_Summary_stats.csv',sep=''))
  cat('files',paste(out.name,'.rda',sep=''),'and', paste(out.name,'.csv',sep=''),'are written to the following directory:','\n',getwd(),'\n')
}

return(SNPEFF_VARIANT_DATA)
}


##################################################________END________################################################





                           
                                           
 





