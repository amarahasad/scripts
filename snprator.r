### extracting SNP information for 1001 genomes data genewise, tair10 annotation 

#load('anno_1135.rda')
#load('tair10.rda')
#load('snps_2029.rda')

#need to be in a folder where the X_1135 files are 
#could change script to use non_imputed data as well .
#gene<-'AT4G04740'



snporator<-function(gene) {

## getting start, stop and chr of the respective gene 
a1<-tair10[which(tair10[ ,5]==gene),2]
a2<-tair10[which(tair10[ ,5]==gene),3]
chr<-tair10[which(tair10[ ,5]==gene),1]

## subsetting anno
D<-anno[which(anno[,1]==chr&anno[,2]>a1&anno[,2]<a2),]
# D contains only SNPs of the gene of interest
h<-1

# because dim(anno)>dim(SNPs)
 for ( r in 1:nrow(D)) {
if(D[h,3]%in%SNPs$SNP==FALSE) { h=h+1} else {break}}
## loading the genotype data we need 
load(paste('X_1135_',SNPs[which(SNPs$SNP==D[h,3]),5],'.rda.gz',sep=''))
# extracting genotype data which are in the gene 
XX<-X[,colnames(X)%in%D$SNP]
rm(X)
j<-nrow(D)
 for ( r in 1:nrow(D)) {
if(D[j,3]%in%SNPs$SNP==FALSE) { j=j-1} else {break}}

### load respective genotype data
 if(SNPs[which(SNPs$SNP==D[j,3]),5]!=SNPs[which(SNPs$SNP==D[h,3]),5]) {
 
load(paste('X_1135_',SNPs[which(SNPs$SNP==D[nrow(D),3]),5],'.rda.gz',sep=''))
XX<-cbind(XX,X[,colnames(X)%in%D$SNP])
rm(X) }
## extracting SNPs for which we have genotype data 
D_<-subset(D,D$SNP%in%colnames(XX))

D_$count<-apply(XX,2,sum)
# get non synonymous SNPs
 Dns<-D_[grep('NON_S',D_[,6]),]
 
A<-data.frame(SNP=c('length','total','non','syn','start','stop','low','moderate','high'),Number=c(a2-a1,nrow(D_),length(grep('NON_S',D_[,6])),(length(grep('SYN',D_[,6]))-length(grep('NON_S',D_[,6]))),length(grep('START',D_[,6])),length(grep('STOP',D_[,6])),length(grep('LOW',D_[,6])),length(grep('MODERATE',D_[,6])),length(grep('HIGH',D_[,6]))))
 
V<-strsplit(Dns[,6],split="|",fixed=T)
for ( i in 1: nrow(Dns)) {
 Dns[i,8]<-paste(unique(V[[i]][grep('NON_S',V[[i]])+3]),collapse=',')}

Dns<-Dns[,-6]

X2<-XX[,colnames(XX)%in%Dns[,3]]

colnames(Dns)[7]<-'AA exchange'
## get accessions with alternative Allele 
Dns$accessions<-NA
 for ( z in 1:nrow(Dns)) {
 Dns[z,8]<-paste(rownames(X2)[which(X2[,z]==1)],collapse=',')}

out1<-paste(gene,'snp_summary.txt')
out2<-paste(gene,'non_syn_snps.txt')

write.table(A,file=out1,row.names=FALSE)
write.table(Dns,file=out2,row.names=FALSE)
cat('files',out1,out2,'have been generated!','\n')
}

