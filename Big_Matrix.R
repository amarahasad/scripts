library(dplyr)

acc <- read.csv("accession_centric.csv", sep = ",")
d_  <- read.csv2("gene_data.csv", sep = ",")
d_ <- d_[,2:3]
colnames(d_)<-c('Gene', 'Accession')
acc <- acc[,2:4]
c<-d_$Gene
setwd("/home/ama42bf/Downloads/Matrix3/")
### So now i uploaded the stuff we need.


Tim <- Sys.time()
### the first if condition is for the begninng to create the matrix A with the first ten lines
### it should run when R does not have the Matrix A in its environment, if you look futher down
### you can see that i delete the Matrix A after saving it so that R doesnt have to allocated it everytime
### so after the first run the second loop will come to action and its resulting Matrix "B"
### will be rbinded to the Matrix A 
for(s in seq(1,9991,by=10)){
  Tim <- Sys.time()
for (s in 1:100) {
   A <- matrix(NA, nrow = 100,ncol = 9999)
   colnames(A)<-c
    for(i in s:(s+9)){
      for (k in 1:9999){
        A[i-s+1,k] <- length(which(strsplit(as.character(d_$Accession[i]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
      }
    }  }
   Sys.time()-Tim

   Tim <- Sys.time()
 
     A <- matrix(NA, nrow = 100,ncol = 9999)
     colnames(A)<-c
     for (s in 1:100) {
       for (k in 1:9999){
         A[s,k] <- length(which(strsplit(as.character(d_$Accession[s]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
         
       }
       cat(s,'\n')
     }  
   Sys.time()-Tim
   
   
   
   
   
      ### second is performed whenever s is not 1 so that i can just rbind matrix b to a over and over again   
}else if (s<9991 && s>1 && (exists("A")==TRUE)){
  B <- matrix(NA, nrow = 10,ncol = 9999)
  colnames(B)<-c
  for(i in s:(s+9)){
    for (k in 1:9999){
      B[i-s+1,k] <- length(which(strsplit(as.character(d_$Accession[i]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
    }
  }
A <- rbind(A,B)  
  ### the last window is only 9 rows long so i just an extra else if for that one
}else if  (s==9991) {
  B <- matrix(NA, nrow = 9,ncol = 9999)
  colnames(B)<-c
  for(i in s:(s+8)) {
    for (k in 1:9999) {
      B[i-s+1,k] <- length(which(strsplit(as.character(d_$Accession[i]),",")[[1]]%in%strsplit(as.character(d_$Accession[k]),",")[[1]]))
    }
}
A <- rbind(A,B)
}
  
#### Here if the nrow of the Matrix Exceeds 500 i delete it and save it
### therefore in the end we will end up with 200 matricies but we can bind them to the desierd 
  ## 9999*9999 matrix later to save computing power.
if((dim(A)[1])>499 && s<9991){
  filename <- paste("Matrix_",s,".txt",sep="")
  write.table(A,filename)
  rm(A)
}
if (s == 9991){
  filename <- paste("Matrix_",s,".txt",sep="")
  write.table(A,filename)
  rm(A)
  rm(B)
}
}
End <- Sys.time()- Tim


