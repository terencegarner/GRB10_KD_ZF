library(ggplot2)
library(gplots)
library(BioQC)
library(RColorBrewer)

#Set working directory and load data
setwd("R:/Terry/Zebrafish")
ctrl<-read.table("ZF Control Full.txt",sep="\t",stringsAsFactors = F)
morph<-read.table("ZF Morph Full.txt",sep="\t",stringsAsFactors = F)
#Load variable lists 
pvals_ctrl<-read.table("ZF_control_pvals.txt",sep="\t",stringsAsFactors = F,header=T)
pvals_morph<-read.table("ZF_morph_pvals.txt",sep="\t",stringsAsFactors = F,header=T)

iter<-1000 #Set number of iterations

#Control----
#Format data object
cd<-ctrl[6:nrow(ctrl),6:ncol(ctrl)]
cd<-apply(cd,2,as.numeric)
rownames(cd)<-ctrl[-c(1:5),1]
colnames(cd)<-ctrl[5,-c(1:5)]

list_res<-list()

for(i in 1:iter){
  
  rand<-sample(pvals_ctrl$id,size = 86,replace = F) #Select random genes, the number of which is the same as the hypernetwork central cluster to be compared
  
  cd_cor<-cor(t(cd[match(rand,rownames(cd)),]),t(cd[-match(rand,rownames(cd)),])) #Correlate randomly selected genes against all others
  neg_cd_bin<-abs(cd_cor)  #Convert values to absolutes 
  neg_cd_bin[which(cd_cor<(-sd(cd_cor)))]<-1 #binarize negative correlations greater (more negative) than sd of all correlations 
  neg_cd_bin[which(neg_cd_bin!=1)]<-0
  
  pos_cd_bin<-abs(cd_cor)#Convert values to absolutes
  pos_cd_bin[which(cd_cor>sd(cd_cor))]<-1 #binarize positive correlations greater (more positive) than sd of all correlations 
  pos_cd_bin[which(pos_cd_bin!=1)]<-0
  
  cd_bin_comb<-cbind(pos_cd_bin,neg_cd_bin)#Combine positive and negative variables to generate the hypernetwork incidence matrix
  
  hyp<-cd_bin_comb%*%t(cd_bin_comb)#Generate the hypernetwork adjacency matrix from incidence matrix
  
  #Calculate hypernetwork metrics
  it<-rep(i,nrow(hyp))
  Sum<-apply(hyp,1,sum)
  Mean<-apply(hyp,1,mean)
  Entropy<-apply(hyp,2,entropy)
  res<-data.frame(it,Sum,Mean,Entropy/log2(86))
  
  colnames(res)<-c("Iteration","Sum","Mean","Entropy")
  
  list_res[[i]]<-res
}

res_c <- do.call(rbind,list_res)

#Morpholino----
#Format data object
md<-morph[6:nrow(morph),6:ncol(morph)]
md<-apply(md,2,as.numeric)
rownames(md)<-morph[-c(1:5),1]
colnames(md)<-morph[5,-c(1:5)]

list_res<-list()

for(i in 1:iter){
  rand<-sample(pvals_ctrl$id,size = 65,replace = F)#Select random genes, the number of which is the same as the hypernetwork central cluster to be compared
  
  md_cor<-cor(t(md[match(rand,rownames(md)),]),t(md[-match(rand,rownames(md)),])) #Correlate randomly selected genes against all others
  neg_md_bin<-abs(md_cor)#Convert values to absolutes
  neg_md_bin[which(md_cor<(-sd(md_cor)))]<-1#binarize negative correlations greater (more negative) than sd of all correlations 
  neg_md_bin[which(neg_md_bin!=1)]<-0
  
  pos_md_bin<-abs(md_cor)#Convert values to absolutes
  pos_md_bin[which(md_cor>sd(md_cor))]<-1#binarize positive correlations greater (more positive) than sd of all correlations 
  pos_md_bin[which(pos_md_bin!=1)]<-0
  
  md_bin_comb<-cbind(pos_md_bin,neg_md_bin)#Combine positive and negative variables to generate the hypernetwork incidence matrix
  
  hyp<-md_bin_comb%*%t(md_bin_comb)#Generate the hypernetwork adjacency matrix from incidence matrix
  
  #Calculate hypernetwork metrics
  it<-rep(i,nrow(hyp))
  Sum<-apply(hyp,1,sum)
  Mean<-apply(hyp,1,mean)
  Entropy<-apply(hyp,2,entropy)
  res<-data.frame(it,Sum,Mean,Entropy/log2(65))
  
  colnames(res)<-c("Iteration","Sum","Mean","Entropy")
  
  list_res[[i]]<-res
}

res_m <- do.call(rbind,list_res)
