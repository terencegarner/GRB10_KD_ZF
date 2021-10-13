library(ggplot2)
library(gplots)
library(BioQC)
library(RColorBrewer)

#Set working directory and load data

setwd("R:/Terry/Zebrafish") 
#Load data
ctrl<-read.table("ZF Control Full.txt",sep="\t",stringsAsFactors = F) 
morph<-read.table("ZF Morph Full.txt",sep="\t",stringsAsFactors = F)
#Load Variable Lists
var_ctrl<-read.table("ZF_Control_var_list.txt",sep="\t",stringsAsFactors = F,header=T)
var_morph<-read.table("ZF_Morph_var_list.txt",sep="\t",stringsAsFactors = F,header=T)

#Control----
#Format data object
cd<-ctrl[6:nrow(ctrl),6:ncol(ctrl)]
cd<-apply(cd,2,as.numeric)
rownames(cd)<-ctrl[-c(1:5),1]
colnames(cd)<-ctrl[5,-c(1:5)]

#Generate Hypernetwork Incidence Matrix
cd_cor<-cor(t(cd[match(var_ctrl$id,rownames(cd)),]),t(cd[-match(var_ctrl$id,rownames(cd)),])) #Correlate target genes against non-target genes
cd_bin_neg<-abs(cd_cor) #Convert values to absolutes 
cd_bin_neg[which(cd_cor<(-sd(cd_cor)))]<-1 #binarize negative correlations greater (more negative) than sd of all correlations 
cd_bin_neg[which(cd_bin_neg!=1)]<-0

cd_bin_pos<-abs(cd_cor) #Convert values to absolutes
cd_bin_pos[which(cd_cor>sd(cd_cor))]<-1 #binarize positive correlations greater (more positive) than sd of all correlations 
cd_bin_pos[which(cd_bin_pos!=1)]<-0

colnames(cd_bin_neg)<-paste("neg_",colnames(cd_bin_neg),sep="") #convert names of non-target genes to highlight positive vs negative correlates
colnames(cd_bin_pos)<-paste("pos_",colnames(cd_bin_pos),sep="")

cd_bin_comb<-cbind(cd_bin_pos,cd_bin_neg) #Combine positive and negative variables to generate the hypernetwork incidence matrix

hyp<-cd_bin_comb%*%t(cd_bin_comb) #Generate the hypernetwork adjacency matrix from incidence matrix

#Generate hypernetwork heatmap and identify clusters
hm<-heatmap.2(hyp,trace="none")
dend<-as.hclust(hm$rowDendrogram)
k<-3
ct<- cutree(dend, k)

#Identify central cluster (here cluster 1)
hyp_1<-hyp[which(ct==1),which(ct==1)]
#Calculate hypernetwork metrics
conn_1<-apply(hyp_1,1,mean)
entropy_1<-apply(hyp_1,1,entropy)
stats_c<-data.frame(conn_1,entropy_1/log2(86))#Normalize Entropy to the size of the list (here 86)

#Cluster Composition
c1<-names(ct[which(ct==1)])

#Identify Complete Subgraph
cd_gal_1<-cd_bin_comb[which(ct==1),]
cd_gal_1<-cd_gal_1[,which(colSums(cd_gal_1)==nrow(cd_gal_1))] #Identify which non-target genes were correlated with all target genes
cd_gal_1_names<-colnames(cd_gal_1)
cd_gal_1_names_unique<-cd_gal_1_names
cd_gal_1_names_unique<-substr(cd_gal_1_names_unique,5,nchar(cd_gal_1_names_unique))
cd_gal_1_names_unique<-unique(cd_gal_1_names_unique)


#Morpholino----
#Format data object
md<-morph[6:nrow(morph),6:ncol(morph)]
md<-apply(md,2,as.numeric)
rownames(md)<-morph[-c(1:5),1]
colnames(md)<-morph[5,-c(1:5)]

#Generate Hypernetwork Incidence Matrix
md_cor<-cor(t(md[match(var_morph$id,rownames(md)),]),t(md[-match(var_morph$id,rownames(md)),])) #Correlate target genes against non-target genes
md_bin_neg<-abs(md_cor) #Convert values to absolutes 
md_bin_neg[which(md_cor<(-sd(md_cor)))]<-1 #binarize negative correlations greater (more negative) than sd of all correlations 
md_bin_neg[which(md_bin_neg!=1)]<-0

md_bin_pos<-abs(md_cor) #Convert values to absolutes
md_bin_pos[which(md_cor>sd(md_cor))]<-1 #binarize positive correlations greater (more positive) than sd of all correlations 
md_bin_pos[which(md_bin_pos!=1)]<-0

colnames(md_bin_neg)<-paste("neg_",colnames(md_bin_neg),sep="") #convert names of non-target genes to highlight positive vs negative correlates
colnames(md_bin_pos)<-paste("pos_",colnames(md_bin_pos),sep="")

md_bin_comb<-cbind(md_bin_pos,md_bin_neg) #Combine positive and negative variables to generate the hypernetwork incidence matrix

hyp<-md_bin_comb%*%t(md_bin_comb) #Generate the hypernetwork adjacency matrix from incidence matrix

#Generate hypernetwork heatmap and identify clusters
hm<-heatmap.2(hyp,trace="none")
dend<-as.hclust(hm$rowDendrogram)
k<-3
ct<- cutree(dend, k)

#Identify central cluster (here cluster 1)
hyp_1<-hyp[which(ct==1),which(ct==1)]
#Calculate hypernetwork metrics
conn_1<-apply(hyp_1,1,mean)
entropy_1<-apply(hyp_1,1,entropy)
stats_m<-data.frame(conn_1,entropy_1/log2(65))#Normalize Entropy to the size of the list (here 65)

#Cluster Composition
c1<-names(ct[which(ct==1)])

#Identify Complete Subgraph
md_gal_1<-md_bin_comb[which(ct==1),]
md_gal_1<-md_gal_1[,which(colSums(md_gal_1)==nrow(md_gal_1))] #Identify which non-target genes were correlated with all target genes
md_gal_1_names<-colnames(md_gal_1)
md_gal_1_names_unique<-md_gal_1_names
md_gal_1_names_unique<-substr(md_gal_1_names_unique,5,nchar(md_gal_1_names_unique))
md_gal_1_names_unique<-unique(md_gal_1_names_unique)