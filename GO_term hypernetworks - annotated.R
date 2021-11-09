library(corrplot)
library(gplots)
library(ggplot2)
library(biomaRt)
library(BioQC)

#Set working directory and load data
setwd("~/PEC/Terry/Zebrafish")
ctrl<-read.table("ZF Control Full.txt",sep="\t",stringsAsFactors = F)
morph<-read.table("ZF Morph Full.txt",sep="\t",stringsAsFactors = F)
#Load variable lists
pvals_ctrl<-read.table("ZF_control_pvals.txt",sep="\t",stringsAsFactors = F,header=T)
pvals_morph<-read.table("ZF_morph_pvals.txt",sep="\t",stringsAsFactors = F,header=T)

#Format data objects
cd<-ctrl[6:nrow(ctrl),6:ncol(ctrl)]
cd<-apply(cd,2,as.numeric)
rownames(cd)<-ctrl[-c(1:5),1]
colnames(cd)<-ctrl[5,-c(1:5)]
cd<-cd[na.omit(match(pvals_ctrl$id,rownames(cd))),]

md<-morph[6:nrow(morph),6:ncol(morph)]
md<-apply(md,2,as.numeric)
rownames(md)<-morph[-c(1:5),1]
colnames(md)<-morph[5,-c(1:5)]
md<-md[na.omit(match(pvals_morph$id,rownames(md))),]

#Load GO term lists and convert human genes to ZF homologues
setwd("~/PEC/Terry/Zebrafish/Paper/Aug 2021/Ontology Lists/hsapiens/")
lists<-list.files(pattern=".txt")
lists<-lapply(lists,read.delim,header=F)
genelists<-lapply(lists,"[[",3)
genelists<-genelists[which(unlist(lapply(genelists,length))>=15)]
convert<-function(x){temp<-getBM(attributes=c('external_gene_name','ensembl_gene_id','drerio_homolog_ensembl_gene','drerio_homolog_associated_gene_name'), filters = 'external_gene_name', values = x, mart = ensembl)
temp<-unique(temp$drerio_homolog_associated_gene_name);temp<-temp[temp!=""]
}
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
genelists<-lapply(genelists,convert)


#Set up hypernetwork iteration parameters
n<-1000 #number of iterations
hyp_size<-10 #number of genes to use per iteration

#Generate hypernetworks using Control data----
dat<-cd

#Format results data frame before analysis
res<-as.data.frame(matrix(ncol=3,nrow=length(genelists)*n))
colnames(res)<-c("Process","Entropy","Network_size")
res$Process<-rep(c(),each=n)
res$List_length<-rep(unlist(lapply(genelists,length)),each=n)

#Populate GO term field with imported GO terms
res$Process<-rep(c("Actin Filament Based Movement",
                   "Amino Acid Transport",
                   "Cellular Modified Amino Acid Biosynthetic Process",
                   "Collagen Fibril Organization",
                   "Dicarboxylic Acid Transport",
                   "Extracellular Structure Organization",
                   "Fat Soluble Vitamin Metabolic Process",
                   "Glutathione Netabolic Process",
                   "Metanephros Development",
                   "Negative Regulation Of Extrinsic Apoptotic Signaling Pathway",
                   "Platelet Activation",
                   "Positive Regulation Of Lipid Localization",
                   "Positive Regulation Of Wound Healing",
                   "Regulation Of mRNA Processing",
                   "Regulation Of RNA Splicing",
                   "Regulation Of Sterol Transport",
                   "Ribosome Assembly"),each=n)

#Iterate hypernetworks separately for each GO term----

#GC Response
ent<-numeric()
lens<-numeric()

target<-genelists[[1]]
for(i in 1:n){

  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])  #Generate a unique list of genes associated with the given GO term, which is also found in the dataset.
  genes<-sample(genes,size = hyp_size,replace=F) #Select a random set of these genes to generate the hypernetwork from
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))#Generate the correlation matrix between randomly selected target genes and all other genes 
  bin<-abs(cordat) #Convert values to absolutes
  
  bin[which(bin>sd(cordat)*1)]<-1 # Binarize the correlation values around the standard deviation of the correlation values
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin) #Generate the hypernetwork adjacency matrix from incidence matrix
  hist(hyp)

  hm<-heatmap.2(hyp,trace="none") #Generate hypernetwork heatmap
  ent[i]<-entropy(hyp)/log(length(hyp)) #Calculate hypernetwork metrics
  lens[i]<-length(hyp)
}
res$Entropy[1:n]<-ent #Populate results dataframe
res$Network_size[1:n]<-lens

#Growth
ent<-numeric()
lens<-numeric()
target<-genelists[[2]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[(n+1):(2*n)]<-ent
res$Network_size[(n+1):(2*n)]<-lens

#Gamete Generation
ent<-numeric()
lens<-numeric()
target<-genelists[[3]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((2*n)+1):(3*n)]<-ent
res$Network_size[((2*n)+1):(3*n)]<-lens

#Aging
ent<-numeric()
lens<-numeric()
target<-genelists[[4]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((3*n)+1):(4*n)]<-ent
res$Network_size[((3*n)+1):(4*n)]<-lens

#Regulation of Metabolic Processes
ent<-numeric()
lens<-numeric()
target<-genelists[[5]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((4*n)+1):(5*n)]<-ent
res$Network_size[((4*n)+1):(5*n)]<-lens


#Regulation of Hormone Metabolic Processes
ent<-numeric()
lens<-numeric()
target<-genelists[[6]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((5*n)+1):(6*n)]<-ent
res$Network_size[((5*n)+1):(6*n)]<-lens

# #Locomotion
ent<-numeric()
lens<-numeric()
target<-genelists[[7]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((6*n)+1):(7*n)]<-ent
res$Network_size[((6*n)+1):(7*n)]<-lens

#Biomineralization
ent<-numeric()
lens<-numeric()
target<-genelists[[8]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((7*n)+1):(8*n)]<-ent
res$Network_size[((7*n)+1):(8*n)]<-lens

#Regulation of DNA-binding transcription factor activity
ent<-numeric()
lens<-numeric()
target<-genelists[[9]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((8*n)+1):(9*n)]<-ent
res$Network_size[((8*n)+1):(9*n)]<-lens


#Protein folding
ent<-numeric()
lens<-numeric()
target<-genelists[[10]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((9*n)+1):(10*n)]<-ent
res$Network_size[((9*n)+1):(10*n)]<-lens


#Microtubule-based process
ent<-numeric()
lens<-numeric()
target<-genelists[[11]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((10*n)+1):(11*n)]<-ent
res$Network_size[((10*n)+1):(11*n)]<-lens

#Cell adhesion
ent<-numeric()
lens<-numeric()
target<-genelists[[12]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((11*n)+1):(12*n)]<-ent
res$Network_size[((11*n)+1):(12*n)]<-lens

#Response to insulin
ent<-numeric()
lens<-numeric()
target<-genelists[[13]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((12*n)+1):(13*n)]<-ent
res$Network_size[((12*n)+1):(13*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[14]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((13*n)+1):(14*n)]<-ent
res$Network_size[((13*n)+1):(14*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[15]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((14*n)+1):(15*n)]<-ent
res$Network_size[((14*n)+1):(15*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[16]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((15*n)+1):(16*n)]<-ent
res$Network_size[((15*n)+1):(16*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[17]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((16*n)+1):(17*n)]<-ent
res$Network_size[((16*n)+1):(17*n)]<-lens



#Generate violin plot of entropy across each process
res$Process<-as.factor(res$Process)

res_cd<-res
ggplot(res_cd,aes(x=Process,y=Entropy))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Repeat process using morpholino data----
dat<-md

res<-as.data.frame(matrix(ncol=3,nrow=length(genelists)*n))
colnames(res)<-c("Process","Entropy","Network_size")
res$Process<-rep(c(),each=n)
res$List_length<-rep(unlist(lapply(genelists,length)),each=n)

res$Process<-rep(c("Actin Filament Based Movement",
                   "Amino Acid Transport",
                   "Cellular Modified Amino Acid Biosynthetic Process",
                   "Collagen Fibril Organization",
                   "Dicarboxylic Acid Transport",
                   "Extracellular Structure Organization",
                   "Fat Soluble Vitamin Metabolic Process",
                   "Glutathione Netabolic Process",
                   "Metanephros Development",
                   "Negative Regulation Of Extrinsic Apoptotic Signaling Pathway",
                   "Platelet Activation",
                   "Positive Regulation Of Lipid Localization",
                   "Positive Regulation Of Wound Healing",
                   "Regulation Of mRNA Processing",
                   "Regulation Of RNA Splicing",
                   "Regulation Of Sterol Transport",
                   "Ribosome Assembly"),each=n)


#GC Response
ent<-numeric()
lens<-numeric()

target<-genelists[[1]]
for(i in 1:n){
  
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[1:n]<-ent
res$Network_size[1:n]<-lens

#Growth
ent<-numeric()
lens<-numeric()
target<-genelists[[2]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[(n+1):(2*n)]<-ent
res$Network_size[(n+1):(2*n)]<-lens

#Gamete Generation
ent<-numeric()
lens<-numeric()
target<-genelists[[3]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((2*n)+1):(3*n)]<-ent
res$Network_size[((2*n)+1):(3*n)]<-lens

#Aging
ent<-numeric()
lens<-numeric()
target<-genelists[[4]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((3*n)+1):(4*n)]<-ent
res$Network_size[((3*n)+1):(4*n)]<-lens

#Regulation of Metabolic Processes
ent<-numeric()
lens<-numeric()
target<-genelists[[5]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((4*n)+1):(5*n)]<-ent
res$Network_size[((4*n)+1):(5*n)]<-lens


#Regulation of Hormone Metabolic Processes
ent<-numeric()
lens<-numeric()
target<-genelists[[6]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((5*n)+1):(6*n)]<-ent
res$Network_size[((5*n)+1):(6*n)]<-lens

# #Locomotion
ent<-numeric()
lens<-numeric()
target<-genelists[[7]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((6*n)+1):(7*n)]<-ent
res$Network_size[((6*n)+1):(7*n)]<-lens

#Biomineralization
ent<-numeric()
lens<-numeric()
target<-genelists[[8]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((7*n)+1):(8*n)]<-ent
res$Network_size[((7*n)+1):(8*n)]<-lens

#Regulation of DNA-binding transcription factor activity
ent<-numeric()
lens<-numeric()
target<-genelists[[9]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((8*n)+1):(9*n)]<-ent
res$Network_size[((8*n)+1):(9*n)]<-lens


#Protein folding
ent<-numeric()
lens<-numeric()
target<-genelists[[10]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((9*n)+1):(10*n)]<-ent
res$Network_size[((9*n)+1):(10*n)]<-lens


#Microtubule-based process
ent<-numeric()
lens<-numeric()
target<-genelists[[11]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((10*n)+1):(11*n)]<-ent
res$Network_size[((10*n)+1):(11*n)]<-lens

#Cell adhesion
ent<-numeric()
lens<-numeric()
target<-genelists[[12]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((11*n)+1):(12*n)]<-ent
res$Network_size[((11*n)+1):(12*n)]<-lens

#Response to insulin
ent<-numeric()
lens<-numeric()
target<-genelists[[13]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((12*n)+1):(13*n)]<-ent
res$Network_size[((12*n)+1):(13*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[14]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_morph$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_morph$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_morph$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((13*n)+1):(14*n)]<-ent
res$Network_size[((13*n)+1):(14*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[15]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((14*n)+1):(15*n)]<-ent
res$Network_size[((14*n)+1):(15*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[16]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((15*n)+1):(16*n)]<-ent
res$Network_size[((15*n)+1):(16*n)]<-lens

#Insulin secretion
ent<-numeric()
lens<-numeric()
target<-genelists[[17]]
for(i in 1:n){   
  genes<-unique(target[target %in% pvals_ctrl$Gene.Symbol])
  genes<-sample(genes,size = hyp_size,replace=F)
  cordat<-cor(t(dat[na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]),t(dat[-na.omit(match(genes,pvals_ctrl$Gene.Symbol)),]))
  bin<-abs(cordat)
  
  bin[which(bin>sd(cordat)*1)]<-1
  bin[which(bin!=1)]<-0
  
  hyp<-bin%*%t(bin)
  hist(hyp)
  hm<-heatmap.2(hyp,trace="none")
  ent[i]<-entropy(hyp)/log(length(hyp))
  lens[i]<-length(hyp)
}
res$Entropy[((16*n)+1):(17*n)]<-ent
res$Network_size[((16*n)+1):(17*n)]<-lens

#Plot
res$Process<-as.factor(res$Process)

res_md<-res
ggplot(res_md,aes(x=Process,y=Entropy))+
  geom_violin()+
  geom_boxplot(width=0.2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#Generate violin plot comparing entropy between datasets

res_comb<-rbind(res_cd,res_md)
res_comb$Group<-as.factor(c(rep("Ctrl",nrow(res_cd)),rep("Morph",nrow(res_md))))

ggplot(res_comb,aes(x=Process,y=Entropy,col=Group))+
  geom_violin(position = position_dodge(width = 0.9),width=1)+
  geom_boxplot(position = position_dodge(width = 0.9),width=0.1)+
  theme_grey(base_size=22)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))