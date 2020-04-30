########Loading packages##################################
library(regioneR)
library(rtracklayer)
library(stringr)
library("seqinr")
library(BiocGenerics)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(tidyverse)
library("reader")
library(ggplot2)
library(pheatmap)
library(ChIPseeker)
library(clusterProfiler)
library(data.table)
library(fitdistrplus)
library(logspline)

##########Loading peak-called ATAC/CHip-seq data##########
#setwd("~/OneDrive/uclPHD/meta-data/experiment/geo/geodatadescription")
#setwd("C:/Users/Mengjia Li/OneDrive/uclPHD/meta-data/experiment/geo/geodatadescription")

#collect<-dir()
#store<-list()
#ID1<-list()
#description<-data.frame(matrix(NA, nrow=length(collect), ncol=7))
#colnames(description)<-c("ID","ERX","title","antibody","name","tissue","genotype")


#replacena<-function(x){
  #if (length(x)==0){
    #x<-as.character("NA")
  #}else
  #{x<-x}
#}


#for (i in 1:length(collect)){
  
#try({  

  #get the file address
  #address<-paste(getwd(), "/", collect[i], sep='')
  
  #get the file (NCBI) ID
  #ID<-str_extract(collect[i], "^[\\s\\S]*?[.]")
  #ID1[[i]]<-str_extract(collect[i], "^[\\s\\S]*?[.]")
  
  #read the description
  #store[[i]]<-read.table(address,header = FALSE)
  
  
  #get the ERX description
  #GSM<-as.vector(sapply(store[[i]], function(x) str_extract(x, "GSM[0-9]*:[\\s\\S]*")))
  #GSM<-as.character(GSM[complete.cases(GSM)])
  #GSM1<-replacena(GSM)
  
  #try to get the object1
  #title<-as.vector(sapply(store[[i]], function(x) str_extract(x, "Title=[\\s\\S]*")))
  #title<-as.character(title[complete.cases(title)])
  #title1<-replacena(title)
  
  #try to get the object2
  #antibody<-as.vector(sapply(store[[i]], function(x) str_extract(x, "antibody=[\\s\\S]*")))
  #antibody<-as.character(antibody[complete.cases(antibody)])
  #antibody1<-replacena(antibody)
  #antibody1<-paste(antibody1, sep="",collapse="")
  
  #try to get the onject3
  #name<-as.vector(sapply(store[[i]], function(x) str_extract(x, "name=[\\s\\S]*")))
  #name<-as.character(name[complete.cases(name)])
  #name1<-replacena(name)
  #name1<-paste(name1, sep="",collapse="")
  
  #try to get the tissue
  #tissue<-as.vector(sapply(store[[i]], function(x) str_extract(x, "tissue=[\\s\\S]*")))
  #tissue<-as.character(tissue[complete.cases(tissue)])
  #tissue1<-replacena(tissue)
  
  #try to get the genotype
  #genotype<-as.vector(sapply(store[[i]], function(x) str_extract(x, "genotype=[\\s\\S]*")))
  #genotype<-as.character(genotype[complete.cases(genotype)])
  #genotype1<-replacena(genotype)
  
  #},silent=FALSE, outFile=stdout()
#)
  #combine then togetger
 #description[i,]<-c(ID,GSM1,title1,antibody1,name1,tissue1,genotype1)
#cat(i, "\n")
  
#}
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo")
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo")
raw<-read.csv("dmtable.csv")
description<-raw[,c(1:5)]
colnames(description)[1]<-"ID"
####################load the bed file and convert to gRange#############################
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo/geodatacollection")
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo/geodatacollection")

collect2<-dir()
fileid<-list()
filetable<-list()
filetableGrange<-list()
interval<-list()

for (i in 1:length(collect2)){
  
  try({ 
    
    
    #get the file address
    address<-paste(getwd(), "/", collect2[i], sep='')
    
    #get the file (NCBI) ID
    fileid[[i]]<-str_extract(collect2[i], "^[\\s\\S]*?[.]")
    
    #read the gff file as data.frame
    filetableGrange[[i]]<-toGRanges(address)
    
    #get the numbers of interval
    interval[[i]]<-length(filetableGrange[[i]])
    
    
  }, silent=FALSE, outFile=stdout())
  
  cat(i, "\n")
  
}



##################basic#################################################
setwd("~/OneDrive/uclPHD/meta-data/experiment/foxo/results of coordination")
setwd("C:/Users/Mengjia Li/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/foxo/results of coordination")

pt2<-list()
dmgenome<-characterToBSGenome("dm3")

foxo<-read.csv("foxofb1.csv")

foxonoNA<-foxo[complete.cases(foxo),]

foxopadj<-subset(foxo, foxo$padj<0.1)

foxopadjnoNA<-foxopadj[complete.cases(foxopadj),]

foxopadjnoNA$chromosome<-droplevels(foxopadjnoNA$chromosome)

foxoopen<-subset(foxopadjnoNA, foxopadjnoNA$log2FoldChange>0)

foxogRange<-makeGRangesFromDataFrame(foxopadjnoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxofullgRange<-makeGRangesFromDataFrame(foxonoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

#foxoopengRange<-makeGRangesFromDataFrame(foxoopen,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

###############################################################
#permutation test
##############################################################
i<-1 

for (i in 1:length(collect2)){
  
  
  try(pt2[[i]]<-overlapPermTest(A=foxogRange, B=filetableGrange[[i]], ntime=2000,genome=dmgenome, alternative = "auto"))
  
  cat(i, "\n")
  
}

psig<-list()
observed<-list()
Zscore<-list()


for (i in 1:length(pt2)){
  
  psig[[i]]<-pt2[[i]]$numOverlaps$pval
  observed[[i]]<-pt2[[i]]$numOverlaps$observed
  Zscore[[i]]<-pt2[[i]]$numOverlaps$zscore
  
}

fileid1<-as.vector(unlist(fileid))
fileid2<-str_replace(fileid1,".$","")

psig<-as.vector(unlist(psig))
observed<-as.vector(unlist(observed))
Zscore<-as.vector(unlist(Zscore))


padj<-p.adjust(psig, method="BH")


result0<-as.data.frame(cbind(fileid1,psig, padj, observed, Zscore))
colnames(result0)<-c("ID","psig","padj","observed","Zscore")
result<-left_join(description, result0, by="ID")
result<-result[complete.cases(result),]

##############output###########################################
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo")
setwd("C:/Users/Mengjia Li/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo")

write.csv(result, file=paste(getwd(), "/", "1000_onlyup_geofirsttry_permtest.csv", sep=''))

############modify previous experiment##########################
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo result")
setwd("C:/Users/Mengjia Li/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo result")

result<-read.csv(paste(getwd(), "/", "2000_all_geo_permtest.csv", sep=''))#4046 obejcts
result$ID<-str_replace(result2$ID,".$","") #only need when read those have ./ at the end of ID
result1<-left_join(description, result, by="ID")
result2<-result1[!is.na(result1$psig),] #only take those instruciton and experiment exist #need to know why the left experiment disappear
result3<-result2[,c(1,3,4,5,9,13,14,15,16)]
write.csv(result3,paste(getwd(),"/","2000_all_geo_permtest_updateinfo.csv",sep=''))


#################################################################
#FURTHER ANALYSIS OF THE GEO DATASET
#################################################################


#DELETED THE CONTROL SIG
setwd("C:/Users/Mengjia Li/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/geo")
result<-read.csv("5000_onlyup_geo_permtest.csv")#4048 object
control<-read.csv("1000_control_geo_permtest.csv")#4047 object

control$ID<-str_replace(control$ID,".$","")#only need when read those have ./ at the end of ID
setdiff(result$ID,control$ID) #SRX1531766 padj=1, so ignore

controlsig<-subset(control,control$padj<0.04) #223
resultsig<-subset(result,result$padj<0.1) #219

bothsig<-intersect(controlsig$ID,resultsig$ID)
subset(resultsig,resultsig$ID %in% bothsig) #trx*1,Trl*6,Creb1*3,Pol2*3,Input/unclassified*2, Lilli*1, Rpb1*1
remain<-setdiff(result$ID,bothsig)

result1<-subset(result,result$ID %in% remain)
rm(control,controlsig,result,resultsig,bothsig,remain)

#Check the constitution in sig
resultsig<-subset(result1,result1$padj<0.1) #202
nrow(subset(resultsig,resultsig$type %in% c("Unclassified","No description","DNase-seq","Input control")))#34
nrow(subset(resultsig,resultsig$type == "Histone"))#57
nrow(subset(resultsig,resultsig$type == "TFs and others"))#58
nrow(subset(resultsig,resultsig$type == "RNA polymerase"))#53

#Check the constitution in all
nrow(subset(result1,result1$type %in% c("Unclassified","No description","DNase-seq","Input control")))#1539
nrow(subset(result1,result1$type == "Histone"))#1009
nrow(subset(result1,result1$type == "TFs and others"))#1207
nrow(subset(result1,result1$type == "RNA polymerase"))#276

#plot the distribution of psig and pajd
descdist(result1$psig, discrete = FALSE)
fit_uniform<-fitdistrplus::fitdist(result1$psig,"unif")
plot(fit_uniform)

###################################################
#HISTONE FIRST
##################################################
histone<-subset(result1,result1$type == "Histone")
histone$object<-factor(histone$object)

positivecode<-c("H2BK5ac", "H2Bubi", "H2Bubi", "H3K18ac", "H3K23ac", "H3K27ac", "H3K27me1", "H3K36me2", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3",
                "H3K79me1", "H3K79me2", "H3K79me3", "H3K9ac", "H3K9acS10P", "H4acTetra", "H4K12ac", "H4K16ac", "H4K20me1", "H4K5ac", "H4K8ac","H2Bub","H3ac","H3K9K14ac","H4ac")

negativecode<-c("H2Av", "H3K27me2", "H3K27me3", "H3K36me1", "H3K9me1", "H3K9me2", "H3K9me3","H2A.V","H2AK118ub","H2AK119ub1")

transcriptionside<-vector()

for (i in 1:nrow(histone)){
  if (histone$object[i] %in% positivecode){
    transcriptionside[i]<-"up"
  } else if (histone$object[i] %in% negativecode){
    transcriptionside[i]<-"repress"
  } else {
    transcriptionside[i]<-"others"
  }
}

histone1<-cbind(histone, transcriptionside)

nrow(subset(histone1,histone1$transcriptionside =="up"))#480 activation histone in all
nrow(subset(histone1,histone1$transcriptionside =="repress"))#491 repression histone in all
nrow(subset(histone1,histone1$transcriptionside =="others")) #38 others histone in all

sighistone<-subset(histone1,histone1$padj<0.1)

nrow(subset(sighistone,sighistone$transcriptionside =="up"))#55 activation histone in padj<0.1
nrow(subset(sighistone,sighistone$transcriptionside =="repress"))#2 repression histone in padj<0.1
nrow(subset(sighistone,sighistone$transcriptionside =="others")) #0 others histone in padj<0.1

#compare distribution of sig vs all
histonemat <- matrix(c(55,480,2, 491, 0, 38), ncol=3, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("activation", "repression", "others")))
fisher.test(histonemat) # p-value = 6.174e-14
chisq.test(histonemat) # X-squared = 51.657, df = 2, p-value = 6.063e-12, Warning message: In chisq.test(histonemat) : Chi-squared近似算法有可能不准, (due to one predicted observed<5)

################################################
#HISTONE PLOT
################################################

histonetable<-data.frame(marks<-c(rep(c("Activation","Repression", "Others"),2)),
                            padj<-c(rep("Significant",3),rep("All",3)),
                          value<-c(96.5,3.5,0,47.6,48.7,3.7))

colnames(histonetable)<-c("Functions","padj","value")

#Create Another matrix for significance
histoneper<-data.frame(as.character(padj<-c("Significant")), value<-c(as.numeric(102))) #add the P on top of the bar, may need to change 10 depend on the figure
colnames(histoneper)<-c("padj","value")

#Create a third matrix to show the percentage of each percentage and the total number.
histonenumber<-data.frame(padj<-c(rep("Significant",3),rep("All",3)),value<-c(3.5+0.5*96.5,0.5*3.5,0, 48.7+3.7+47.6*0.5, 3.7*0.5+48.7, 48.7*0.5))#sigall #sigup  #sigother #sig down #all count #allup #all other #all down
colnames(histonenumber)<-c("padj","value")

#plot the dataset
p=ggplot(histonetable, aes(padj,value))+ 
  geom_bar(stat="identity", size=0.1, position="stack",width = 0.4, aes(fill=Functions))+ #careful here, don't put fill on ggplot if want to combine with other geom_text
  ylab("Percentage of data set")+xlab("Overlapping significance")+
  scale_y_continuous(breaks = c(20,40,60,80,100))+
  scale_fill_manual(values = c("#ff7761","#00dffc","#6d819c"))


p1=p+geom_text(data=histoneper,label="fisher exact test: p=6.174e-14", size=5, show.legend = F)

p2=p1+geom_text(data=histonenumber,label=c("96.5%","3.5%","","47.6%","3.7%","48.7%"), size=4.5, show.legend = F)


p2+theme(axis.line=element_line(size=2),axis.text=element_text(size=20),axis.title=element_text(size=20,face="plain",color="deepskyblue4"),
         legend.text = element_text(size=20), legend.title = element_text(size=20))

#remove necessary
rm(histone,histonemat,histonenumber,histonetable,otherhistone,p,p1,p2,sighistone,i,marks,negativecode,otherhistonelevels,padj,positivecode,transcriptionside,value,histoneper)

###############################################################
#NEXT WOULD BE TRANSCRIPTION FACTOR AND CHROMATIN MODIFIER
###############################################################

cm_tf<-subset(result1,result1$type == "TFs and others")

setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/report")
objectanno<-read.csv("allobjannotated2.csv") #main difference will be HDAC1/RPD3
tf<-subset(objectanno,objectanno$category %in% "tf")
cm<-subset(objectanno,objectanno$category %in% "cm")

category<-vector()

for (i in 1:nrow(cm_tf)){
  if (tolower(cm_tf$object[i]) %in% tolower(tf$object)){
    category[i]<-"tf"
  } else if (tolower(cm_tf$object[i]) %in% tolower(cm$object)){
    category[i]<-"cm"
  } else {
    category[i]<-"NA"
  }
}

cm_tf_1<-cbind(cm_tf,category)

cm_tf_na<-subset(cm_tf_1,cm_tf_1$category=="NA")
cm_tf_na$object<-factor(cm_tf_na$object)
remain<-levels(cm_tf_na$object)