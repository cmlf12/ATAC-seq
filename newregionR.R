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
#library(BSgenome.Dmelanogaster.UCSC.dm6)
library(tidyverse)
library("reader")
#library(EnrichedHeatmap)
#library("GenometriCorr")
#library(ComplexHeatmap)
library(ggplot2)
library(pheatmap)
library(ChIPseeker)
library(clusterProfiler)



##########Loading peak-called ATAC/CHip-seq data##########
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/datacollection")
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/datacollection")


collect<-dir()
filelist<-list()
filetable<-list()
filetableGrange<-list()
prerelease<-data.frame(matrix(NA, nrow=length(collect), ncol=60))
ID<-list()
Sequence<-list()
object<-list()
stage<-list()
tissue<-list()
chip<-list()
interval<-list()


for (i in 1:length(collect)){
  
  
  #get the file address
  address<-paste(getwd(), "/", collect[i], sep='')

  #get the file ID
  ID[[i]]<-str_extract(address, "modENCODE_[0-9]*")
  
  #get the object/
  
  object[[i]]<-str_extract(collect[i], "^[\\s\\S]*?_")
  
  #get the stage
  stage[[i]]<-str_extract(collect[i], "Developmental-Stage=[\\s\\S]*?_")
  
  #get the tissue type
  tissue[[i]]<-str_extract(collect[i], "Tissue=[\\s\\S]*?_")
  
  #chip-chip or chip-seq
  chip[[i]]<-str_extract(collect[i], "ChIP-[\\s\\S]*_")
  

  try({ 
    
    #read the gff file as data.frame
    filetable[[i]]<-read.table(address,  stringsAsFactors = FALSE, fill=TRUE,skipNul = TRUE)
    
    if (ncol(filetable[[i]])==9){
      
      colnames(filetable[[i]])<-c("chr", "source", "type","starter","ender","score","strander","phase","attribute")
      
    }
    
    else if (ncol(filetable[[i]])==8){
      
      colnames(filetable[[i]])<-c("chr", "source", "type","starter","ender","score","strander","attribute")
      
    }
    
    else{
      colnames(filetable[[i]])<-c("chr", "source", "type","starter","ender","score","attribute")
    }
    
    
    #rename the chromosome to match the chain file
    c<-filetable[[i]]$chr
    chr<-sapply(c, function(x) paste("chr", x, sep=''))
    filetable[[i]]<-filetable[[i]][,-1]
    filetable[[i]]<-cbind(chr,filetable[[i]])
    
    #identify the genome release
    prerelease[i,]<-readLines(address, n=60, skipNul = TRUE, ok=TRUE)
    
    #convert to grangers file
    filetableGrange[[i]]<-makeGRangesFromDataFrame(filetable[[i]],  keep.extra.columns = TRUE, seqnames.field = "chr", start.field = "starter", end.field = "ender", strand.field = "strander")
    
    #get the numbers of interval
    interval[[i]]<-length(filetableGrange[[i]])
    
    
 }, silent=FALSE, outFile=stdout())
  
  cat(i, "\n")
  
}

define<-data.frame(sapply(prerelease, function(x) str_extract(x, "r[0-9]")))

minus1<-function(x){
  x<-sub("_", "", x)
}


convertline<-function(x){
  x<-as.vector(unlist(x))
  x<-sapply(x, function(x) minus1(x))
  unname(x)
}

object1<-convertline(object)
tissue1<-convertline(tissue)
stage1<-convertline(stage)
chip1<-convertline(chip)
ID<-as.vector(unlist(ID))
interval1<-as.vector(unlist(interval))

k<-cbind(object1,tissue1,stage1,chip1, interval1)


#####identify the range################
a<-list()
for (i in 1:nrow(define)){
  a[[i]]<-define[i,][!is.na(define[i,])]
}

carry<-function(x){
  x[1]
}

releasefinal<-data.frame(matrix(NA, nrow=1248, ncol=2))

i<-1
for (i in 1:length(a)){
  b<-carry(a[[i]])
  releasefinal[i,]<-cbind(ID[[i]], b)
}

problemrelease<-subset(releasefinal, releasefinal$X2 == "r5")
NArelease<-releasefinal[which(is.na(releasefinal$X2)),]


#############identify the foxo data###########################################
setwd("~/OneDrive/uclPHD/meta-data/experiment/foxo/results of coordination")

#foxoupdated<-read.csv("fbfoxo coordination.csv")
#foxo<-read.csv("FBallResults.csv")
#foxoupdated1<-foxoupdated[-c(1:2),]
#colnames(foxoupdated1)<-c("rownames.Res.", "chr", "conversion results" )
#foxo1<-left_join(foxo, foxoupdated1, by="rownames.Res.")
#write.csv(foxo1,file=paste(getwd(), "/", "foxofb.csv", sep=''))

#############convert foxo to formal expression##################
foxo<-read.csv("foxofb.csv")
chromosome<-sapply(foxo$chr, function(x) str_extract(x, "^.*:"))
start<-sapply(foxo$chr, function(x) str_extract(x, ":[0-9]+.."))
end<-sapply(foxo$chr, function(x) str_extract(x, "..[0-9]+$"))
chromosome<-sapply(chromosome, function(x) str_remove(x, ":"))
chromosome<-sapply(chromosome, function(x) paste("chr", x, sep=''))
start<-sapply(start, function(x) str_remove(x, ":"))
start<-sapply(start, function(x) str_remove(x, "..$"))
end<-sapply(end, function(x) str_remove(x, ".."))
a<-data.frame(cbind(chromosome, start, end))
foxo1<-cbind(foxo,a)
#write.csv(foxo1,file=paste(getwd(), "/", "foxofb1.csv", sep=''))

export.gff3(foxogRange, paste(getwd(), "/", ID[[i]], "foxosig.gff3", sep=''))


#############load genomer5#########################################
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/datacollection")

collect<-dir()
pt2<-list()
dmgenome<-characterToBSGenome("dm3")

########convert foxo to gRange list###########################

setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/foxo/results of coordination")
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/foxo/results of coordination")

foxo<-read.csv("foxofb1.csv")

foxonoNA<-foxo[complete.cases(foxo),]
foxopadj<-subset(foxo, foxo$padj<0.1)
foxopadjnoNA<-foxopadj[complete.cases(foxopadj),]
foxopadjnoNA$chromosome<-droplevels(foxopadjnoNA$chromosome)
#foxoopen<-subset(foxopadjnoNA, foxopadjnoNA$log2FoldChange>0)

foxogRange<-makeGRangesFromDataFrame(foxopadjnoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)
#foxofullgRange<-makeGRangesFromDataFrame(foxonoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)
#foxoopengRange<-makeGRangesFromDataFrame(foxoopen,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

#export.gff3(foxogRange, paste(getwd(), "/", ID[[i]], "foxosig.gff3", sep=''))

rm(foxo,foxoNA,foxopadj,foxopadjnoNA,foxonoNA)
#######################################permutation test##################################
i<-1 

for (i in 1:length(collect)){
  
  
  try(pt2[[i]]<-overlapPermTest(A=foxofullgRange, B=filetableGrange[[i]], ntime=1500,genome=dmgenome, alternative = "auto"))
  
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

ID<-as.vector(unlist(ID))
psig<-as.vector(unlist(psig))
observed<-as.vector(unlist(observed))
Zscore<-as.vector(unlist(Zscore))

ID.num<-as.vector(as.character(1:length(ID)))


padj<-p.adjust(psig, method="BH")


result<-as.data.frame(cbind(ID, psig, padj, observed, Zscore,k))


##############output###########################################
setwd("~/OneDrive/uclPHD/meta-data/experiment")

write.csv(result, file=paste(getwd(), "/", "10000_onlyup_permtest.csv", sep=''))
##########################remove all unnecessary thing for further analysis####
rm(chip, define,filelist,filetable,interval,object,prerelease,Sequence,stage,tissue,address,c,chip1,chr,collect,i,interval1,stage1,tissue1,convertline,minus1)


##########load it back###############################################
setwd("~/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap")
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap") 

result<-read.csv(paste(getwd(), "/", "10000_onlyup_aftercontrol.csv", sep='')) #1077 object
#control<-read.csv(paste(getwd(), "/", "1000_onlyupcontrol_permtest.csv", sep=''))

###########plot something############################################
plot(result$padj)
plot(result$psig)

##########select against control###################################
upcontrolresult$padj<-as.numeric(as.character(upcontrolresult$padj))
controlsig<-subset(upcontrolresult, upcontrolresult$padj<0.1)
foxosig<-subset(result, result$padj<0.01)

normalize<-intersect(foxosig$ID, controlsig$ID)

foxoaftercontrol<-subset(result, !result$ID %in% normalize)



##############################################################
##########separate the histone and the other 
##############################################################
histone<-filter(result, result$object1 %in% str_extract(result$object1, "H[0-9][\\s\\S]*"))

positivecode<-c("H2BK5ac", "H2Bubi", "H2Bubi", "H3K18ac", "H3K23ac", "H3K27ac", "H3K27me1", "H3K36me2", "H3K36me3", "H3K4me1", "H3K4me2", "H3K4me3",
"H3K79me1", "H3K79me2", "H3K79me3", "H3K9ac", "H3K9acS10P", "H4acTetra", "H4K12ac", "H4K16ac", "H4K20me1", "H4K5ac", "H4K8ac")

negativecode<-c("H2Av", "H3K27me2", "H3K27me3", "H3K36me1", "H3K9me1", "H3K9me2", "H3K9me3")

transcriptionside<-vector()

for (i in 1:nrow(histone)){
  if (histone$object1[i] %in% positivecode){
    transcriptionside[i]<-"up"
  } else if (histone$object1[i] %in% negativecode){
    transcriptionside[i]<-"repress"
  } else {
    transcriptionside[i]<-"NA"
  }
}

histone1<-cbind(histone, transcriptionside)

sighistone<-subset(histone1, histone1$padj<0.01)
nosighistone<-subset(histone1, histone1$padj>0.01)

#################################################################
##Seperate tf and and chromatin modifier
#################################################################
setwd("C:/Users/cmlf/OneDrive/uclPHD/meta-data/switch_experiment/atac_overlap/report")
objectanno<-read.csv("allobjannotated2.csv") #main difference will be HDAC1/RPD3
tf<-subset(objectanno,objectanno$category %in% "tf")
cm<-subset(objectanno,objectanno$category %in% "cm")

nonhistone<-subset(result,result$X %in% setdiff(result$X, histone$X))
category<-vector()

for (i in 1:nrow(nonhistone)){
  if (nonhistone$object1[i] %in% tf$object){
    category[i]<-"tf"
  } else if (nonhistone$object1[i] %in% cm$object){
    category[i]<-"cm"
  } else {
    category[i]<-"NA"
  }
}

nonhistone1<-cbind(nonhistone,category)

tfresult<-subset(nonhistone1,nonhistone1$category %in% "tf")
cmresult<-subset(nonhistone1,nonhistone1$category %in% "cm")

trx<-as.character(subset(cm$object,cm$Trx_polycomb %in% "TRX"))
polycomb<-as.character(subset(cm$object,cm$Trx_polycomb %in% c("polycomb","Polycomb")))
other<-as.character(subset(cm$object,cm$Trx_polycomb %in% "N"))
group<-vector()

for (i in 1:nrow(cmresult)){
  if (cmresult$object1[i] %in% trx){
    group[i]<-"trx"
  } else if (cmresult$object1[i] %in% polycomb){
    group[i]<-"polycomb"
  } else {
    group[i]<-"others"
  }
}

cmresult1<-cbind(cmresult,group)
sigcmresult<-subset(cmresult1,cmresult$padj<0.01)


################################################################
#TEST HISTONE SIGNIFICANT
################################################################

n<-nrow(histone1) #Number of groups
us<-nrow(subset(sighistone, sighistone$transcriptionside == "up")) #Number of groups padj<0.01 and have activation marks #86
ds<-nrow(subset(sighistone, sighistone$transcriptionside == "repress")) #Number of groups padj<0.01 and have repression marks #16
ns<-nrow(subset(sighistone, sighistone$transcriptionside == "NA")) #Number of groups padj<0.01 and have mild marks #9
uns<-nrow(subset(nosighistone, nosighistone$transcriptionside == "up")) #Number of groups padj>0.01 and have activation marks #268
dns<-nrow(subset(nosighistone, nosighistone$transcriptionside == "repress")) #Number of groups padj>0.01 and have repression marks #106
nns<-nrow(subset(nosighistone, nosighistone$transcriptionside == "NA")) #Number of groups padj>0.01 and have mild marks #27
uas<-us+uns #354
das<-ds+dns #122
nas<-ns+nns #36
aas<-uas+das+nas #512
as<-us+ds+ns #111

  
#compare distribution of nonsig vs sig
histonemat <- matrix(c(us,uns, ds, dns, ns, nns), ncol=3, nrow=2,
              dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("activation", "repression", "neutral")))
fisher.test(histonemat) # p-value = 0.02418
chisq.test(histonemat) #X-squared = 6.9293, df = 2, p-value = 0.03128

#compare distribution of activation of sig vs nonsig
histonemat1 <- matrix(c(us,uns, ds+ns, dns+nns), ncol=2, nrow=2,
                      dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("activation", "non-activation")))
fisher.test(histonemat1) # p-value = 0.0364
chisq.test(histonemat1) #X-squared = 4.1313, df = 1, p-value = 0.0421

#compare distribution of repression of sig vs nonsig
histonemat2<- matrix(c(ds,dns, us+ns, uns+nns), ncol=2, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("repressio", "non-repressio")))
fisher.test(histonemat2) # p-value = 0.008114
chisq.test(histonemat2) #X-squared = 6.2733, df = 1, p-value = 0.01226


#compare distribution of sig vs all
histonemat <- matrix(c(us,uas, ds, das, ns, nas), ncol=3, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("activation", "repression", "neutral")))
fisher.test(histonemat) # p-value = 0.08261
chisq.test(histonemat) # X-squared = 4.6937, df = 2, p-value = 0.09567

#compare distribution of activation of sig vs all
histonemat1 <- matrix(c(us,uas, ds+ns, das+nas), ncol=2, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("activation", "non-activation")))
fisher.test(histonemat1) # p-value = 0.08548
chisq.test(histonemat1) #X-squared = 2.6675, df = 1, p-value = 0.1024

#compare distribution of repression of sig vs all
histonemat2<- matrix(c(ds,das, us+ns, uas+nas), ncol=2, nrow=2,
                      dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("repressio", "non-repressio")))
fisher.test(histonemat2) # p-value = 0.03197
chisq.test(histonemat2) #X-squared = 4.1579, df = 1, p-value = 0.04144

#compare distribution of others of sig vs all
histonemat3<- matrix(c(ns,nas, us+ds, uas+das), ncol=2, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), Histone_marks = c("others", "non-others")))
fisher.test(histonemat3) # p-value = 0.6867
chisq.test(histonemat3) #X-squared = 0.038058, df = 1, p-value = 0.8453


###################################################################################################
#PLOT THE HISTONE NUMBER
###################################################################################################
#Create matrix to plot the data#
histonetable<-data.frame(marks<-c(rep("activation", 2),rep("repression",2), rep("others",2)),
                         padj<-rep(c("<0.01","all"),3),
                         value<-c(us,uas, ds, das, ns, nas))
colnames(histonetable)<-c("functions","padj","value")

#Create Another matrix for significance
histoneper<-data.frame(as.character(padj<-c("<0.01")), value<-c(as.numeric(as+50))) #add the P on top of the bar, may need to change 10 depend on the figure
colnames(histoneper)<-c("padj","value")

#Create a third matrix to show the percentage of each percentage and the total number.
histonenumber<-data.frame(padj<-c(rep("<0.01",4),rep("all",4)),value<-c(as+20,as-0.5*us, ds+0.5*ns, 0.5*ds,401+20,401-0.5*uns, dns+0.5*nns, 0.5*dns))#sigall #sigup  #sigother #sig down #all count #allup #all other #all down
colnames(histonenumber)<-c("padj","value")

#plot the dataset
p=ggplot(histonetable, aes(padj,value))+ 
  geom_bar(stat="identity", size=0.1, position="stack",width = 0.4, aes(fill=functions))+ #careful here, don't put fill on ggplot if want to combine with other geom_text
  ylab("Number of data set")+xlab("Overlapping significance")


p1=p+geom_text(data=histoneper,label="Fisher exact test: p=0.085", size=5, show.legend = F)

p2=p1+geom_text(data=histonenumber,label=c("111 data set","77.5%","8.1%","14.4%","401 data set","66.8%","6.7%","26.5%"), size=4.5, show.legend = F)


p2+theme(axis.line=element_line(size=2),axis.text=element_text(size=20),axis.title=element_text(size=20,face="plain",color="deepskyblue4"),
        legend.text = element_text(size=20), legend.title = element_text(size=20))

###################################################################################################
#PLOT THE HISTONE PERCENTAGE
###################################################################################################
#Create matrix to plot the data

histonetable2<-data.frame(marks<-c(rep(c("Activation","Other","Repression"),2)),
                         padj<-c(rep("Significant",3),rep("Unsignificant",3)),
                         value<-c(77.5,8.1,14.4,66.8,6.7,26.5))

colnames(histonetable2)<-c("Functions","padj","value")

#Create Another matrix for significance
histoneper2<-data.frame(as.character(padj<-c("Significant")), value<-c(as.numeric(102))) #add the P on top of the bar, may need to change 10 depend on the figure
colnames(histoneper2)<-c("padj","value")

#Create a third matrix to show the percentage of each percentage and the total number.
histonenumber2<-data.frame(padj<-c(rep("Significant",3),rep("Unsignificant",3)),value<-c(8.1+14.4+0.5*77.5,14.4+0.5*8.1, 14.4*0.5, 6.7+26.5+66.8*0.5, 26.5+0.5*6.7, 26.5*0.5))#sigall #sigup  #sigother #sig down #all count #allup #all other #all down
colnames(histonenumber2)<-c("padj","value")

#plot the dataset
p=ggplot(histonetable2, aes(padj,value))+ 
  geom_bar(stat="identity", size=0.1, position="stack",width = 0.4, aes(fill=Functions))+ #careful here, don't put fill on ggplot if want to combine with other geom_text
  ylab("Percentage of data set")+xlab("Overlapping significance")+
  scale_y_continuous(breaks = c(20,40,60,80,100))+
  scale_fill_manual(values = c("#ff7761","#00dffc","#6d819c"))


p1=p+geom_text(data=histoneper2,label="Pearson's Chi-squared test: p=0.0313", size=5, show.legend = F)

p2=p1+geom_text(data=histonenumber2,label=c("77.5%","8.1%","14.4%","66.8%","6.7%","26.5%"), size=4.5, show.legend = F)


p2+theme(axis.line=element_line(size=2),axis.text=element_text(size=20),axis.title=element_text(size=20,face="plain",color="deepskyblue4"),
         legend.text = element_text(size=20), legend.title = element_text(size=20))

###############################################################################
#TRX AND POLYCOMB
###############################################################################
#1.Use self annotated result
st<-nrow(subset(sigcmresult,sigcmresult$group == "trx")) #9
sp<-nrow(subset(sigcmresult,sigcmresult$group == "polycomb"))#13
so<-nrow(subset(sigcmresult,sigcmresult$group == "others"))#26
at<-nrow(subset(cmresult1,cmresult1$group == "trx")) #79
ap<-nrow(subset(cmresult1,cmresult1$group == "polycomb")) #67
ao<-nrow(subset(cmresult1,cmresult1$group == "others")) #225

#compare distribution of nonsig vs sig
cmmatrix <- matrix(c(st,at, sp, ap, so, ao), ncol=3, nrow=2,
                     dimnames = list(significant = c("significant", "unsignificant"), chromatin_group = c("trx", "polycomb", "others")))
fisher.test(cmmatrix) # p-value = 0.5768
chisq.test(cmmatrix) #X-squared = 0.95792, df = 2, p-value = 0.6194

#2.Use the define flybase Gene groups
st<-nrow(subset(sigcmresult,sigcmresult$group == "trx")) #17
sp<-nrow(subset(sigcmresult,sigcmresult$group == "polycomb"))#5
so<-nrow(subset(sigcmresult,sigcmresult$group == "others"))#26
at<-nrow(subset(cmresult1,cmresult1$group == "trx")) #70
ap<-nrow(subset(cmresult1,cmresult1$group == "polycomb")) #50
ao<-nrow(subset(cmresult1,cmresult1$group == "others")) #251

#compare distribution of nonsig vs sig
cmmatrix <- matrix(c(st,at, sp, ap, so, ao), ncol=3, nrow=2,
                   dimnames = list(significant = c("significant", "unsignificant"), chromatin_group = c("trx", "polycomb", "others")))
fisher.test(cmmatrix) # p-value = 0.03914
chisq.test(cmmatrix) #X-squared = 7.0785, df = 2, p-value = 0.02903

#3.Use the define flybase Gene groups,and deleted HDAC1
st<-nrow(subset(sigcmresult,sigcmresult$group == "trx"))-8 #9
sp<-nrow(subset(sigcmresult,sigcmresult$group == "polycomb"))#5
so<-nrow(subset(sigcmresult,sigcmresult$group == "others"))#26
at<-nrow(subset(cmresult1,cmresult1$group == "trx"))-17 #53
ap<-nrow(subset(cmresult1,cmresult1$group == "polycomb")) #50
ao<-nrow(subset(cmresult1,cmresult1$group == "others")) #251

cmmatrix <- matrix(c(st,at, sp, ap, so, ao), ncol=3, nrow=2,
                   dimnames = list(significant = c("significant", "unsignificant"), chromatin_group = c("trx", "polycomb", "others")))
fisher.test(cmmatrix) # p-value = 0.4447
chisq.test(cmmatrix) #X-squared = 1.5405, df = 2, p-value = 0.4629

cmmatrix <- matrix(c(st,at, sp+so, so+ao), ncol=2, nrow=2, #compare trx with other in sig vs all
                   dimnames = list(significant = c("significant", "unsignificant"), chromatin_group = c("trx", "others")))
fisher.test(cmmatrix) # p-value = 0.3678
chisq.test(cmmatrix) #X-squared = 0.64913, df = 1, p-value = 0.4204

########################################################################
#PLOT TRX/POLYCOMB PERCENTAGE
#########################################################################
cmtable<-data.frame(marks<-c(rep(c("Trx","Others","Polycomb"),2)),
                          padj<-c(rep("Significant",3),rep("All",3)),
                          value<-c(35.4,54.2,10.4,18.9,67.4,13.7))

colnames(cmtable)<-c("Functions","padj","value")

#Create Another matrix for significance
cmper<-data.frame(as.character(padj<-c("Significant")), value<-c(as.numeric(102))) #add the P on top of the bar, may need to change 10 depend on the figure
colnames(cmper)<-c("padj","value")

#Create a third matrix to show the percentage of each percentage and the total number.
cmnumber<-data.frame(padj<-c(rep("Significant",3),rep("All",3)),value<-c(10.4+0.5*54.2+35.4,0.5*10.4+35.4, 35.4*0.5, 0.5*67.4+13.7+18.9, 18.9+0.5*13.7, 18.9*0.5))#sigall #sigup  #sigother #sig down #all count #allup #all other #all down
colnames(cmnumber)<-c("padj","value")

#plot the dataset
p=ggplot(cmtable, aes(padj,value))+ 
  geom_bar(stat="identity", size=0.1, position="stack",width = 0.4, aes(fill=Functions))+ #careful here, don't put fill on ggplot if want to combine with other geom_text
  ylab("Percentage of data set")+xlab("Overlapping significance")+
  scale_y_continuous(breaks = c(20,40,60,80,100))+
  scale_fill_manual(values = c("#ff7761","#00dffc","#6d819c"))


p1=p+geom_text(data=cmper,label="Pearson's Chi-squared test: p=0.02903", size=5, show.legend = F)

p2=p1+geom_text(data=cmnumber,label=c("54.2%","10.4%","35.4%","67.4%","13.7%","18.9%"), size=4.5, show.legend = F)


p2+theme(axis.line=element_line(size=2),axis.text=element_text(size=20),axis.title=element_text(size=20,face="plain",color="deepskyblue4"),
         legend.text = element_text(size=20), legend.title = element_text(size=20))

#################################################################################
#deep into the overlap pattern
###############################################################################
#get the list of object which is significant after control#############
result1<-result

sigresult<-subset(result1, result1$padj<0.01)

whichgrange<-as.vector(as.numeric(sigresult$X))

#############################################################
#get the list of interesting object
#############################################################
########select the specific complex component to plot################

#complex<-c("nejire","ISWI","ACF1","NURF301","CtBP","MLE","MOF","MSL-1","HDAC1","dMi-2","Pho","dSFMBT","WDS","JHDM1","Psc","dRING","PIWI","Rhino","SNR1","brahma","Enhancer-of-zeste","Pc","Su(var)3-9","HP1a","JIL-1","LSD1","HDAC1","Su(var)3-7") #ISWI,MSL,NURD,PHORC,

#lets find the trx group proteins!
complex<-c("nejire","ISWI","ACF1","NURF301","CtBP","SNR1","brahma","HDAC1","dMi-2","WDS") #iswi, #INO80, #SWI/SNF #NURD #Compass/TRR/TRX(only WDS), #TAC1.#NURF #ASH1
#Then find the polycomb group proteins
complex<-c("JHDM1","Psc","dRING","Pho","dSFMBT","Enhancer-of-zeste","Pc")

component<-subset(result1, result1$object1 %in% complex) 
whichcomponent<-as.vector(as.numeric(component$X))

#############run the findOverlaps to find the overlap in foxo ###############
overlap<-list()
ID2<-vector()
object2<-vector()

for (i in 1:length(whichcomponent)){           #choose either whichgrange or whichcomponent, and change three position
  num<-whichcomponent[i]
  overlapped<-findOverlaps(foxogRange, filetableGrange[[num]])
  overlap[[i]]<-as.vector(as.character(overlapped@from))
  ID2[i]<-ID[num]
  object2[i]<-object1[num]
}

overlapdata<-data.frame(matrix(NA, nrow=length(whichcomponent), ncol=length(foxogRange)))
a<-as.vector(rep(0, times=length(foxogRange)))


for (i in 1:length(overlap)){
  overlapdata[i,]<-replace(a, as.numeric(overlap[[i]]), 1)
}

colnames(overlapdata)<-c(1:length(foxogRange))
overlapdata<-cbind(ID2,object2,overlapdata) 

overlapdata1<-overlapdata[order(pmatch(overlapdata$object2,complex, duplicates.ok = TRUE)),, drop=FALSE]


#########################################################################################
##############Use chipseeker to explore the feature of foxo sig peak
###plot the foxo sig region#############################
covplot(foxogRange)


##input the genome annotation file#####################
setwd("~/OneDrive/uclPHD/meta-data/drosophila genome reference")
degnomeannomation <-import("dmel-all-r5.9.gff")
txdb<-makeTxDbFromGRanges(degnomeannomation)

foxogRange1<-foxogRange
seqlevelsStyle(foxogRange1)<-"Ensembl"


########foxo peak with promoter################
promoter<-getPromoters(TxDb=txdb, upstream = 3000, downstream = 3000 )
promotermatrix<-getTagMatrix(foxogRange, windows=promoter)

######heatmap of the promoter
tagHeatmap(promotermatrix, xlim=c(-3000,3000))

plotAvgProf(promotermatrix, xlim=c(-3000,3000), xlab="Genomic region (5-3)", ylab="Read count frequency")

#####do some annotaiton#######
peakANNO<-annotatePeak(foxogRange1,tssRegion=c(-3000,3000), TxDb=txdb)
foxopeakid<-peakANNO@anno$geneId
peakannoresult<-as.data.frame(as.GRanges(peakANNO))
#write.csv(peakannoresult, file=paste(getwd(), "/", "foxosigpeak.csv", sep=''))

######plot the annotation #####################
plotAnnoPie(peakANNO)


upsetplot(peakANNO, vennpie=TRUE)


#############################################################################
#specify some condition on object and plot heatmap
############################################################################

forheatmap<-as.matrix(overlapdata1[,-c(1,2)])
rownames(forheatmap)<-overlapdata1$object2
colnames(forheatmap)<-foxopeakid


##############################################################
#ggplot for heatmap
##############################################################
#reshape the data

ggheatmap<-data.frame(matrix(NA, nrow=12705, ncol=3))

for (i in 1:nrow(forheatmap)){
  
  objectname<-rownames(forheatmap)[i]
  
  for (k in 1:ncol(forheatmap)){
    
    position<-colnames(forheatmap)[k]
    
    value<-forheatmap[i,k]
    
    z=((as.numeric(i)-1)*77+as.numeric(k))   ##attention the number
    
    ggheatmap[z,]<-cbind(objectname,position,value)
    
    k<-k+1
  }
  
  i<-i+1
  
}

colnames(ggheatmap)<-c("object", "foxopeak", "overlap")
cols<-c("1"="violetred2", "0"="gray70")

#heatmap by ggplot
ggplot(data = ggheatmap, aes(x = object, y = foxopeak)) +
geom_tile(aes(fill = overlap), colour="steelblue")+
scale_fill_manual(values=cols)




##################################################
#heatmap by pheatmap

pheatmap(forheatmap, color=c( "gray70","violetred2"), fontsize=9, fontsize_row=6, cellwidth = 7, cellheight = 5, cluster_cols=TRUE, cluster_rows=FALSE)

