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
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(tidyverse)
library("reader")
library(EnrichedHeatmap)
library("GenometriCorr")
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(ggplot2)
library(pheatmap)
library(ChIPseeker)
library(clusterProfiler)



##########Loading peak-called ATAC/CHip-seq data##########
setwd("~/OneDrive/uclPHD/meta-data/experiment/datacollection")

collect<-dir()
filelist<-list()
filetable<-list()
filetableGrange<-list()
prerelease<-data.frame(matrix(NA, nrow=length(collect), ncol=60))
filetowork<-list()
sr<-list()
mr<-list()
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


###########gather the gff file to a new folder##########################
#setwd("~/OneDrive/uclPHD/meta-data")
#if(!dir.exists("gff3tobe")){
# dir.create("gff3tobe")
#}

#mydir<-
#for(file in filelist) {
# See ?file.copy for more options
# file.copy(file, "gff3tobe")
#}

#setwd("~/OneDrive/uclPHD/meta-data/experiment")


#############identify the foxo data###########################################
setwd("~/OneDrive/uclPHD/meta-data/experiment/foxo/results of coordination")

#foxoupdated<-read.csv("fbfoxo coordination.csv")
#foxo<-read.csv("FBallResults.csv")
#foxoupdated1<-foxoupdated[-c(1:2),]
#colnames(foxoupdated1)<-c("rownames.Res.", "chr", "conversion results" )
#foxo1<-left_join(foxo, foxoupdated1, by="rownames.Res.")
#write.csv(foxo1,file=paste(getwd(), "/", "foxofb.csv", sep=''))



######################convert foxo to formal expression##################
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



#############convert the release version to the newest #####################
chain<-import.chain("~/OneDrive/uclPHD/meta-data/experiment/chain/dm3ToDm6.over.chain")
updatedfile<-lapply(filetableGrange, function(x) liftOver(x, chain))
updatedfile.df<-lapply(updatedfile, as.data.frame)

##############output the file################################
setwd("~/OneDrive/uclPHD/meta-data/experiment/updated")

i<-1

for (i in 1:length(updatedfile.df)){
  export.gff3(updatedfile.df[[i]], paste(getwd(), "/", ID[[i]], "updated.gff3", sep=''))
}


##################loading the updated gff file back###########
setwd("~/OneDrive/uclPHD/meta-data/experiment/updated")

#get the file address

collect<-dir()
i<-1

for (i in 1:length(collect)){
  filetowork[[i]]<-import(collect[i])
  names(filetowork[[i]])<-collect[i]
}



###############check whether sequence identical after converting####################################################
i<-1

for (i in 1:length(filelist)){
  
  Sequence[[i]]<-BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6, filetowork[[i]])
  
}

BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, filetableGrange[[i]])

BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, names=filetable[[i]]$chr, start=filetable[[i]]$starter, end=filetable[[i]]$ender, strand=filetable[[i]]$strander, as.character=TRUE)





#################check the updated, and if anything missing#########
i<-1

for (i in 1:length(collect)){
  
  #Find site split to multiple regions
  
  idd<-data.frame(table(filetowork[[i]]$attribute))
  sr[[i]]<-filetowork[[i]][filetowork[[i]]$attribute %in% idd$Var1[idd$Freq > 1],]
  
  
  #Find merge split
  mrid<-setdiff(filetableGrange[[i]]$attribute, filetowork[[i]]$attribute)
  mr[[i]]<-filetableGrange[[i]][filetableGrange[[i]]$attribute %in% mrid,]
  
}

##############make a tabel for number of split and merge#####################
split<-sapply(sr, function(x) length(x))
merge<-sapply(mr,function(x) length(x))
ID<-1:length(collect)

smr<-data.frame(ID, split, merge)








##################basic#################################################
setwd("~/OneDrive/uclPHD/meta-data/experiment/datacollection")
collect<-dir()
pt2<-list()
dmgenome<-characterToBSGenome("dm3")

################get file#############################################
######################convert foxo to gRange list)

setwd("~/OneDrive/uclPHD/meta-data/experiment/foxo/results of coordination")

foxo<-read.csv("foxofb1.csv")

foxonoNA<-foxo[complete.cases(foxo),]

foxopadj<-subset(foxo, foxo$padj<0.1)

foxopadjnoNA<-foxopadj[complete.cases(foxopadj),]

foxopadjnoNA$chromosome<-droplevels(foxopadjnoNA$chromosome)

foxoopen<-subset(foxopadjnoNA, foxopadjnoNA$log2FoldChange>0)
foxoclose<-subset(foxopadjnoNA, foxopadjnoNA$log2FoldChange<0)
foxoallopen<-subset(foxonoNA, foxonoNA$log2FoldChange>0)
foxoallclose<-subset(foxonoNA, foxonoNA$log2FoldChange<0)


foxogRange<-makeGRangesFromDataFrame(foxopadjnoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxofullgRange<-makeGRangesFromDataFrame(foxonoNA,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxoopengRange<-makeGRangesFromDataFrame(foxoopen,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxoclosegRange<-makeGRangesFromDataFrame(foxoclose,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxofullopengRange<-makeGRangesFromDataFrame(foxoallopen,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

foxofullclosegRange<-makeGRangesFromDataFrame(foxoallclose,  keep.extra.columns = TRUE, seqnames.field = "chromosome", start.field = "start", end.field = "end", ignore.strand=TRUE)

#export.gff3(foxogRange, paste(getwd(), "/", ID[[i]], "foxosig.gff3", sep=''))





#######################################permutation test##################################
i<-1 

for (i in 1:length(collect)){
  
  
  try(pt2[[i]]<-overlapPermTest(A=foxofullclosegRange, B=filetableGrange[[i]], ntime=1000,genome=dmgenome, alternative = "auto"))
  
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

write.csv(result, file=paste(getwd(), "/", "1000_onlyupaftercontrol_permtest.csv", sep=''))



##########load it back###############################################
setwd("~/OneDrive/uclPHD/meta-data/experiment")

result<-read.csv(paste(getwd(), "/", "1000_onlyupaftercontrol_permtest.csv", sep=''))


###########plot something############################################
plot(result$padj)
plot(result$psig)

##########select against control###################################
upcontrolresult$padj<-as.numeric(as.character(upcontrolresult$padj))
controlsig<-subset(upcontrolresult, upcontrolresult$padj<0.1)
foxosig<-subset(result, result$padj<0.01)

normalize<-intersect(foxosig$ID, controlsig$ID)

foxoaftercontrol<-subset(result, !result$ID %in% normalize)

##########separate the histone and the other ########################
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


nonhistone<-setdiff(result, histone)





#################################################################################
#################deep into the overlap pattern#################################
###############################################################################
#######get the list of object which is significant after control#############
sigresult<-subset(result, result$padj<0.01)

whichgrange<-as.vector(as.numeric(sigresult$X))

######get the list of interesting object#################
########select the specific complex component to plot################
complex<-c("nejire","ISWI","ACF1","NURF301","CtBP","MLE","MOF","MSL-1","HDAC1","dMi-2","Pho","dSFMBT","WDS","JHDM1","Psc","dRING","PIWI","Rhino","SNR1","brahma","Enhancer-of-zeste","Pc","Su(var)3-9","HP1a","JIL-1","LSD1","HDAC1","Su(var)3-7","H3K9me3") #ISWI,MSL,NURD,PHORC,
component<-subset(result, result$object1 %in% complex)

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




#########################################################################################
##############Use chipseeker to explore the feature of foxo sig peak#####################
########################################################################################
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
peakANNO<-annotatePeak(foxogRange,tssRegion=c(-3000,3000), TxDb=txdb)
foxopeakid<-peakANNO@anno$geneId

peakannoresult<-as.data.frame(as.GRanges(peakANNO))
write.csv(peakannoresult, file=paste(getwd(), "/", "foxosigpeak.csv", sep=''))

######plot the annotation #####################
plotAnnoPie(peakANNO)


upsetplot(peakANNO, vennpie=TRUE)


#############################################################################
######specify some condition on object and plot heatmap######################
############################################################################
######only histone##############
histoneheatmap<-filter(overlapdata, overlapdata$object2 %in% str_extract(overlapdata$object2, "H[0-9][\\s\\S]*"))


#####not histone#############################
nonhistoneheatmap<-setdiff(overlapdata, histoneheatmap)

#rownames(overlapdata)<-overlapdata[,1]



################################################################
######create the heatmap#######################################
##############################################################
###################create the simple heatmap##################
forheatmap<-as.matrix(overlapdata[,-c(1,2)])
rownames(forheatmap)<-object2
colnames(forheatmap)<-foxopeakid

#forheatmaphistone<-as.matrix(histoneheatmap[,-c(1,2)])
#rownames(forheatmaphistone)<-as.vector(as.character(histoneheatmap$object2))


##################ggplot for heatmap##################################
###########reshape the data#######################################

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

##########heatmap by ggplot##################################
ggplot(data = ggheatmap, aes(x = object, y = foxopeak)) +
geom_tile(aes(fill = overlap), colour="steelblue")+
scale_fill_manual(values=cols)

###########heatmap by pheatmap####################
pheatmap(forheatmap, color=c( "gray70","violetred2"), fontsize=9, fontsize_row=6, cellwidth = 7, cellheight = 5, cluster_cols=TRUE)

###################maybe caculate the average effect for each object########
level<-as.character(levels(nonhistone$object1))

averagepval<-list()
number<-list()

for (i in 1:length(level)){
  obj<-level[i]
  objcollection<-subset(nonhistone, nonhistone$object1 %in% obj)
  averagepval[i]<-sum(as.numeric(as.character(objcollection$psig)))/nrow(objcollection)
  number[i]<-nrow(objcollection)
}

averagepval<-as.vector(unlist(averagepval))
number<-as.vector(unlist(number))
averagepadj<-p.adjust(averagepval, method='fdr')

objectalone<-as.data.frame(cbind(level, averagepval, averagepadj, number))

#########try other method(genometric correlation)######################

genomecorelation<-list()
i<-1 

for (i in 1:length(collect)){
try(genomecorelation[[i]]<- GenometriCorrelation(filetableGrange[[i]], foxogRange,  permut.number = 500,keep.distributions = TRUE, showProgressBar = FALSE))
cat(i, "\n")

}


#############basic overlap summary######################
overlapRegions(foxogRange, filetableGrange[[i]], min.pctA= , get.pctA=TRUE, get.bases=TRUE) #Minimum percentage of A being overlapped///how much should be fine?

#############basic graph summary########################
i<-
overlapGraphicalSummary(foxogRange, filetableGrange[[i]], regions.labels=c("foxo", ID[i]), regions.colors=c("black","forestgreen","red"))


#############combine signal and target gRange to a matrix########
mat1 = normalizeToMatrix(foxogRange, filetableGrange[1], value_column = "coverage", 
                         extend = 5000, mean_mode = "w0", w = 50)

