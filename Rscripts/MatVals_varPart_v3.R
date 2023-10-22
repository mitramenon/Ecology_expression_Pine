##Estimation of maternal values using variancePart package with normalisation and weights######
##############################################################################################
###2-18-2020#####
###Mitra Menon


library(edgeR)
library(variancePartition)
library(BiocParallel)
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)
library(lme4)


##Data loading and prep###############

expISO<-read.table("../Out_trinity/RSEM.isoform.counts.matrix", header=T, row.names=1, com='', check.names=F)
colnames(expISO)[179:180]<-c("6150_rsem","13740_rsem")
IDs<-colnames(expISO)
IDs<-data.frame(Garden_Tag=sapply(strsplit(IDs,"_"),'[',1))

samplesID<-read.table("../Out_trinity/samplesID_RNAseq.txt",header=T,sep="\t")
samplesID<-samplesID[match(IDs$Garden_Tag,samplesID$Garden_Tag), ] #REORDER TO MATCH THE RNASEQ DATASET

##Prep specific for RNAseq#################################

rnaseqMatrix = as.matrix(expISO)
rnaseqMatrix = round(rnaseqMatrix)
exp_study = DGEList(counts=rnaseqMatrix, group=samplesID$group)

keep <- filterByExpr(exp_study)
y <- exp_study[keep, , keep.lib.sizes=FALSE]
y_ISO<-calcNormFactors(y)
eff.lib.size<-y_ISO$samples$lib.size * y_ISO$samples$norm.factors
save(y_ISO,"exp_ISO.rds")

###Prep for running variancePart#######################
samplesID$Garden<-as.factor(samplesID$Garden)
samplesID$SY<-as.factor(samplesID$SY)
samplesID$Pop<-as.factor(samplesID$Pop)
samplesID$TreeID<-as.factor(samplesID$TreeID)
samplesID$group2<-paste(samplesID$Pop,samplesID$Garden,sep = ":")
rownames(samplesID)<-samplesID$Garden_Tag
colnames(y_ISO$counts)<-samplesID$Garden_Tag

samplesID$SY<-as.numeric(samplesID$SY)
samplesID$Note<-as.numeric(samplesID$Note)


######################################################################
###Run varpart and extract componenets for each gene#####
formFull <- ~ Note + SY+  (0+Garden|Pop/TreeID)

# estimate weights using linear mixed model of dream
vobjDream= voomWithDreamWeights(y_ISO, formFull , samplesID)
results<- fitVarPartModel(vobjDream, formFull,samplesID,showWarnings = F)


##Pulling out various variance components

geneResults<-lapply(results,function(df) return(ranef(df)))
results_sum<-lapply(results,function(df) return(summary(df)))

GI<-lapply(results_sum,function(t) t[["coefficients"]][1])

FamVarBS<-lapply(results_sum,function(t) t[["varcor"]]$`TreeID:Pop`[1,"GardenBS"])
FamVarWP<-lapply(results_sum,function(t) t[["varcor"]]$`TreeID:Pop`[2,"GardenWP"])

PopVarBS<-lapply(results_sum,function(t) t[["varcor"]]$Pop[1,"GardenBS"])
PopVarWP<-lapply(results_sum,function(t) t[["varcor"]]$Pop[2,"GardenWP"])

VarComps<-data.frame(Trans=names(geneResults),FamVar_WP= unlist(FamVarWP),FamVar_BS=unlist(FamVarBS),
                     PopVar_WP=unlist(PopVarWP),PopVar_BS=unlist(PopVarBS))
VarComps$FamSlope<-VarComps$FamVar_WP-VarComps$FamVar_BS
VarComps$PopSlope<-VarComps$PopVar_WP-VarComps$PopVar_BS

write.table(VarComps,file="VarComps_Qst.txt",sep="\t",row.names = F,quote=F)


##Obtaining maternal values#################3
  
MatVals<-vector("list",length(geneResults))
names(MatVals)<-names(geneResults)
PopVals<-vector("list",length(geneResults))
names(PopVals)<-names(geneResults)

for (i in 1:length(MatVals)){
  
  
  ID<-names(geneResults)[i]
  Trans<-geneResults[[i]]
  
  FamE_BS<-Trans$`TreeID:Pop`$GardenBS+ GI[[i]]
  FamE_WP<-Trans$`TreeID:Pop`$GardenWP+ GI[[i]]
  PopE_BS<-Trans$Pop$GardenBS
  PopE_WP<-Trans$Pop$GardenWP
  PopE_BS<-rep(PopE_BS,each=3)
  PopE_WP<-rep(PopE_WP,each=3)
  
  FamVal_BS<-FamE_BS+PopE_BS
  FamVal_WP<-FamE_WP+PopE_WP
  
  MatVals[[i]]<-data.frame(locs=rownames(Trans$`TreeID:Pop`),BS=FamVal_BS,WP=FamVal_WP)
  PopVals[[i]]<-data.frame(locs=rownames(Trans$`TreeID:Pop`),BS= PopE_BS+ GI[[i]], WP=PopE_WP+ GI[[i]])
  
  #outF<-paste0("MatVal_",ID,".txt")
  #write.table(MatVals[[i]],file=outF,sep="\t",row.names = F,quote=F)
  
}

FamID<-sapply(strsplit(as.character(MatVals[[1]]$locs),":"),"[",1)
PopID<-sapply(strsplit(as.character(MatVals[[1]]$locs),":"),"[",2)
  
FamValBS<-lapply(MatVals,function(df) return(df[ ,2]))
FamValBS<-do.call(cbind,FamValBS)
FamValBS<-cbind(Fam=FamID,Pop=PopID,FamValBS)

FamValWP<-lapply(MatVals,function(df) return(df[ ,3]))
FamValWP<-do.call(cbind,FamValWP)
FamValWP<-cbind(Fam=FamID,Pop=PopID,FamValWP)

PopValBS<-lapply(PopVals,function(df) return(df[ ,2]))
PopValBS<-do.call(cbind,PopValBS)
PopValBS<-cbind(Pop=PopID,PopValBS)

PopValWP<-lapply(PopVals,function(df) return(df[ ,3]))
PopValWP<-do.call(cbind,PopValWP)
PopValWP<-cbind(Pop=PopID,PopValWP)

