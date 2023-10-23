##################Go annotations and enrichment analysis##################
####################################################################
#Mitra Menon: mitra.menon28@gmail.com
#last modifid: 10/21/2023

library(data.table)
library(seqinr)
##############################
#Function to prep table for input to agriGO

InAgriGO<-function(entapOut){
  
  ##entapOut is a list containing the name of the transcript as the header and the GO terms under it as contents
  
  agriGo<-NULL
  goTerm<-NULL
  for (i in 1:length(entapOut)){
    agriGo[[i]]<-rep(names(entapOut[i]),length(entapOut[[i]]))
    goTerm[[i]]<-entapOut[[i]]
  }
  
  InAgri<-data.frame(Trans=unlist(agriGo),GO=unlist(goTerm))
  return(InAgri)
}

#Load annotations obtained from EnTAP
annot<-fread("final_annotations_no_contam_lvl4.tsv",sep="\t",header=T,data.table=F)


Bio<-strsplit(annot[ ,"GO Biological"],",")
names(Bio)<-annot[ ,1]
Bio<-lapply(Bio,function(L) sapply(strsplit(L,"-"),"[",1))
Bio<-lapply(Bio,function(L) grep("^GO", L, value = TRUE)) #the grep will keep only stuff starting with GO

Cell<-strsplit(annot[ ,"GO Cellular"],",")
names(Cell)<-annot[ ,1]
Cell<-lapply(Cell,function(L) sapply(strsplit(L,"-"),"[",1))
Cell<-lapply(Cell,function(L) grep("^GO", L, value = TRUE))


Mol<-strsplit(annot[ ,"GO Molecular"],",")
names(Mol)<-annot[ ,1]
Mol<-lapply(Mol,function(L) sapply(strsplit(L,"-"),"[",1))
Mol<-lapply(Mol,function(L) grep("^GO", L, value = TRUE))

funcs<-NULL
for (i in 1:length(Bio)){
  funcs[[i]]<-c(Bio[[i]],Cell[[i]],Mol[[i]])
}

funcs<-lapply(funcs,function(X) return(unique(X)))
names(funcs)<-annot[ ,1]

#annot<-fread("~/Desktop/Chapter3_macDavis/Go_enrich/final_annotations_no_contam_lvl4.tsv",sep="\t",header=T,data.table=F)
entap<-read.fasta("final_annotated.faa")
funcs_entap<-funcs[names(funcs)%in%names(entap)]


#################################
#Subset of each category on which we would perform GO enrichment
##############################################
Qst<-read.table("VarComps_Qst.txt",header=T,sep="\t")
load("categories.RData")

BS_adap_annot<-funcs_entap[names(funcs_entap)%in%outBS_sigEST$Trans]
WP_adap_annot<-funcs_entap[names(funcs_entap)%in%outWP_sigEST$Trans]


AdapBOTH<-funcs[names(funcs)%in%adap$Trans]
AdapPLA<-funcs_entap[names(funcs_entap)%in%AdaptPlast$Trans]




Ref<-InAgriGO(funcs_entap)
Ref$GO<-as.character(Ref$GO)
Ref$filter<-sapply(strsplit(Ref$GO,"GO:"),"[",2)
Ref<-Ref[complete.cases(Ref$filter), ]
#Ref<-Ref[Ref$Trans%in%Qst$Trans, ]

Bs_annot<-InAgriGO(BS_adap_annot)
Bs_annot$GO<-as.character(Bs_annot$GO)
Bs_annot$filter<-sapply(strsplit(Bs_annot$GO,"GO:"),"[",2)
Bs_annot<-Bs_annot[complete.cases(Bs_annot$filter), ]

WP_annot<-InAgriGO(WP_adap_annot)
WP_annot$GO<-as.character(WP_annot$GO)
WP_annot$filter<-sapply(strsplit(WP_annot$GO,"GO:"),"[",2)
WP_annot<-WP_annot[complete.cases(WP_annot$filter), ]

Adap_annot<-InAgriGO(AdapPLA)
Adap_annot$GO<-as.character(Adap_annot$GO)
Adap_annot$filter<-sapply(strsplit(Adap_annot$GO,"GO:"),"[",2)
Adap_annot<-Adap_annot[complete.cases(Adap_annot$filter), ]

##############################################################################################
###Using GofuncR##########multiple testing accounted via family-wise error rates
################################################################################################

########Using all annotations########
library(GOfuncR)

Bs_annot$Cat<-rep(1,nrow(Bs_annot))
Back<-Ref[!(Ref$Trans%in%Bs_annot$Trans), ]
Back$Cat<-rep(0,nrow(Back))

BS_test<-rbind(Bs_annot[ ,c(1,4)],Back[ ,c(1,4)])#0 for background and 1 for candidate
#gene_lt<-entap[names(entap)%in%BS_test$Trans]#length of each transcript #maybe not use this if we get issues with having to use TaxDB
CAnnot_entap<-Ref[ ,c(1,2)] #df with geneID in col1 and GO annots in col2, if there are multiple could try using get_parent_node
#CAnnot_entap<-CAnnot_entap[match(BS_test$Trans,CAnnot_entap$Trans), ] #acting odd due to duplicate trans IDs

GoEnrich_BP_BS<-go_enrich(genes=BS_test, test = 'hyper', n_randsets = 1000,
                          annotations = CAnnot_entap,organismDb = NULL,silent=TRUE)
test<-GoEnrich_BP_BS$results #enriched terms
testOV<-test[test$FWER_overrep<=0.05, ]
testUN<-test[test$FWER_underrep<=0.05, ]
write.table(rbind(testOV,testUN),file="BS_CondAdap_Qstquantile.txt",sep="\t",row.names = F,quote=F)


WP_annot$Cat<-rep(1,nrow(WP_annot))
Back<-Ref[!(Ref$Trans%in%WP_annot$Trans), ]
Back$Cat<-rep(0,nrow(Back))
WP_test<-rbind(WP_annot[ ,c(1,4)],Back[ ,c(1,4)])
CAnnot_entapWP<-Ref[ ,c(1,2)]
#CAnnot_entapWP<-CAnnot_entapWP[match(WP_test$Trans,CAnnot_entapWP$Trans), ]
GoEnrich_BP_WP<-go_enrich(genes=WP_test, test = 'hyper', n_randsets = 1000,annotations = CAnnot_entapWP,organismDb = NULL,silent = T)
testWP<-GoEnrich_BP_WP$results #enriched terms
testWP_ov<-testWP[testWP$FWER_overrep<0.05, ]
testWP_und<-testWP[testWP$FWER_underrep<0.05, ]

write.table(rbind(testWP_ov,testWP_und),file="WP_CondAdap_Qstquantile.txt",sep="\t",row.names = F,quote=F)


Adap_annot$Cat<-rep(1,nrow(Adap_annot))
Back<-Ref[!(Ref$Trans%in%Adap_annot$Trans), ]
Back$Cat<-rep(0,nrow(Back))
Adap_test<-rbind(Adap_annot[ ,c(1,4)],Back[ ,c(1,4)])
CAnnot_entapAdap<-Ref[ ,c(1,2)]
#CAnnot_entapAdap<-CAnnot_entapAdap[match(Adap_test$Trans,CAnnot_entapAdap$Trans), ]
GoEnrich_BP_Adap<-go_enrich(genes=Adap_test, test = 'hyper', n_randsets = 1000,annotations = CAnnot_entapAdap,organismDb = NULL,silent = T)
testAdap<-GoEnrich_BP_Adap$results #enriched terms
testAdap_ov<-testAdap[testAdap$FWER_overrep<=0.05, ]
testAdap_und<-testAdap[testAdap$FWER_underrep<=0.05, ]
write.table(rbind(testAdap_ov,testAdap_und),file="AdapPlast_Qstquantile.txt",sep="\t",row.names = F,quote=F)

