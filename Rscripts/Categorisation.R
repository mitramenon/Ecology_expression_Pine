###############################################
####Categorising transcripts using Fst-qst. 
#mitra menon: mitra.menon28@gmail.com
#last modified: 10/21/2023
##BS= Cold garden
##WP= Warm garden



#######################
library(data.table)

#load simulated data of qst obtained the transcriptome using the script 

Qst_BSsims<-fread("QstBS_combined.txt",sep="\t",header=T,data.table=F)
Qst_WPsims<-fread("QstWP_sim.txt",sep="\t",header=T,data.table=F)

#load Fst for the same set of maternal trees and decide threshold for neutral transcripts
#these are obtained with ddRAD-seq using varcomp
Fst<-read.csv("Fst_rnaSeqSamples.csv")
cutoff=quantile(Fst$Fst,probs = 0.95)

#point estimates for qst obtained from `MatVals_varPart.R`
Qst<-read.table("VarComps_Qst.txt",header=T,sep="\t")
Qst$diff<-abs(Qst$Qst_BS-Qst$Qst_WP)

Qst_BS_complete<-apply(Qst_BSsims,2, function(X) return(sum(complete.cases(X))))
Qst_WP_complete<-apply(Qst_WPsims,2, function(X) return(sum(complete.cases(X))))



#Define categories
QstBS_adap<-apply(Qst_BSsims,2,function(X) return(ifelse(quantile(X,probs =0.025,na.rm=T) >= cutoff , "A" , "N" )))
QstBS_adap<-QstBS_adap[QstBS_adap=="A"] 

QstWP_adap<-apply(Qst_WPsims,2,function(X) return(ifelse(quantile(X,probs =0.025,na.rm=T) >=cutoff, "A" , "N" )))
QstWP_adap<-QstWP_adap[QstWP_adap=="A"] 


#(1)Conditionally adaptive in Cold garden (BS)
outBS_sigEST<-Qst[Qst$Trans%in%names(QstBS_adap), ]
outBS_sigEST<-outBS_sigEST[!(outBS_sigEST$Trans%in%names(QstWP_adap)), ] 
outBS_sigEST<-outBS_sigEST[outBS_sigEST$diff>0, ] 
outBS_sigEST<-outBS_sigEST[complete.cases(outBS_sigEST$diff), ] 

mean(outBS_sigEST$Qst_BS)
median(outBS_sigEST$Qst_BS)
sd(outBS_sigEST$Qst_BS)


#(2)Conditionally adaptive in warm garden (WP)
outWP_sigEST<-Qst[Qst$Trans%in%names(QstWP_adap), ]
outWP_sigEST<-outWP_sigEST[!(outWP_sigEST$Trans%in%names(QstBS_adap)), ] 
outWP_sigEST<-outWP_sigEST[outWP_sigEST$diff>0, ] 
outWP_sigEST<-outWP_sigEST[complete.cases(outWP_sigEST$diff), ]

mean(outWP_sigEST$Qst_WP)
median(outWP_sigEST$Qst_WP)
sd(outWP_sigEST$Qst_WP)

#(3)Adaptive in both gardens but not plastic 
adap<-Qst[Qst$Trans%in%names(QstWP_adap), ]
adap<-adap[adap$Trans%in%names(QstBS_adap), ] 
nonPlastic<-adap[adap$diff==0, ]
nonPlastic<-nonPlastic[complete.cases(nonPlastic$diff), ] 

#(4)Adaptive plasticity 
result<-NULL
for (t in 1:nrow(adap)){
trans<-adap$Trans[t]
CI_bs<-Qst_BSsims[ ,trans]
CI_wp<-Qst_WPsims[ ,trans]
result[[t]]<-t.test(CI_bs,CI_wp,paired=TRUE)$p.value
}

adap<-cbind(adap,t_test=unlist(result))
adap<-adap[order(adap$t_test), ]
AdaptPlast<-adap[adap$t_test<=0.05, ] 

AdaptPlast<-AdaptPlast[complete.cases(AdaptPlast$diff), ] 

#Save the categories
save(outWP_sigEST, outBS_sigEST, AdaptPlast, file = "categories.RData")