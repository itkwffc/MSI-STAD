rm(list = ls())
library(survival)
library(survminer)
library(ROCR)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\05.modelValidation\\05-2.ROC")  
bioROC=function(inputFile=null,rocFile=null){
	
	rt=read.table(inputFile,header=T,sep="\t")
	col<-rainbow(ncol(rt)-3)
	col[2]<-"#FFB9B9"
	pdf(file=rocFile,width=5,height=5)
	for(i in 4:ncol(rt))
	{
	  temp<-i
	  i<-colnames(rt)[i]
	  data<-rt[,c(i,"fustat")]
	  data<-data[data[,1]!="Unknow",]
	  data[,1]<-as.numeric(data[,1])
	  pred <- prediction(data[,1], data[,2])
	  perf <- performance(pred,"tpr","fpr")
	  aucvalue<-performance(pred,"auc") # shows calculated AUC for model
	  if(temp==4){plot(perf,colorize=F, col=col[temp-3],lwd=2) }else{plot(perf,colorize=F, col=col[temp-3],lwd=2,add=T) }
	  text(x=0.8,y=(0+temp/20-0.2),col =col[temp-3] ,paste(i,"AUC:",round(aucvalue@y.values[[1]],2)))
	}
	lines(c(0,1),c(0,1),col = "black", lty = 4 )
	dev.off()
}

bioROC(inputFile="tcgaRiskCli-STAD.txt",rocFile="tcga.ROC-allCli.pdf")
#END