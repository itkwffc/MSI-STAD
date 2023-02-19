rm(list=ls())
library(survival)
library(survminer)
pFilter=0.05                                                    
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\03.geneSurvivalAnalysis")     
rt=read.table("TCGA-expTime.txt",header=T,sep="\t",check.names=F,row.names=1)
outTab=data.frame()
sigGenes=c("futime","fustat")
rt<-rt[rt$futime!="Unknow",]
rt$futime<-as.numeric(rt$futime)
#Unicox
for(i in colnames(rt[,11:ncol(rt)])){
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<pFilter){
	    sigGenes=c(sigGenes,i)
			outTab=rbind(outTab,
			             cbind(id=i,
			             HR=coxSummary$conf.int[,"exp(coef)"],
			             HR.95L=coxSummary$conf.int[,"lower .95"],
			             HR.95H=coxSummary$conf.int[,"upper .95"],
			             pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
			             )
	}
}
#输出Unicox结果Table S03
write.table(outTab,file="tcga.uniCox.Result.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
#输出这些有生存意义的gene表达矩阵后续分析
write.table(uniSigExp,file="tcga.uniSigExp.txt",sep="\t",row.names=F,quote=F)
#对有生存意义的lncRNA绘制生存曲线Figure 2
rt$futime=rt$futime/365
rt<-rt[,c("futime","fustat",outTab$id)]
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  if(sd(rt[,i])<0.001){next}
  #KM分析
  group=ifelse(rt[,i]>median(rt[,i]),"high","low")
  diff=survdiff(Surv(futime, fustat) ~group,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outVector=cbind(i,pValue)
  outTab=rbind(outTab,outVector)
  #对p<0.05的基因绘制生存曲线
  if(pValue<pFilter){
    #绘制生存曲线	    
    if(pValue<0.001){pValue="p<0.001"
    }else{pValue=paste0("p=",sprintf("%.03f",pValue))
    }
    fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
    surPlot=ggsurvplot(fit, 
                       data=rt,
                       pval=pValue,
                       conf.int=T,
                       pval.size=5,
                       risk.table=TRUE,
                       legend.labs=c("High", "Low"),
                       legend.title=i,
                       xlab="Time(years)",
                       break.time.by = 1,
                       risk.table.title="",
                       palette=c("red", "blue"),
                       risk.table.height=.25)
    pdf(file=paste0("sur.",i,".pdf"),onefile = FALSE,width = 6.5,height =5.5)
    print(surPlot)
    dev.off()
  }
}
#END