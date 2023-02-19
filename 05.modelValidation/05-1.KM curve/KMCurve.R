rm(list=ls())
library(survival)
library(survminer)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\05.modelValidation\\05-1.KM curve")        

bioSurvival=function(inputFile=null,outFile=null){
		rt=read.table(inputFile,header=T,sep="\t")                  
		#比较高低风险组生存差异，得到显著性p值
		diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
		pValue=1-pchisq(diff$chisq,df=1)
		pValue=signif(pValue,4)
		pValue=format(pValue, scientific = TRUE)
		fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	
		surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=paste0("p=",pValue),
		           pval.size=5,
		           legend.labs=c("High risk", "Low risk"),
		           legend.title="Risk",
		           xlab="Time(years)",
		           break.time.by = 1,
		           risk.table.title="",
		           palette=c("red", "blue"),
		           risk.table=F,
		           risk.table.height=.25)
		pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
		print(surPlot)
		dev.off()
}
bioSurvival(inputFile="tcgaRisk-STAD.txt",outFile="tcga.survival-STAD.pdf")
bioSurvival(inputFile="geoRisk-STAD.txt",outFile="geo.survival-STAD.pdf")
#END