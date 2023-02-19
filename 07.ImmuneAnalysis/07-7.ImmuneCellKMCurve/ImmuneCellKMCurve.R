rm(list = ls())
library(survival)
library(survminer)
kmPfilter=0.05           #KM方法显著性过滤标准
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-7.ImmuneCellKMCurve")
rt=read.table("ImmuneCellTime-GEO.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取输入文件
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.001){next}
	#KM分析
	group=ifelse(rt[,i]>median(rt[,i]),"high","low")
	diff=survdiff(Surv(futime, fustat) ~group,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)

	#对p<0.05的免疫细胞绘制生存曲线
	if(pValue<kmPfilter){
		outVector=cbind(i,pValue)
		outTab=rbind(outTab,outVector)
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
				           risk.table.height=.35)
		pdf(file=paste0("sur.",i,".pdf"),onefile = FALSE,width = 8,height =7)
		print(surPlot)
		dev.off()
	}
}
#输出基因和p值表格文件
#write.table(outTab,file="immuneSur.result-GEO.xls",sep="\t",row.names=F,quote=F)
#END