rm(list = ls())
library(pheatmap)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\05.modelValidation\\05-4.riskPlot")         
bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null){
		rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)   #读取输入文件
		rt=rt[order(rt$riskScore),]    
	
		riskClass=rt[,"risk"]
		lowLength=length(riskClass[riskClass=="low"])
		highLength=length(riskClass[riskClass=="high"])
		lowMax=max(rt$riskScore[riskClass=="low"])
		line=rt[,"riskScore"]
		line[line>10]=10
		pdf(file=riskScoreFile,width = 8,height = 6)
		plot(line, type="p", pch=20,
		     xlab="Patients (increasing risk socre)", ylab="Risk score",
		     col=c(rep("blue",lowLength),rep("red",highLength)) )
		abline(h=lowMax,v=lowLength,lty=2)
		legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
		dev.off()
		
	
		color=as.vector(rt$fustat)
		color[color==1]="red"
		color[color==0]="blue"
		pdf(file=survStatFile,width = 8,height = 6)
		plot(rt$futime, pch=19,
		     xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		     col=color)
		legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
		abline(v=lowLength,lty=2)
		dev.off()
}
bioRiskPlot(inputFile="geoRisk-STAD.txt",riskScoreFile="geo.riskScore.pdf",survStatFile="geo.survStat.pdf")
bioRiskPlot(inputFile="tcgaRisk-STAD.txt",riskScoreFile="tcga.riskScore.pdf",survStatFile="tcga.survStat.pdf")
#END