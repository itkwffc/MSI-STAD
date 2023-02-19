rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggExtra)
corFilter=0.0                    #相关系数过滤条件
pvalueFilter=0.05             #相关性检验p值过滤条件
geneRiskFile="ChemokineTime-TCGA.txt"          #gene风险文件
immuneRiskFile="TCGA-expTime.txt"      #immune风险文件
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-6.CoExpressionplot")
#读取gene风险文件
geneRisk=read.table(geneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
geneRisk=geneRisk[,3:(ncol(geneRisk))]
#读取immune风险文件
immuneRisk=read.table(immuneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
immuneRisk=immuneRisk[,11:(ncol(immuneRisk))]

#sift感兴趣的immunecell和chemokine
immunesift<-read.table("HublncRNA.txt",sep = "\t",check.names = F,header = T)
chemokinesift<-read.table("HubChemokine.txt",sep = "\t",check.names = F,header = T)
geneRisk<-geneRisk[,colnames(geneRisk)%in%chemokinesift[,1]]
immuneRisk<-immuneRisk[,colnames(immuneRisk)%in%immunesift[,1],drop=F]
#样品取交集
sameSample=intersect(row.names(geneRisk),row.names(immuneRisk))
geneRisk=geneRisk[sameSample,]
immuneRisk=immuneRisk[sameSample,,drop=F]

#对基因进行循环
outTab=data.frame()
for(gene in colnames(geneRisk)){
    #对免疫细胞进行循环
	for(immune in colnames(immuneRisk)){
		x=as.numeric(geneRisk[,gene])
		y=as.numeric(immuneRisk[,immune])
		corT=cor.test(x,y)
		cor=corT$estimate
		pvalue=corT$p.value
		outTab=rbind(outTab,cbind(gene,immune,cor=cor,corPval=pvalue))
		if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
			#绘制相关性曲线
			df1=as.data.frame(cbind(x,y))
			p1=ggplot(df1, aes(x, y)) + 
					 xlab(gene)+ylab(immune)+
					 geom_point()+ 
					 geom_smooth(method="lm",formula=y~x) + 
					 theme_bw()+
					 stat_cor(aes(x =x, y =y))
			pdfFile=paste0(gene,"_",immune,".pdf")
			pdf(file=pdfFile,width=5,height=5)
			print(p1)
			dev.off()
		}
	}
}
#输出肿瘤类型和相关性表格文件
#write.table(file="corResult.xls",outTab,sep="\t",quote=F,row.names=F)
#END