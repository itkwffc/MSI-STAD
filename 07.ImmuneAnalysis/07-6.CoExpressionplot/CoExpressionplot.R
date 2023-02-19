rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggExtra)
corFilter=0.0                    #���ϵ����������
pvalueFilter=0.05             #����Լ���pֵ��������
geneRiskFile="ChemokineTime-TCGA.txt"          #gene�����ļ�
immuneRiskFile="TCGA-expTime.txt"      #immune�����ļ�
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-6.CoExpressionplot")
#��ȡgene�����ļ�
geneRisk=read.table(geneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
geneRisk=geneRisk[,3:(ncol(geneRisk))]
#��ȡimmune�����ļ�
immuneRisk=read.table(immuneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
immuneRisk=immuneRisk[,11:(ncol(immuneRisk))]

#sift����Ȥ��immunecell��chemokine
immunesift<-read.table("HublncRNA.txt",sep = "\t",check.names = F,header = T)
chemokinesift<-read.table("HubChemokine.txt",sep = "\t",check.names = F,header = T)
geneRisk<-geneRisk[,colnames(geneRisk)%in%chemokinesift[,1]]
immuneRisk<-immuneRisk[,colnames(immuneRisk)%in%immunesift[,1],drop=F]
#��Ʒȡ����
sameSample=intersect(row.names(geneRisk),row.names(immuneRisk))
geneRisk=geneRisk[sameSample,]
immuneRisk=immuneRisk[sameSample,,drop=F]

#�Ի������ѭ��
outTab=data.frame()
for(gene in colnames(geneRisk)){
    #������ϸ������ѭ��
	for(immune in colnames(immuneRisk)){
		x=as.numeric(geneRisk[,gene])
		y=as.numeric(immuneRisk[,immune])
		corT=cor.test(x,y)
		cor=corT$estimate
		pvalue=corT$p.value
		outTab=rbind(outTab,cbind(gene,immune,cor=cor,corPval=pvalue))
		if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
			#�������������
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
#����������ͺ�����Ա����ļ�
#write.table(file="corResult.xls",outTab,sep="\t",quote=F,row.names=F)
#END