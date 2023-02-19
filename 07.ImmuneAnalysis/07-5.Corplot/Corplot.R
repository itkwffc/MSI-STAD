setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-5.Corplot")     
#TCGA
rm(list=ls())
library(corrplot)                    #���ð�
geneRiskFile="ChemokineTime-TCGA.txt"          #gene�����ļ�
immuneRiskFile="ImmuneCellTime-TCGA.txt"      #immune�����ļ�
#��ȡgene�����ļ�
geneRisk=read.table(geneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
geneRisk=geneRisk[,3:(ncol(geneRisk))]
#��ȡimmune�����ļ�
immuneRisk=read.table(immuneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
immuneRisk=immuneRisk[,3:(ncol(immuneRisk))]
#ɸѡ����Ȥ��immunecell��chemokine
immunesift<-read.table("HubimmuneCell.txt",sep = "\t",check.names = F,header = T)
chemokinesift<-read.table("HubChemokine.txt",sep = "\t",check.names = F,header = T)
geneRisk<-geneRisk[,colnames(geneRisk)%in%chemokinesift[,1]]
immuneRisk<-immuneRisk[,colnames(immuneRisk)%in%immunesift[,1]]
#��Ʒȡ����
sameSample=intersect(row.names(geneRisk),row.names(immuneRisk))
geneRisk=geneRisk[sameSample,]
immuneRisk=immuneRisk[sameSample,]
data=cbind(geneRisk,immuneRisk)
#���������ͼ��
pdf("corrplot-TCGA.pdf",height=14.5,width=14.5)              #����ͼƬ���ļ�����
par(oma=c(0.5,0.5,0.5,1.2))
data=data[,colMeans(data)>0]
M=cor(data)
corrplot(M,
         order="original",
         method = "color",
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

#GEO
rm(list=ls())
library(corrplot)                    #���ð�
geneRiskFile="ChemokineTime-GEO.txt"          #gene�����ļ�
immuneRiskFile="ImmuneCellTime-GEO.txt"      #immune�����ļ�
#��ȡgene�����ļ�
geneRisk=read.table(geneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
geneRisk=geneRisk[,3:(ncol(geneRisk))]
#��ȡimmune�����ļ�
immuneRisk=read.table(immuneRiskFile,header=T,sep="\t",check.names=F,row.names=1)
immuneRisk=immuneRisk[,3:(ncol(immuneRisk))]
#ɸѡ����Ȥ��immunecell��chemokine
immunesift<-read.table("HubimmuneCell.txt",sep = "\t",check.names = F,header = T)
chemokinesift<-read.table("HubChemokine.txt",sep = "\t",check.names = F,header = T)
geneRisk<-geneRisk[,colnames(geneRisk)%in%chemokinesift[,1]]
immuneRisk<-immuneRisk[,colnames(immuneRisk)%in%immunesift[,1]]
#��Ʒȡ����
sameSample=intersect(row.names(geneRisk),row.names(immuneRisk))
geneRisk=geneRisk[sameSample,]
immuneRisk=immuneRisk[sameSample,]
data=cbind(geneRisk,immuneRisk)
#���������ͼ��
pdf("corrplot-GEO.pdf",height=14.5,width=14.5)              #����ͼƬ���ļ�����
par(oma=c(0.5,0.5,0.5,1.2))
data=data[,colMeans(data)>0]
M=cor(data)
corrplot(M,
         order="original",
         method = "color",
         addCoef.col = "black",
         diag = TRUE,
         tl.col="black",
         col=colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
#END