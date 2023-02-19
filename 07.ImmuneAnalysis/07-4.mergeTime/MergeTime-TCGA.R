rm(list=ls())
library(limma)
immuneFile="TCGACIBERSORT-Results.txt"     #�����ļ�
surFile="TCGATime.txt"                     #�ٴ������ļ�
pFilter=0.05                           #CIBERSORT�����������
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-4.mergeTime")    
#��ȡ���߽���ļ����������ݽ�������
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
#ɾ��������Ʒ
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data=avereps(data)
#��ȡ��������
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#���ݺϲ���������
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ImmuneCellTime-TCGA.txt",sep="\t",row.names=F,quote=F)

rm(list=ls())
immuneFile="chemoexp-TCGA.txt"     #�����ļ�
surFile="TCGATime.txt"                     #�ٴ������ļ�
#��ȡ���߽���ļ����������ݽ�������
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
data=as.matrix(t(immune))
#ɾ��������Ʒ
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data=avereps(data)
#��ȡ��������
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#���ݺϲ���������
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ChemokineTime-TCGA.txt",sep="\t",row.names=F,quote=F)
#END