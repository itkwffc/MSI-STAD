rm(list=ls())
library(limma)
immuneFile="GEOCIBERSORT-Results.txt"     #�����ļ�
surFile="GEOTime.txt"                     #�ٴ������ļ�
pFilter=0.05                           #CIBERSORT�����������
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-4.mergeTime")  

#��ȡ���߽���ļ����������ݽ�������
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
data=avereps(data)
#��ȡ��������
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#���ݺϲ���������
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ImmuneCellTime-GEO.txt",sep="\t",row.names=F,quote=F)

rm(list=ls())
immuneFile="chemoexp-GEO.txt"     #�����ļ�
surFile="GEOTime.txt"                     #�ٴ������ļ�
#��ȡ���߽���ļ����������ݽ�������
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
data=as.matrix(t(immune))
data=avereps(data)
#��ȡ��������
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#���ݺϲ���������
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ChemokineTime-GEO.txt",sep="\t",row.names=F,quote=F)
#END