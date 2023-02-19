rm(list=ls())
library(limma)
immuneFile="TCGACIBERSORT-Results.txt"     #输入文件
surFile="TCGATime.txt"                     #临床数据文件
pFilter=0.05                           #CIBERSORT结果过滤条件
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-4.mergeTime")    
#读取免疫结果文件，并对数据进行整理
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data=avereps(data)
#读取生存数据
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ImmuneCellTime-TCGA.txt",sep="\t",row.names=F,quote=F)

rm(list=ls())
immuneFile="chemoexp-TCGA.txt"     #输入文件
surFile="TCGATime.txt"                     #临床数据文件
#读取免疫结果文件，并对数据进行整理
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
data=as.matrix(t(immune))
#删除正常样品
group=sapply(strsplit(row.names(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[group==0,]
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(data))
data=avereps(data)
#读取生存数据
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ChemokineTime-TCGA.txt",sep="\t",row.names=F,quote=F)
#END