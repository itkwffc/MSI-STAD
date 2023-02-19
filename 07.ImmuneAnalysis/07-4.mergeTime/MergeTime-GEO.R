rm(list=ls())
library(limma)
immuneFile="GEOCIBERSORT-Results.txt"     #输入文件
surFile="GEOTime.txt"                     #临床数据文件
pFilter=0.05                           #CIBERSORT结果过滤条件
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-4.mergeTime")  

#读取免疫结果文件，并对数据进行整理
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
data=avereps(data)
#读取生存数据
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ImmuneCellTime-GEO.txt",sep="\t",row.names=F,quote=F)

rm(list=ls())
immuneFile="chemoexp-GEO.txt"     #输入文件
surFile="GEOTime.txt"                     #临床数据文件
#读取免疫结果文件，并对数据进行整理
immune=read.table(immuneFile,sep="\t",header=T,row.names=1,check.names=F)
data=as.matrix(t(immune))
data=avereps(data)
#读取生存数据
surTime=read.table(surFile,header=T,sep="\t",check.names=F,row.names=1)
#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
out=cbind(surTime,data)
out=cbind(id=row.names(out),out)
write.table(out,file="ChemokineTime-GEO.txt",sep="\t",row.names=F,quote=F)
#END