rm(list=ls())
library(limma)

expFile="networkExp.txt"      #表达数据文件
cliFile="TCGA-clinical-data.txt"            #临床数据

setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\02.TCGAmergeTime")    #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)


#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data<-t(data)

#读取生存数据
cli=read.table(cliFile,sep="\t",check.names=F,header=T)     #读取临床文件
cli<-cli[!duplicated(cli$Id),]
rownames(cli)<-cli$Id
cli<-cli[,2:ncol(cli)]
#表数据和生存数据合并输出结果
sameSample=intersect(row.names(data),rownames(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
#输出结果用于后续分析
write.table(out,file="TCGA-expTime.txt",sep="\t",row.names=F,quote=F)
#END