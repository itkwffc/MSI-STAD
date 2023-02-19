rm(list = ls())
library("limma")             #引用包
#TCGA
TCGAexpFile="TCGAexp.txt"         #表达输入文件
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-1.CIBERSORT\\TCGA")
#读取输入文件，并对输入文件整理
rt=read.table(TCGAexpFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#去除低表达基因
data=data[rowMeans(data)>0,]
data=rbind(ID=colnames(data),data)
write.table(data,file="TCGAuniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件
#运行CIBERSORT，得到免疫细胞含量结果
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "TCGAuniq.symbol.txt", perm=1000, QN=TRUE)

#GEO
rm(list = ls())
GEOexpFile="GEOexp.txt"
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-1.CIBERSORT\\GEO")
#读取输入文件，并对输入文件整理
rt=read.table(GEOexpFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#去除低表达基因
data=data[rowMeans(data)>0,]
data=rbind(ID=colnames(data),data)
write.table(data,file="GEOuniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件

#运行CIBERSORT，得到免疫细胞含量结果
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "GEOuniq.symbol.txt", perm=1000, QN=TRUE)
#END