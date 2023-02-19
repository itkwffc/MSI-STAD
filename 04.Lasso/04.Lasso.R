rm(list = ls())  #一键清空~
library("glmnet")
library("survival")
library("reshape")
library("sva")
library("ggplot2")
coxSigFile="tcga.uniSigExp.txt"       #输入文件
geneFile="subceRNAgenes.txt"         #基因文件
GEOFile<-"GEOexpTime.txt"            #GEO输入文件
Lasso1File<-"model-lamda.pdf"        #绘制图片的名
Lasso2File<-"deviance.pdf"          #绘制图片的名
norm<-F
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\04.Lasso")            #设置工作目录
#读取trian组数据
trainMergeTime=read.table(coxSigFile,header=T,sep="\t",row.names=1)              #读取文件
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)          #读取基因文件
trainMergeTime<-trainMergeTime[trainMergeTime$futime>0,]
trainMergeTime$futime=trainMergeTime$futime/365
#读取test组数据
testMergeTime=read.table(GEOFile,header=T,sep="\t",row.names=1)
testMergeTime<-testMergeTime[testMergeTime$futime>0,]
#因为GEO数据生存时间已经统一成以年为单位#testMergeTime$futime=testMergeTime$futime/365

#筛选在tcga和geo数据库都存在的gene
samegene<-intersect(as.vector(geneRT[,1]),colnames(trainMergeTime))
samegene<-intersect(samegene,colnames(testMergeTime))
#samegene<-c("IL1RL1","SPAG16","FAM110B","ANKRD6","ACSS3","CORO2B","TNFAIP8L3")
trainMergeTime=trainMergeTime[,c("futime","fustat",samegene)]
testMergeTime=testMergeTime[,c("futime","fustat",samegene)]
#构建模型
mergeGeneTimeData<-trainMergeTime
x=as.matrix(mergeGeneTimeData[,c(3:ncol(mergeGeneTimeData))])
y=data.matrix(Surv(mergeGeneTimeData$futime,mergeGeneTimeData$fustat))
fit=glmnet(x, y, family = "cox" ,maxit = 1e+05)#LASSO
#绘制LASSO profile
pdf(file=Lasso1File,width =8 ,height =6 )
plot(fit, xvar="lambda", label=T)
dev.off()
#绘制deviance图
cvfit=cv.glmnet(x, y, family="cox", maxit = 1e+05,type.measure = "deviance")
pdf(file=Lasso2File,width =8 ,height =6 )
plot(cvfit)
dev.off()

#输出相关基因系数-geneCoef-STAD.txt
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef-STAD.txt",sep="\t",quote=F,row.names=F)

#输出train组风险值
trainFinalGeneExp=trainMergeTime[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(trainMergeTime[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="tcgaRisk-STAD.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
testFinalGeneExp=testMergeTime[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(testMergeTime[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="geoRisk-STAD.txt",sep="\t",quote=F,row.names=F)
#END