rm(list = ls())  #һ�����~
library("glmnet")
library("survival")
library("reshape")
library("sva")
library("ggplot2")
coxSigFile="tcga.uniSigExp.txt"       #�����ļ�
geneFile="subceRNAgenes.txt"         #�����ļ�
GEOFile<-"GEOexpTime.txt"            #GEO�����ļ�
Lasso1File<-"model-lamda.pdf"        #����ͼƬ����
Lasso2File<-"deviance.pdf"          #����ͼƬ����
norm<-F
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\04.Lasso")            #���ù���Ŀ¼
#��ȡtrian������
trainMergeTime=read.table(coxSigFile,header=T,sep="\t",row.names=1)              #��ȡ�ļ�
geneRT=read.table(geneFile,header=F,sep="\t",check.names=F)          #��ȡ�����ļ�
trainMergeTime<-trainMergeTime[trainMergeTime$futime>0,]
trainMergeTime$futime=trainMergeTime$futime/365
#��ȡtest������
testMergeTime=read.table(GEOFile,header=T,sep="\t",row.names=1)
testMergeTime<-testMergeTime[testMergeTime$futime>0,]
#��ΪGEO��������ʱ���Ѿ�ͳһ������Ϊ��λ#testMergeTime$futime=testMergeTime$futime/365

#ɸѡ��tcga��geo���ݿⶼ���ڵ�gene
samegene<-intersect(as.vector(geneRT[,1]),colnames(trainMergeTime))
samegene<-intersect(samegene,colnames(testMergeTime))
#samegene<-c("IL1RL1","SPAG16","FAM110B","ANKRD6","ACSS3","CORO2B","TNFAIP8L3")
trainMergeTime=trainMergeTime[,c("futime","fustat",samegene)]
testMergeTime=testMergeTime[,c("futime","fustat",samegene)]
#����ģ��
mergeGeneTimeData<-trainMergeTime
x=as.matrix(mergeGeneTimeData[,c(3:ncol(mergeGeneTimeData))])
y=data.matrix(Surv(mergeGeneTimeData$futime,mergeGeneTimeData$fustat))
fit=glmnet(x, y, family = "cox" ,maxit = 1e+05)#LASSO
#����LASSO profile
pdf(file=Lasso1File,width =8 ,height =6 )
plot(fit, xvar="lambda", label=T)
dev.off()
#����devianceͼ
cvfit=cv.glmnet(x, y, family="cox", maxit = 1e+05,type.measure = "deviance")
pdf(file=Lasso2File,width =8 ,height =6 )
plot(cvfit)
dev.off()

#�����ػ���ϵ��-geneCoef-STAD.txt
coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef-STAD.txt",sep="\t",quote=F,row.names=F)

#���train�����ֵ
trainFinalGeneExp=trainMergeTime[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(trainMergeTime[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="tcgaRisk-STAD.txt",sep="\t",quote=F,row.names=F)

#���test�����ֵ
testFinalGeneExp=testMergeTime[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(testScore),"high","low"))
outTab=cbind(testMergeTime[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="geoRisk-STAD.txt",sep="\t",quote=F,row.names=F)
#END