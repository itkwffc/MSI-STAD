rm(list=ls())
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
corFilter=0.0            #相关系数过滤标准
pvalueFilter=1       #p值过滤标准
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\06.CorrelationAnalysis\\06-2.HubGeneCor")     #设置工作目录
expFile="expTime-withImmuCheckPoint.txt"            
hubgeneFile<-"HubGene.txt"
hubCernaFile<-"TCGA-expTime.txt"
immuCheckPoint<-c("PDCD1","CD274","CTLA4")
#获得meta数据，rawcount数据和矫正后的表达数据
load("metaDataoffline.RData")
#输出感兴趣的immu checkpoint exp文件
#转录组数据id转化
Ensembl_ID <- rownames(rnaExpr)
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
rnaExpr<-rnaExpr[rownames(rnaExpr)%in%gene_symbol$ENSEMBL,]
rnaExpr<-rnaExpr[gene_symbol$ENSEMBL,]
rownames(rnaExpr)<-gene_symbol$SYMBOL

#筛选checkpoint表达数据
rnaExpr<-rnaExpr[immuCheckPoint,]
rnaExpr<-as.data.frame(t(rnaExpr))

#删除rnaExpr正常样品
group=sapply(strsplit(rownames(rnaExpr),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rnaExpr=rnaExpr[group==0,]
rownames(rnaExpr)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(rnaExpr))

#合并cerna和checkpoint的exp
hubCerna<-read.table(file=hubCernaFile,check.names = F,sep = "\t",header = T)
rownames(hubCerna)<-hubCerna$id
samesample<-intersect(rownames(rnaExpr),rownames(hubCerna))
rnaExpr<-rnaExpr[samesample,]
hubCerna<-hubCerna[samesample,]
hubCerna<-cbind(hubCerna,rnaExpr)

#输出,之后还需要整理
write.table(hubCerna,file = expFile,sep = "\t",quote = F,row.names = F)

#读取表达文件，并对输入文件整理
hubgene<-read.table(file=hubgeneFile,sep="\t",check.names = F,header = F)
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt<-rt[,c(hubgene[,1],immuCheckPoint)]

#相关性检验
outTab=data.frame()
for(i in hubgene[,1]){
  for(j in immuCheckPoint){
    gene1=i
    gene2=j
    if(gene1!=gene2){
      x=as.numeric(rt[,i])
      y=as.numeric(rt[,j])
      corT=cor.test(x,y,method="spearman")#spearman对数据来源分布不做要求属于无参检验，但是检验效力小于pearson
      cor=corT$estimate
      pvalue=corT$p.value
      if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(gene1=gene1,gene2=gene2,cor=cor,corPval=pvalue))
        #绘制相关性曲线
        df1=as.data.frame(cbind(x,y))
        p1=ggplot(df1, aes(x, y)) + 
          xlab(gene1)+ylab(gene2)+
          geom_point()+ geom_smooth(method="lm",formula=y~x) + theme_bw()+
          stat_cor(method = 'spearman', aes(x =x, y =y))
        pdf(file=paste0("cor.",gene1,"_",gene2,".pdf"),width=5,height=5)
        print(p1)
        dev.off()
      }
    }
  }
}
#write.table(file="corResult.xls",outTab,sep="\t",quote=F,row.names=F)
#END