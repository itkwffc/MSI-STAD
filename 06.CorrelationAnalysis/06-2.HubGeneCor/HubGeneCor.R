rm(list=ls())
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
corFilter=0.0            #���ϵ�����˱�׼
pvalueFilter=1       #pֵ���˱�׼
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\06.CorrelationAnalysis\\06-2.HubGeneCor")     #���ù���Ŀ¼
expFile="expTime-withImmuCheckPoint.txt"            
hubgeneFile<-"HubGene.txt"
hubCernaFile<-"TCGA-expTime.txt"
immuCheckPoint<-c("PDCD1","CD274","CTLA4")
#���meta���ݣ�rawcount���ݺͽ�����ı�������
load("metaDataoffline.RData")
#�������Ȥ��immu checkpoint exp�ļ�
#ת¼������idת��
Ensembl_ID <- rownames(rnaExpr)
gene_symbol <- bitr(Ensembl_ID, fromType="ENSEMBL", toType=c("SYMBOL"), OrgDb="org.Hs.eg.db")
rnaExpr<-rnaExpr[rownames(rnaExpr)%in%gene_symbol$ENSEMBL,]
rnaExpr<-rnaExpr[gene_symbol$ENSEMBL,]
rownames(rnaExpr)<-gene_symbol$SYMBOL

#ɸѡcheckpoint��������
rnaExpr<-rnaExpr[immuCheckPoint,]
rnaExpr<-as.data.frame(t(rnaExpr))

#ɾ��rnaExpr������Ʒ
group=sapply(strsplit(rownames(rnaExpr),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rnaExpr=rnaExpr[group==0,]
rownames(rnaExpr)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(rnaExpr))

#�ϲ�cerna��checkpoint��exp
hubCerna<-read.table(file=hubCernaFile,check.names = F,sep = "\t",header = T)
rownames(hubCerna)<-hubCerna$id
samesample<-intersect(rownames(rnaExpr),rownames(hubCerna))
rnaExpr<-rnaExpr[samesample,]
hubCerna<-hubCerna[samesample,]
hubCerna<-cbind(hubCerna,rnaExpr)

#���,֮����Ҫ����
write.table(hubCerna,file = expFile,sep = "\t",quote = F,row.names = F)

#��ȡ�����ļ������������ļ�����
hubgene<-read.table(file=hubgeneFile,sep="\t",check.names = F,header = F)
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt<-rt[,c(hubgene[,1],immuCheckPoint)]

#����Լ���
outTab=data.frame()
for(i in hubgene[,1]){
  for(j in immuCheckPoint){
    gene1=i
    gene2=j
    if(gene1!=gene2){
      x=as.numeric(rt[,i])
      y=as.numeric(rt[,j])
      corT=cor.test(x,y,method="spearman")#spearman��������Դ�ֲ�����Ҫ�������޲μ��飬���Ǽ���Ч��С��pearson
      cor=corT$estimate
      pvalue=corT$p.value
      if((abs(cor)>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(gene1=gene1,gene2=gene2,cor=cor,corPval=pvalue))
        #�������������
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