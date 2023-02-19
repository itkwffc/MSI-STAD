rm(list=ls())
library(limma)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
expFile="TCGAexp.txt"        
riskFile="tcgaRisk-STAD.txt"      
fdrFilter=0.05         
logFCfilter=1             
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\08.EnrichmentAnalysis\\08-2.riskDiffRelated")  
#DEG Analysis
rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]
colnames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data))
data=avereps(t(data))
data=t(data)
risk=read.table(riskFile,header=T,sep="\t",row.names=1,check.names=F)
sameSample=intersect(colnames(data),row.names(risk))
data=data[,sameSample]
risk=risk[sameSample,]
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=data[,row.names(riskLow)]
dataHigh=data[,row.names(riskHigh)]
data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep(1,conNum),rep(2,treatNum))
outTab=data.frame()
for(i in row.names(data)){
	geneName=unlist(strsplit(i,"\\|",))[1]
	geneName=gsub("\\/", "_", geneName)
	rt=rbind(expression=data[i,],Type=Type)
	rt=as.matrix(t(rt))
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=log2(treatGeneMeans)-log2(conGeneMeans)
	pvalue=wilcoxTest$p.value
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,lowMean=conGeneMeans,highMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="riskDiff.txt",sep="\t",row.names=F,quote=F)

#KEGG and GO Analysis
pvalueFilter=0.05        
qvalueFilter=1
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
rt=read.table("riskDiff.txt",sep="\t",check.names=F,header=T)    
#KEGG
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG-RiskDiffRelated.txt",sep="\t",quote=F,row.names = F)
showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}
pdf(file="bubble-RiskDiffRelated-KEGG.pdf",width = 10,height = 7)
aa=dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)
print(aa)
dev.off()

#GO
kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}
pdf(file="bubble-RiskDiffRelated-GO.pdf",width = 10,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
#END