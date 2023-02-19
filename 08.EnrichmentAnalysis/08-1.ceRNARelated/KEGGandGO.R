rm(list=ls())
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05        
qvalueFilter=1           


colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\08.EnrichmentAnalysis\\08-1.ceRNARelated")     
rt=read.table("totalcerna.txt",sep="\t",check.names=F,header=T)    


genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]      

#KEGG
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG-MSIRelated.txt",sep="\t",quote=F,row.names = F)
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}
pdf(file="bubble-MSIRelated-KEGG.pdf",width = 10,height = 7)
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
pdf(file="bubble-MSIRelated-GO.pdf",width = 10,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()
#END