rm(list = ls())
library(pheatmap)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\05.modelValidation\\05-5.heatmap")
bioheatmap=function(expFile=null,cliFile=full,outputgile=null){
  #读取数据
  rt=read.table(expFile,header=T,sep="\t",row.names=1,check.names=F)
  rt<-rt[order(rt$riskScore),]
  cli<-read.table(cliFile,header = T,sep = "\t",row.names = 1,check.names = F)
  #取相同数据
  samesample<-intersect(rownames(rt),rownames(cli))
  cli<-cli[samesample,]
  rt<-rt[samesample,]
  identical(row.names(cli),row.names(rt))
  #取表达矩阵
  hmExp<-t(rt[,3:(ncol(rt)-2)])
  #hmExp=log2(hmExp+0.1)
  data<-cbind(cli,rt[,ncol(rt)])
  colnames(data)[ncol(data)]<-"RiskScore"
  data$RiskScore<-factor(data$RiskScore,levels = c("low","high"))
  ann_colors = list(RiskScore = c(low="light blue", high="firebrick"))
  gaps_col=ceiling(nrow(data)/2)
  pdf(file=outputgile,height=8.6,width=10)
  pheatmap(hmExp, 
           annotation=data,annotation_colors = ann_colors, 
           color = colorRampPalette(c("blue", "white", "red"))(50),
           cluster_cols =F,
           show_colnames = F,
           show_rownames = T,
           scale="row",
           gaps_col = c(gaps_col),
           fontsize = 12,
           fontsize_row=10,
           fontsize_col=10)
  dev.off()
}
bioheatmap(expFile="tcgaRisk-STAD.txt" ,cliFile="TCGA-clinical-data-addipsgroup.txt",outputgile="heatmap-STAD-TCGA.pdf")
bioheatmap(expFile="geoRisk-STAD.txt" ,cliFile="GEO-clinical-data-addipsgroup.txt",outputgile="heatmap-STAD-GEO.pdf")
#END