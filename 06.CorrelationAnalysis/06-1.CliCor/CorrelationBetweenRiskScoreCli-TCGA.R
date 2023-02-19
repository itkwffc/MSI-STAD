rm(list=ls())
options(stringsAsFactors=F)
library(limma)
library(ggpubr)

setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\06.CorrelationAnalysis\\06-1.CliCor")
gene="riskScore"           #感兴趣因素

#读取riskscore文件
rt=read.table("tcgaRisk-STAD.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,c(gene)]
exp=data.frame(exp)
colnames(exp)<-gene
#读取包含ips的临床数据文件
cli=read.table("TCGA-clinical-data-addipsgroup.txt",sep="\t",header=T,check.names=F,row.names=1)
#合并数据riskscore和ips
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)
colnames(rt)[1]<-gene
rt[,1]<-as.numeric(rt[,1])
#临床因素包括ips和riskscore相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c(gene,clinical)]
	colnames(data)=c("RiskScore","clinical")
	data=data[(data[,"clinical"]!="Unknow"),]
	#设置比较组
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
    for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	#绘制boxplot
	boxplot=ggboxplot(data, x="clinical", y="RiskScore", color="clinical",
	          xlab=clinical,
	          ylab="RiskSCore",
	          legend.title=clinical,
	          add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0(clinical,".pdf"),width=5.5,height=5)
	print(boxplot)
	dev.off()
}
#END