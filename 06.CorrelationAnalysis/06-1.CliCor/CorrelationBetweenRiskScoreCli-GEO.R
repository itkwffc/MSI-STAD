rm(list=ls())
options(stringsAsFactors=F)
library(limma)
library(ggpubr)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\06.CorrelationAnalysis\\06-1.CliCor")     #修改工作目录
gene="riskScore"           

#读取表达文件，并对输入文件整理
rt=read.table("geoRisk-STAD.txt",sep="\t",header=T,check.names=F) #表达数据文件名称需要根据研究的肿瘤修改
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,c(gene)]
exp=data.frame(exp)
colnames(exp)<-gene
#读取临床数据文件
cli=read.table("GEO-clinical-data-addipsgroup.txt",sep="\t",header=T,check.names=F,row.names=1)
#合并数据
samSample=intersect(row.names(exp),row.names(cli))
exp=exp[samSample,]
cli=cli[samSample,]
rt=cbind(exp,cli)
colnames(rt)[1]<-gene
rt[,1]<-as.numeric(rt[,1])
#临床相关性分析，输出图形结果
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