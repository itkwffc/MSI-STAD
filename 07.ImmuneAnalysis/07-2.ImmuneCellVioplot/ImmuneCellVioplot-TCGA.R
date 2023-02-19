rm(list = ls())
#引用包
library(pheatmap)
library(vioplot)
input="TCGACIBERSORT-Results.txt"      #输入文件
cliFile<-"TCGACliRiskscore.txt"
outpdf="vioplot-TCGA-ImmuneCell.pdf"               #输出图片名称
pFilter=0.05                       #CIBERSORT结果过滤条件
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-2.ImmuneCellVioplot")    #设置工作目录

#读取免疫结果文件，并对数据进行整理
immune=read.table(input,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#读取msi或者risk信息
risk=read.table(cliFile,sep="\t",header=T,row.names = 1,check.names=F)
#取共同样品
group=sapply(strsplit(rownames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data<-data[group==0,]
rownames(data)<-substring(rownames(data),1,12)
samesample<-intersect(rownames(data),rownames(risk))
data<-data[samesample,]
risk<-risk[samesample,]
identical(rownames(data),rownames(risk))

group=risk$risk
data1=data[group=="low",]
data2=data[group=="high",]
normalNum=nrow(data1)
tumorNum=nrow(data2)
rt=rbind(data1,data2)
#输出小提琴图
outTab=data.frame()
pdf(outpdf,height=8,width=13)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，低表达用蓝色表示，高表达用红色表示
for(i in 1:ncol(rt)){
  if(sd(rt[1:normalNum,i])==0){
    rt[1,i]=0.001
  }
  if(sd(rt[(normalNum+1):(normalNum+tumorNum),i])==0){
    rt[(normalNum+1),i]=0.001
  }
  normalData=rt[1:normalNum,i]
  tumorData=rt[(normalNum+1):(normalNum+tumorNum),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=wilcoxTest$p.value
  if(p<pFilter){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("LowRisk", "HighRisk"),
       lwd=5,bty="n",cex=1.5,
       col=c("blue","red"))
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()
#输出免疫细胞和p值表格文件
#write.table(outTab,file="diff-tcga-risk.result.txt",sep="\t",row.names=F,quote=F)
#END