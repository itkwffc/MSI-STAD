rm(list=ls())
#���ð�
library(pheatmap)
library(vioplot)
input="GEOCIBERSORT-Results.txt"      #�����ļ�
cliFile<-"GEOCliRiskscore.txt"
outpdf="vioplot-GEO-ImmuneCell.pdf"               #���ͼƬ����
pFilter=0.05                       #CIBERSORT�����������
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\07.ImmuneAnalysis\\07-2.ImmuneCellVioplot")
#��ȡ���߽���ļ����������ݽ�������
immune=read.table(input,sep="\t",header=T,row.names=1,check.names=F)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
#��ȡmsi����risk��Ϣ
risk=read.table(cliFile,sep="\t",header=T,row.names = 1,check.names=F)
#ȡ��ͬ��Ʒ
samesample<-intersect(rownames(data),rownames(risk))
data<-data[samesample,]
risk<-risk[samesample,]
identical(rownames(data),rownames(risk))
#����ͼƬ
group=risk$risk
data1=data[group=="low",]
data2=data[group=="high",]
normalNum=nrow(data1)
tumorNum=nrow(data2)
rt=rbind(data1,data2)
#���С����ͼ
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
#��ÿ������ϸ��ѭ��������vioplot���ͱ�������ɫ��ʾ���߱����ú�ɫ��ʾ
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
#�������ϸ����pֵ�����ļ�
#write.table(outTab,file="diff-geo-risk.result.txt",sep="\t",row.names=F,quote=F)
#END