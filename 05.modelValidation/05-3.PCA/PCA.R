rm(list=ls())
library(ggplot2)
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\05.modelValidation\\05-3.PCA")   
bioPCA=function(inputFile=null, pcaFile=null){
	rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)
	data=rt[c(3:(ncol(rt)-2))]
	risk=rt[,"risk"]
	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
	#ªÊ÷∆PCAÕº
	pdf(file=pcaFile, height=4.5, width=5.5)   
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
}
bioPCA(inputFile="tcgaRisk-STAD.txt", pcaFile="tcga.PCA.pdf")
bioPCA(inputFile="geoRisk-STAD.txt", pcaFile="geo.PCA.pdf")
#END