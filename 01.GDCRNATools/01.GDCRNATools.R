rm(list = ls())#清空工作空间
library(GDCRNATools)#引用包
library(ggplot2)#引用包

adjpFilter=0.05          #差异分析FDR阈值
logFCfilter.mrna=1          #差异分析logFC阈值mrna
FCfilter.mrna=2^logFCfilter.mrna
logFCfilter.mirna=0.5            #差异分析logFC阈值mirna
FCfilter.mirna=2^logFCfilter.mirna
logFCfilter.lnrna=1            #差异分析logFC阈值lnrna
FCfilter.lnrna=2^logFCfilter.lnrna
hyperPfilter=0.01        #超几何分布检验pvalue阈值
corPfilter=0.01          #相关性检验pvalue阈值
offline<-T              #获取TCGA-STAD数据数据的方式T离线F在线下载
MSIfile<-'MSI-prediction-STAD.txt' #读取MSI数据
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\01.GDCRNATools")        #设置工作目录

cancer<-c("STAD")#肿瘤类型，根据研究的肿瘤修改
project <- paste("TCGA",cancer,sep = "-")       
pros<-paste("TCGA",paste(cancer,collapse = ""),sep = "-")
rnadir <- paste(pros, 'RNAseq', sep='/')
mirdir <- paste(pros, 'miRNAs', sep='/')

###获取TCGA-STAD数据可以选择在线或者离线两种模式
if(offline==T){
    load("metaDataoffline.RData")
    adjpFilter=0.05          #差异分析FDR阈值
    logFCfilter.mrna=1          #差异分析logFC阈值mrna
    FCfilter.mrna=2^logFCfilter.mrna
    logFCfilter.mirna=0.5            #差异分析logFC阈值mirna
    FCfilter.mirna=2^logFCfilter.mirna
    logFCfilter.lnrna=1            #差异分析logFC阈值lnrna
    FCfilter.lnrna=2^logFCfilter.lnrna
    hyperPfilter=0.01        #超几何分布检验pvalue阈值
    corPfilter=0.01          #相关性检验pvalue阈值
}else{
    for (i in 1:length(project)){
        if (i==1){
            metaMatrix.RNA<-temp.metaMatrix.RNA<-gdcParseMetadata(project.id = project[i],
                                                                  data.type  = 'RNAseq', 
                                                                  write.meta = F)
            
        }else{
            temp.metaMatrix.RNA<-gdcParseMetadata(project.id = project[i],
                                                  data.type  = 'RNAseq', 
                                                  write.meta = F)
            metaMatrix.RNA<-rbind(metaMatrix.RNA,temp.metaMatrix.RNA)
            
        }
        
    }
    
    #转录组metadata过滤
    metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
    metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
    meta.RNA.filter = dir(path=rnadir)#读取rnaseq目录中的文件名，将在线meta文件精简
    metaMatrix.RNA<-metaMatrix.RNA[metaMatrix.RNA$file_id%in%meta.RNA.filter,]
    
    #获取miRNA metadata
    for (i in 1:length(project)){
        if (i==1){
            metaMatrix.MIR<-temp.metaMatrix.MIR<-gdcParseMetadata(project.id = project[i],
                                                                  data.type  = 'miRNAs', 
                                                                  write.meta = F)
            
        }else{
            temp.metaMatrix.MIR<-gdcParseMetadata(project.id = project[i],
                                                  data.type  = 'miRNAs', 
                                                  write.meta = F)
            metaMatrix.MIR<-rbind(metaMatrix.MIR,temp.metaMatrix.MIR)
            
        }
        
    }
    
    #miRNA metadata过滤
    metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
    metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)
    meta.MIR.filter = dir(path=mirdir)
    metaMatrix.MIR<-metaMatrix.MIR[metaMatrix.MIR$file_id%in%meta.MIR.filter,]
    
    
    #转录组数据合并
    rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                             path      = rnadir, 
                             organized = FALSE,   ## if target data are in folders
                             data.type = 'RNAseq')
    
    #miRNA数据合并
    mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                             path      = mirdir,
                             organized = FALSE,   ## if target data are in folders
                             data.type = 'miRNAs')
    
    #转录组数据矫正
    rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
    
    #miRNA数据矫正
    mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
}

###将MSI信息插入RNA和MRI的meta信息
MSIinfo<-read.table(file=MSIfile,check.names = F,header = T)
MSIinfo<-MSIinfo[!duplicated(MSIinfo$patient),]
metaMatrix.RNAall<-merge(metaMatrix.RNA,MSIinfo,by="patient",all=T)
metaMatrix.MIRall<-merge(metaMatrix.MIR,MSIinfo,by="patient",all=T)
metaMatrix.RNA<-metaMatrix.RNAall[metaMatrix.RNAall$patient%in%metaMatrix.RNA$patient,]
metaMatrix.MIR<-metaMatrix.MIRall[metaMatrix.MIRall$patient%in%metaMatrix.MIR$patient,]

###过滤na以及非肿瘤样品
metaMatrix.RNA.DEG <- metaMatrix.RNA[!is.na(metaMatrix.RNA$MSI),]#na.omit是删除所有有na的行，仅删除msi有na的
metaMatrix.MIR.DEG<- metaMatrix.MIR[!is.na(metaMatrix.MIR$MSI),]
metaMatrix.RNA.DEG <- metaMatrix.RNA.DEG [!duplicated(metaMatrix.RNA.DEG$patient),]
metaMatrix.MIR.DEG<- metaMatrix.MIR.DEG[!duplicated(metaMatrix.MIR.DEG$patient),]
metaMatrix.RNA.DEG$MSI<-factor(metaMatrix.RNA.DEG$MSI)
metaMatrix.MIR.DEG$MSI<-factor(metaMatrix.MIR.DEG$MSI)
group=sapply(strsplit(colnames(rnaCounts),"\\-"),"[",4)#删除rna counts数据正常样品
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rnaCounts.DEG=rnaCounts[,group==0]
colnames(rnaCounts.DEG)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(rnaCounts.DEG))
samesample<-intersect(colnames(rnaCounts.DEG),metaMatrix.RNA.DEG$patient)
rnaCounts.DEG<-rnaCounts.DEG[,samesample]
metaMatrix.RNA.DEG<-metaMatrix.RNA.DEG[metaMatrix.RNA.DEG$patient%in%samesample,]
group=sapply(strsplit(colnames(mirCounts),"\\-"),"[",4)#删除mir counts数据正常样品
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
mirCounts.DEG=mirCounts[,group==0]
colnames(mirCounts.DEG)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(mirCounts.DEG))
samesample<-intersect(colnames(mirCounts.DEG),metaMatrix.MIR.DEG$patient)
mirCounts.DEG<-mirCounts.DEG[,samesample]
metaMatrix.MIR.DEG<-metaMatrix.MIR.DEG[metaMatrix.MIR.DEG$patient%in%samesample,]

###转录组差异分析，用的是原始count数据
DEGAll <- gdcDEAnalysis(counts     = rnaCounts.DEG, 
                        group      = metaMatrix.RNA.DEG$MSI, 
                        comparison = '1-0', 
                        method     = 'DESeq2')

###miRNA差异分析
degMI <- gdcDEAnalysis(counts     = mirCounts.DEG, 
                        group      = metaMatrix.MIR.DEG$MSI, 
                        comparison = '1-0', 
                        method     = 'DESeq2')

#所有基因差异（采用的是count数据）
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = FCfilter.lnrna, pval = adjpFilter)

#所有基因miRNA,输出差异的miRNA（采用的是count数据）
deMI <- gdcDEReport(deg = degMI, gene.type = 'all', fc = FCfilter.mirna, pval = adjpFilter)
deMIout=cbind(row.names(deMI),deMI)
write.table(deMIout, file='miRNA.diff.txt', sep='\t', quote=F, row.names=F)
#miRNA火山图
allDiff=gdcDEReport(deg = degMI, gene.type = 'all', fc = 0, pval = 1)
Significant=ifelse((allDiff$FDR<adjpFilter & abs(allDiff$logFC)>logFCfilter.mirna), ifelse(allDiff$logFC>logFCfilter.mirna,"Up","Down"), "Not")
p = ggplot(allDiff, aes(logFC, -log10(FDR)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
pdf("miRNA.vol.pdf",width=5.5,height=5)
print(p)
dev.off()

#差异lncRNA
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = FCfilter.lnrna, pval = adjpFilter)
write.table(deLNC, file='lncRNA.diff.txt', sep='\t', quote=F, row.names=F)
#lncRNA火山图
allDiff=gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = 0, pval = 1)
Significant=ifelse((allDiff$FDR<adjpFilter & abs(allDiff$logFC)>logFCfilter.lnrna), ifelse(allDiff$logFC>logFCfilter.lnrna,"Up","Down"), "Not")
p = ggplot(allDiff, aes(logFC, -log10(FDR)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
pdf("lncRNA.vol.pdf",width=5.5,height=5)
print(p)
dev.off()

#差异mRNA
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = FCfilter.mrna, pval = adjpFilter)
write.table(dePC, file='mRNA.diff.txt', sep='\t', quote=F, row.names=F)
#mRNA火山图
allDiff=gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = 0, pval = 1)
Significant=ifelse((allDiff$FDR<adjpFilter & abs(allDiff$logFC)>logFCfilter.mrna), ifelse(allDiff$logFC>logFCfilter.mrna,"Up","Down"), "Not")
p = ggplot(allDiff, aes(logFC, -log10(FDR)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
pdf("mRNA.vol.pdf",width=5.5,height=5)
print(p)
dev.off()

###开始建立ceRNA网络关系（利用是矫正数据）
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(deMI),
                          lnc.targets = 'miRcode',    ###'spongeScan', 'starBase', and 'miRcode'
                          pc.targets  = 'miRcode',    ###'spongeScan', 'starBase', and 'miRcode'
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

#根据设置的条件进行ceRNA网络过滤
ceOutput2 <- ceOutput[ceOutput$hyperPValue<hyperPfilter & 
    ceOutput$corPValue<corPfilter & ceOutput$regSim != 0,]
write.table(ceOutput2, file='ceRNA.score.txt', sep='\t', quote=F, row.names=F) #ceRNA.score.txt文件即为cerna分析结果
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
edges=edges[which(edges[,2] %in% rownames(deMI)),]#需要注意的是分析给出的是所有的mirRNA,需要筛选差异表达的miRNA
nodes=nodes[which(nodes[,1] %in% c(as.vector(edges[,1]),rownames(deMI))),]#筛选差异表达的miRNA之后，要继续修改node文件。把没有和差异表达miRNA连接的node删
#输出ceRNA各个节点表达量用于后续的实验
sameSample=intersect(colnames(rnaExpr),colnames(mirExpr))
rnaExpr=rnaExpr[,sameSample]
mirExpr=mirExpr[,sameSample]
allExp=rbind(rnaExpr,mirExpr)
expOut=allExp[as.vector(nodes[,1]),]
expOut=cbind(Symbol=nodes[,2],expOut)
write.table(expOut, file="networkExp.txt", sep='\t', quote=F, row.names=F)
#END