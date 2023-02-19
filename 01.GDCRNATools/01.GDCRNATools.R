rm(list = ls())#��չ����ռ�
library(GDCRNATools)#���ð�
library(ggplot2)#���ð�

adjpFilter=0.05          #�������FDR��ֵ
logFCfilter.mrna=1          #�������logFC��ֵmrna
FCfilter.mrna=2^logFCfilter.mrna
logFCfilter.mirna=0.5            #�������logFC��ֵmirna
FCfilter.mirna=2^logFCfilter.mirna
logFCfilter.lnrna=1            #�������logFC��ֵlnrna
FCfilter.lnrna=2^logFCfilter.lnrna
hyperPfilter=0.01        #�����ηֲ�����pvalue��ֵ
corPfilter=0.01          #����Լ���pvalue��ֵ
offline<-T              #��ȡTCGA-STAD�������ݵķ�ʽT����F��������
MSIfile<-'MSI-prediction-STAD.txt' #��ȡMSI����
setwd("D:\\Bioinformatics\\immCeRNA-STAD\\github\\01.GDCRNATools")        #���ù���Ŀ¼

cancer<-c("STAD")#�������ͣ������о��������޸�
project <- paste("TCGA",cancer,sep = "-")       
pros<-paste("TCGA",paste(cancer,collapse = ""),sep = "-")
rnadir <- paste(pros, 'RNAseq', sep='/')
mirdir <- paste(pros, 'miRNAs', sep='/')

###��ȡTCGA-STAD���ݿ���ѡ�����߻�����������ģʽ
if(offline==T){
    load("metaDataoffline.RData")
    adjpFilter=0.05          #�������FDR��ֵ
    logFCfilter.mrna=1          #�������logFC��ֵmrna
    FCfilter.mrna=2^logFCfilter.mrna
    logFCfilter.mirna=0.5            #�������logFC��ֵmirna
    FCfilter.mirna=2^logFCfilter.mirna
    logFCfilter.lnrna=1            #�������logFC��ֵlnrna
    FCfilter.lnrna=2^logFCfilter.lnrna
    hyperPfilter=0.01        #�����ηֲ�����pvalue��ֵ
    corPfilter=0.01          #����Լ���pvalue��ֵ
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
    
    #ת¼��metadata����
    metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
    metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
    meta.RNA.filter = dir(path=rnadir)#��ȡrnaseqĿ¼�е��ļ�����������meta�ļ�����
    metaMatrix.RNA<-metaMatrix.RNA[metaMatrix.RNA$file_id%in%meta.RNA.filter,]
    
    #��ȡmiRNA metadata
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
    
    #miRNA metadata����
    metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
    metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)
    meta.MIR.filter = dir(path=mirdir)
    metaMatrix.MIR<-metaMatrix.MIR[metaMatrix.MIR$file_id%in%meta.MIR.filter,]
    
    
    #ת¼�����ݺϲ�
    rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                             path      = rnadir, 
                             organized = FALSE,   ## if target data are in folders
                             data.type = 'RNAseq')
    
    #miRNA���ݺϲ�
    mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                             path      = mirdir,
                             organized = FALSE,   ## if target data are in folders
                             data.type = 'miRNAs')
    
    #ת¼�����ݽ���
    rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
    
    #miRNA���ݽ���
    mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
}

###��MSI��Ϣ����RNA��MRI��meta��Ϣ
MSIinfo<-read.table(file=MSIfile,check.names = F,header = T)
MSIinfo<-MSIinfo[!duplicated(MSIinfo$patient),]
metaMatrix.RNAall<-merge(metaMatrix.RNA,MSIinfo,by="patient",all=T)
metaMatrix.MIRall<-merge(metaMatrix.MIR,MSIinfo,by="patient",all=T)
metaMatrix.RNA<-metaMatrix.RNAall[metaMatrix.RNAall$patient%in%metaMatrix.RNA$patient,]
metaMatrix.MIR<-metaMatrix.MIRall[metaMatrix.MIRall$patient%in%metaMatrix.MIR$patient,]

###����na�Լ���������Ʒ
metaMatrix.RNA.DEG <- metaMatrix.RNA[!is.na(metaMatrix.RNA$MSI),]#na.omit��ɾ��������na���У���ɾ��msi��na��
metaMatrix.MIR.DEG<- metaMatrix.MIR[!is.na(metaMatrix.MIR$MSI),]
metaMatrix.RNA.DEG <- metaMatrix.RNA.DEG [!duplicated(metaMatrix.RNA.DEG$patient),]
metaMatrix.MIR.DEG<- metaMatrix.MIR.DEG[!duplicated(metaMatrix.MIR.DEG$patient),]
metaMatrix.RNA.DEG$MSI<-factor(metaMatrix.RNA.DEG$MSI)
metaMatrix.MIR.DEG$MSI<-factor(metaMatrix.MIR.DEG$MSI)
group=sapply(strsplit(colnames(rnaCounts),"\\-"),"[",4)#ɾ��rna counts����������Ʒ
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rnaCounts.DEG=rnaCounts[,group==0]
colnames(rnaCounts.DEG)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(rnaCounts.DEG))
samesample<-intersect(colnames(rnaCounts.DEG),metaMatrix.RNA.DEG$patient)
rnaCounts.DEG<-rnaCounts.DEG[,samesample]
metaMatrix.RNA.DEG<-metaMatrix.RNA.DEG[metaMatrix.RNA.DEG$patient%in%samesample,]
group=sapply(strsplit(colnames(mirCounts),"\\-"),"[",4)#ɾ��mir counts����������Ʒ
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
mirCounts.DEG=mirCounts[,group==0]
colnames(mirCounts.DEG)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(mirCounts.DEG))
samesample<-intersect(colnames(mirCounts.DEG),metaMatrix.MIR.DEG$patient)
mirCounts.DEG<-mirCounts.DEG[,samesample]
metaMatrix.MIR.DEG<-metaMatrix.MIR.DEG[metaMatrix.MIR.DEG$patient%in%samesample,]

###ת¼�����������õ���ԭʼcount����
DEGAll <- gdcDEAnalysis(counts     = rnaCounts.DEG, 
                        group      = metaMatrix.RNA.DEG$MSI, 
                        comparison = '1-0', 
                        method     = 'DESeq2')

###miRNA�������
degMI <- gdcDEAnalysis(counts     = mirCounts.DEG, 
                        group      = metaMatrix.MIR.DEG$MSI, 
                        comparison = '1-0', 
                        method     = 'DESeq2')

#���л�����죨���õ���count���ݣ�
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all', fc = FCfilter.lnrna, pval = adjpFilter)

#���л���miRNA,��������miRNA�����õ���count���ݣ�
deMI <- gdcDEReport(deg = degMI, gene.type = 'all', fc = FCfilter.mirna, pval = adjpFilter)
deMIout=cbind(row.names(deMI),deMI)
write.table(deMIout, file='miRNA.diff.txt', sep='\t', quote=F, row.names=F)
#miRNA��ɽͼ
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

#����lncRNA
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding', fc = FCfilter.lnrna, pval = adjpFilter)
write.table(deLNC, file='lncRNA.diff.txt', sep='\t', quote=F, row.names=F)
#lncRNA��ɽͼ
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

#����mRNA
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding', fc = FCfilter.mrna, pval = adjpFilter)
write.table(dePC, file='mRNA.diff.txt', sep='\t', quote=F, row.names=F)
#mRNA��ɽͼ
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

###��ʼ����ceRNA�����ϵ�������ǽ������ݣ�
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC),
                          deMIR       = rownames(deMI),
                          lnc.targets = 'miRcode',    ###'spongeScan', 'starBase', and 'miRcode'
                          pc.targets  = 'miRcode',    ###'spongeScan', 'starBase', and 'miRcode'
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)

#�������õ���������ceRNA�������
ceOutput2 <- ceOutput[ceOutput$hyperPValue<hyperPfilter & 
    ceOutput$corPValue<corPfilter & ceOutput$regSim != 0,]
write.table(ceOutput2, file='ceRNA.score.txt', sep='\t', quote=F, row.names=F) #ceRNA.score.txt�ļ���Ϊcerna�������
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
edges=edges[which(edges[,2] %in% rownames(deMI)),]#��Ҫע����Ƿ��������������е�mirRNA,��Ҫɸѡ��������miRNA
nodes=nodes[which(nodes[,1] %in% c(as.vector(edges[,1]),rownames(deMI))),]#ɸѡ��������miRNA֮��Ҫ�����޸�node�ļ�����û�кͲ������miRNA���ӵ�nodeɾ
#���ceRNA�����ڵ���������ں�����ʵ��
sameSample=intersect(colnames(rnaExpr),colnames(mirExpr))
rnaExpr=rnaExpr[,sameSample]
mirExpr=mirExpr[,sameSample]
allExp=rbind(rnaExpr,mirExpr)
expOut=allExp[as.vector(nodes[,1]),]
expOut=cbind(Symbol=nodes[,2],expOut)
write.table(expOut, file="networkExp.txt", sep='\t', quote=F, row.names=F)
#END