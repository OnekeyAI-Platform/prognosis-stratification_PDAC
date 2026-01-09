#library
library(devtools)
library(e1071)
library(preprocessCore)
library(parallel)
library(bseqsc)
library(ggplot2)
library(CIBERSORT)
library(corrplot)
library(vioplot)
library(limma)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggVolcano)

# DEGs
setwd("D:/R/PDAC-test3")
# read.tsvfile

#install.packages("tidyverse")
library(tidyverse)

# read data
data=read.table("CPTAC-3.star_tpm_filtered.tsv", header=T, sep="\t", check.names=F)

# group
group = read.table(file = 'groupby.txt', sep = '\t', header = TRUE )

# data clean
colnames(data) <- gsub("-", ".", colnames(data))
group$id <- gsub("-", ".", group$id)

rownames(data) = data[,1]
rownames(data) <- make.unique(data[,1])
data = data[,-1]

# retain 58 cases
keepSamps <- intersect(colnames(data), group$id)
data       <- data[,keepSamps]
group = group1
group <- group[match(keepSamps,group$id), ]

# matrix
exprList <- list(
  Low = tpm[group$groupby == "Low", , drop = FALSE],
  High  = tpm[group$groupby == "High",  , drop = FALSE]
)

rownames(group) = group[,1]
group = data.frame(Group = group[,-1])

# read clean data
groupname = read.table("groupby.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(group) = groupname[,1]
rownames(group) <- gsub("-", ".", rownames(group))

# statistics
rownames(group) = group[,1]
group$groupby <- gsub("High", "1", group$groupby)
group$groupby <- gsub("Low", "0", group$groupby)
grp <- group[,2]
names(grp) <- row.names(group)
table(grp)

lowNum  <- sum(grp == 0)
highNum <- sum(grp == 1)

Type=c(rep(1,lowNum), rep(2,highNum))
data1 <- tpm[ , grp == 0, drop = FALSE]
data2 <- tpm[ , grp == 1, drop = FALSE]
data  <- cbind(data1, data2)

str(data)
data = as.matrix(data)


#DEGs analysis
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

write.table(outTab,file="TCGA.all.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

logFCfilter=0.5
fdrFilter=0.05
#Output
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & 
                   as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="TCGA.diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)


#Heatmap
#gene display
geneNum=50     
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=6.5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#
pdf(file="vol.pdf", width=5, height=5)
xMax=6
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1.2)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=1.5)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=1.5)
abline(v=0,lty=2,lwd=3)
dev.off()


# CIBERSORT

data(LM22) 
data=read.table("CPTAC-3.star_tpm_filtered用于差异表达.tsv", header=T, sep="\t", check.names=F)
group = read.table(file = 'groupby.txt', sep = '\t', header = TRUE )
colnames(data) <- gsub("-", ".", colnames(data))
group$id <- gsub("-", ".", group$id)

rownames(data) = data[,1]
rownames(data) <- make.unique(data[,1])
data = data[,-1]

# 
keepSamps <- intersect(colnames(data), group$id)
data       <- data[,keepSamps]
group = group1
group <- group[match(keepSamps,group$id), ]


#
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#
results <- cibersort(sig_matrix = LM22, mixture_file = data)

#保存
results=as.matrix(results[,1:(ncol(results)-3)])
results=rbind(id=colnames(results),results)
write.table(results, file="CIBERSORT-Results1.txt", sep="\t", quote=F, col.names=F)

#
immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
immune=as.matrix(immune)
data=t(immune)


col=rainbow(nrow(data),s=0.7,v=0.7)
pdf("CIBERSORT1.pdf",height=10,width=22)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2 = axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()



cellnum <- read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', 
               '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', 
               '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', 
               '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', 
               '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', 
               '#3A6963','#968175')
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
  data4plot <- rbind(
    data4plot,
    cbind(cell.prop[,i],rownames(cell.prop),
          rep(colnames(cell.prop)[i],nrow(cell.prop)
          )
    )
  )
}
colnames(data4plot)<-c('proportion','celltype','sample')
data4plot$proportion <- as.numeric(data4plot$proportion)
pdf(file="CIBERSORT2.pdf",height=10,width=22)
ggplot(data4plot,aes(sample,proportion,fill=celltype))+
  geom_bar(stat="identity",position="fill")+
  scale_fill_manual(values=my36colors)+
  ggtitle("cell portation")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.title.x=element_text(size=1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  guides(fill=guide_legend(title=NULL))
dev.off()

rt=read.table("CIBERSORT-Results.txt", header=T, sep="\t", check.names=F, row.names=1)
group_info = read.table("group.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(group_info) = group_info[,1]
group_info = data.frame(Group = group_info[,-1])
groupname = read.table("group.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(group_info) = groupname[,1]
rownames(group_info) <- gsub("-", ".", rownames(group_info))

group = group_info[rownames(rt), "Group"]
LowNum = sum(group_info == 0)
HighNum = sum(group_info == 1)
Type = c(rep(1, LowNum), rep(2, HighNum))
rt1 = rt[group == 0, ]
rt2 = rt[group == 1, ]
rt = rbind(rt1, rt2)

outTab=data.frame()
pdf(file="vioplot1.pdf", width=13, height=8)
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x, y,
     xlim=c(0,63), ylim=c(min(rt), max(rt)+0.05),
     main="", xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

for(i in 1:ncol(rt)){
  if(sd(rt[1:conNum,i])==0){
    rt[1,i]=0.00001
  }
  if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
    rt[(conNum+1),i]=0.00001
  }
  rt1=rt[1:conNum,i]
  rt2=rt[(conNum+1):(conNum+treatNum),i]
  vioplot(rt1,at=3*(i-1),lty=1,add = T,col = '#ADD8E6')
  vioplot(rt2,at=3*(i-1)+1,lty=1,add = T,col = '#FFB6C1')
  wilcoxTest=wilcox.test(rt1,rt2)
  p=wilcoxTest$p.value
  if(p<0.05){
    cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
    outTab=rbind(outTab,cellPvalue)
  }
  mx=max(c(rt1,rt2))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
}
legend("topright", 
       c("Low Risk", "High Risk"),
       lwd=3,bty="n",cex=1,
       col=c("#ADD8E6","#FFB6C1"))
text(seq(1,64,3),-0.04,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
dev.off()
write.table(outTab,file="Diff1.xls",sep="\t",row.names=F,quote=F)



#GO,KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
R.utils::setOption("clusterProfiler.download.method","auto")

#
pvalueFilter=0.05  
qvalueFilter=0.5       


colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
ontology.col=c("#00AFBB", "#E7B800", "#90EE90")


rt=read.table("TCGA.diff.Wilcoxon.txt", header=T, sep="\t", check.names=F)    
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]

#GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)
showNum=10
if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GObarplot.pdf", width=10, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
pdf(file="GObubble.pdf", width=10, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

data=GO[order(GO$pvalue),]
datasig=data[data$pvalue<0.05,,drop=F]
BP = datasig[datasig$ONTOLOGY=="BP",,drop=F]
CC = datasig[datasig$ONTOLOGY=="CC",,drop=F]
MF = datasig[datasig$ONTOLOGY=="MF",,drop=F]
BP = head(BP,6)
CC = head(CC,6)
MF = head(MF,6)
data = rbind(BP,CC,MF)
main.col = ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(GO=data$ID,start=1,end=max(BgGene))
rownames(df) = df$GO
bed2 = data.frame(GO=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(GO=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(GO=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5


pdf("GO.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',ontology.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
main.legend = Legend(
  labels = c("Biological Process","Cellular Component", "Molecular Function"),  type="points",pch=15,
  legend_gp = gpar(col=ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()

#KEGG
rt=read.table("TCGA.diff.Wilcoxon.txt", header=T, sep="\t", check.names=F)     

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]

#KEGG
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

showNum=30
if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGGbarplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

pdf(file="KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=18
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
  data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

pdf(file="KEGG.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()




# GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)

diff_data <- read.delim("TCGA.diff.Wilcoxon.txt", header = TRUE, stringsAsFactors = FALSE)
head(diff_data)
str(diff_data)

gene_list <- diff_data$log2FC
names(gene_list) <- diff_data$gene_symbol
gene_list <- sort(gene_list, decreasing = TRUE)

gene_list <- gene_list[!is.na(gene_list)]

c2_gene_sets <- msigdbr(
  species = "Homo sapiens",  
  category = "C2"           
)


table(c2_gene_sets$gs_subcat)

gene_sets_list <- split(
  c2_gene_sets$entrez_gene, 
  c2_gene_sets$gs_name
)

gene_symbols <- names(gene_list)

gene_mapping <- mapIds(
  org.Hs.eg.db,
  keys = gene_symbols,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_list_entrez <- gene_list
names(gene_list_entrez) <- gene_mapping[names(gene_list)]

gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]

set.seed(123)  
gsea_result <- GSEA(
  geneList = gene_list_entrez,
  exponent = 1,          
  minGSSize = 10,        
  maxGSSize = 500,       
  pvalueCutoff = 0.05,   
  pAdjustMethod = "BH",  
  TERM2GENE = data.frame(
    term = c2_gene_sets$gs_name,
    gene = c2_gene_sets$entrez_gene
  ),
  verbose = TRUE,
  seed = TRUE
)

if (nrow(gsea_result) > 0) {
  head(gsea_result@result, 6)
  write.csv(gsea_result@result, "GSEA_c2_results.csv", row.names = FALSE)

  if (nrow(gsea_result) >= 5) {

    if (nrow(gsea_result) >= 3) {
      p3 <- gseaplot2(gsea_result, 
                      geneSetID = 1:3, 
                      pvalue_table = TRUE)
      ggsave("GSEA_multiple_pathways.png", p3, width = 12, height = 8, dpi = 300)
    }
    
    p4 <- ridgeplot(gsea_result, showCategory = 15) + 
      labs(x = "enrichment distribution")
    ggsave("GSEA_ridgeplot.png", p4, width = 10, height = 8, dpi = 300)
  }
  
} 
