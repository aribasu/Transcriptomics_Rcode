#this script creates detailed heatmaps from your differentially expressed genes
#first part of script contains basic heatmap creation by reading in a file of a priori genes
#second part of script uses all the differentially expressed genes you identified previously

###################################################################
#first, make a basic heatmap from a small set of gene expression data
###################################################################
myData <- read.delim("balelab_miR_diff2foldFDR.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE, row.names=1)
class(myData)
head(myData)
myData.matrix <- as.matrix(myData)
#carry out hclust on the collapsed data matrix to generate a distance matrix for clustering
hr <- hclust(as.dist(1-cor(t(myData.matrix), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(myData.matrix, method="spearman")), method="average") #cluster columns by spearman correlation

# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force k number of clusters
mycl <- cutree(hr, k=2)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 

#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
library(gplots)
myheatcol <- greenred(75)
library(RColorBrewer)

#plot the hclust results as a heatmap
heatmap.2(myData.matrix, Rowv=as.dendrogram(hr), Colv=NA, col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc, cexRow=1.5, cexCol=1, margins=c(10,10), key = T, keysize = 1.05) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.



###################################################################################################
# generate a heatmap of differentially expressed transcripts using the entire dataset (B6 and balb)
###################################################################################################
#averaging the replicate arrays for all the data, then heatmap script below
head(diffData)
library(limma)
colnames(diffData) <- factorial
rownames(diffData) <- diffSymbols
head(diffData)
diffData.AVG <- avearrays(diffData)
head(diffData.AVG)
dim(diffData.AVG)

hr <- hclust(as.dist(1-cor(t(diffData.AVG), method="pearson")), method="complete") #cluster rows by pearson correlation
hc <- hclust(as.dist(1-cor(diffData.AVG, method="spearman")), method="complete") #cluster columns by spearman correlation
# Cut the resulting tree and create color vector for clusters.  Vary the cut height to give more or fewer clusters, or you the 'k' argument to force n number of clusters
mycl <- cutree(hr, k=7)
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9) 
mycolhc <- mycolhc[as.vector(mycl)] 
#load the gplots package for plotting the heatmap
library(gplots) 
#assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]).
myheatcol <- greenred(75)
library(RColorBrewer)

#plot the hclust results as a heatmap
heatmap.2(diffData.AVG, Rowv=as.dendrogram(hr), Colv=NA, col=myheatcol, scale="row", labRow=NA, density.info="none", trace="none", RowSideColors=mycolhc, cexRow=1, cexCol=1, margins=c(8,30)) # Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) # Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.

###############################################################################################
#select sub-clusters of co-regulated transcripts for downstream analysis
###############################################################################################
clid <- c(7); ysub <- diffData.AVG[names(mycl[mycl%in%clid]),]; hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
clusterIDs <- data.frame(Labels=rev(hrsub$labels[hrsub$order]))
clusterIDs <- as.vector(t(clusterIDs))
heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=NA, col=myheatcol, scale="row", labRow=getSYMBOL(clusterIDs, "lumiMouseAll.db"), labCol=sampleLabels.ALL, density.info="none", trace="none", RowSideColors=mycolhc[mycl%in%clid], margins=c(8,30)) # Create heatmap for chosen sub-cluster.

#retrieve gene symbols and entrezIDs for selected cluster and print out to an excel spreadsheet for downstream applications (i.e. GO enrichment in DAVID)
myCluster <- cbind(getSYMBOL(clusterIDs, "lumiMouseAll.db"), getEG(clusterIDs, "lumiMouseAll.db"))
write.table(myCluster, "Cluster7.xls", sep="\t", quote=FALSE)

