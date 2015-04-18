#goal of this script is to using multivariate statisical approaches to explore the structure of your data

###############################################################################################
# set up your experimental design by reading in a targets file that explains treatments, conditions, etc
###############################################################################################
library(limma)
#read in a tab-delimited "targets" file with the study design
targets <- readTargets("balelab_studyDesign_miR.txt", sep="\t")
targets
factorial <- factor(paste(targets$tissue, targets$cellType, sep="."))
factorial
design <- model.matrix(~0+factorial)
colnames(design) <- levels(factorial)
design#all samples
#sampleLabels <- read.delim("studyDesign.txt", sep="\t", stringsAsFactors = FALSE)
#sampleLabels.ALL <- paste(sampleLabels$treatment, sampleLabels$mouse, sep=".")
sampleLabels.ALL <- c("naive.1", "naive.2", "naive.3", "naive.4", "memory.1", "memory.2", "memory.3", "memory.4")
distance.ALL <- dist(t(filtered.matrix),method="maximum")
clusters.ALL <- hclust(distance.ALL, method = "complete") 
plot(clusters.ALL, label = sampleLabels.ALL)


###############################################################################################
#carry out hierarchical clustering on filtered data
###############################################################################################
#all samples
#sampleLabels <- read.delim("AriTreg_studyDesign.txt", sep="\t", stringsAsFactors = FALSE)
#sampleLabels.ALL <- paste(sampleLabels$treatment, sampleLabels$mouse, sep=".")
sampleLabels.ALL <- c("naive.1", "naive.2", "naive.3", "naive.4", "memory.1", "memory.2", "memory.3", "memory.4")
distance.ALL <- dist(t(filtered.matrix),method="maximum")
clusters.ALL <- hclust(distance.ALL, method = "complete") 
plot(clusters.ALL, label = sampleLabels.ALL)




#FIX THIS
#assign distinct colors based on treatment
# color.map.treatment <- function(factorial) { if (factorial=="naive") "#FF0000" else "#0000FF" }
# cols.treatment <- unlist(lapply(factorial, color.map.treatment))

###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
pca.res <- prcomp(t(filtered.matrix), scale.=F, retx=T)
head(pca.res$rotation) #$rotation gives you the eigenvectors or 'loadings'
head(pca.res$x) #$x gives you the 'scores'
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per
plot.pch <- (1:length(levels(factorial)))[as.numeric(factorial)]
plot(pca.res$x, col=1, pch=plot.pch, las=1, cex=3, xlab=paste("PC1 (", pc.per [1], "%)", sep=""),ylab=paste("PC2 (", pc.per [2], "%)", sep=""))
text(pca.res$x, sampleLabels.ALL, pos= 1 )
grid()
legend(-40, 0, levels(factorial), pch=1:length(levels(factorial)), pt.cex=1.5)
summary(pca.res) # Prints variance summary for all principal components.
#3d PCA plot of PC1, PC2 and PC3
library(rgl)
rgl.open(); offset <- 50; par3d(windowRect=c(offset, offset, 640+offset, 640+offset)); rm(offset); rgl.clear(); rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1); spheres3d(pca.res$x[,1], pca.res$x[,2], pca.res$x[,3], color=cols.treatment, radius=5, alpha=1, shininess=20); aspect3d(1, 1, 1); axes3d(col='black'); title3d("", "", "PC1", "PC2", "PC3", col='black'); bg3d("white") 

#take a look at the loadings for all the samples, examining one principal component at a time
# loadings <- as.data.frame(pca.res$x)
# library(ggplot2)
# loadings$var <- factor(targets$Name, as.character(targets$Name))
# qplot(x=var, y=PC3, data=loadings, geom="bar", stat="identity", fill=cols.treatment)
#note: PC1 is strongly influenced by infection vs naive, while PC2 is entirely due to mouse strain


