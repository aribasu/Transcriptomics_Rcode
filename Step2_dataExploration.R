#goal of this script is to using multivariate statisical approaches to explore the structure of your data
#begin by taking a look at the text expression matrix you created at the end of the last class
head(filtered.matrix)

###############################################################################################
# set up your experimental design by reading in a targets file that explains treatments, conditions, etc
###############################################################################################
library(limma)
#read in a tab-delimited "targets" file with the study design
targets <- readTargets("sunyerMouse_studyDesign.txt", sep="\t")
targets
myGroups <- factor(paste(targets$cellType, targets$treatment, targets$phenotype, sep="."))
myGroups
design <- model.matrix(~0+myGroups)
colnames(design) <- levels(myGroups)
design


###############################################################################################
#carry out hierarchical clustering on filtered data
###############################################################################################
#make some sample labels 
sampleLabels <- paste(targets$cellType, targets$treatment, targets$phenotype, targets$rep, sep=".")
distance <- dist(t(filtered.matrix),method="maximum")
clusters <- hclust(distance, method = "complete") 
plot(clusters, label = sampleLabels)


###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
#assign distinct colors based on treatment
color.map.treatment <- function(myGroups) { 
  if (myGroups=="naive") "#FF0000" else 
    if (myGroups=="exposed") "#FF0000""#0000FF" }

cols.treatment <- unlist(lapply(myGroups, color.map.treatment))

pca.res <- prcomp(t(filtered.matrix), scale.=F, retx=T)
head(pca.res$rotation) #$rotation gives you the eigenvectors or 'loadings'
head(pca.res$x) #$x gives you the 'scores'
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per
plot.pch <- (1:length(levels(factorial)))[as.numeric(factorial)]
plot(pca.res$x, col=1, pch=plot.pch, las=1, cex=3, 
     xlab=paste("PC1 (", pc.per [1], "%)", sep=""),
     ylab=paste("PC2 (", pc.per [2], "%)", sep=""))
text(pca.res$x, sampleLabels.ALL, pos= 1 )
grid()
legend(-40, 0, levels(factorial), pch=1:length(levels(factorial)), pt.cex=1.5)
summary(pca.res) # Prints variance summary for all principal components.

#3d PCA plot of PC1, PC2 and PC3
library(rgl)
rgl.open() 
offset <- 50 
par3d(windowRect=c(offset, offset, 640+offset, 640+offset)) 
rm(offset); 
rgl.clear(); 
rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1); 
spheres3d(pca.res$x[,1], pca.res$x[,2], pca.res$x[,3], 
          color=cols.treatment, radius=5, alpha=1, shininess=20); 
aspect3d(1, 1, 1) 
axes3d(col='black') 
title3d("", "", "PC1", "PC2", "PC3", col='black')
bg3d("white") 

#take a look at the loadings for all the samples, examining one principal component at a time
# loadings <- as.data.frame(pca.res$x)
# library(ggplot2)
# loadings$var <- factor(targets$Name, as.character(targets$Name))
# qplot(x=var, y=PC3, data=loadings, geom="bar", stat="identity", fill=cols.treatment)
#note: PC1 is strongly influenced by infection vs naive, while PC2 is entirely due to mouse strain


