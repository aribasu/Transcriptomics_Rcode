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
plot(clusters, label = sampleLabels, hang = -1)


###############################################################################################
#Principal component analysis of the filtered data matrix
###############################################################################################
pca.res <- prcomp(t(filtered.matrix), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
head(pca.res$rotation) #$rotation shows you how much each GENE influenced each PC (callled 'eigenvalues', or loadings)
head(pca.res$x) #$x shows you how much each SAMPLE influenced each PC (called 'scores')
plot(pca.res, las=1)
pc.var<-pca.res$sdev^2 #sdev^2 gives you the eigenvalues
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per

#make some graphs to visualize your PCA result
library(ggplot2)
#lets first plot any two PCs against each other
#turn your scores for each gene into a data frame
data.frame <- as.data.frame(pca.res$x)
ggplot(data.frame, aes(x=PC1, y=PC2, colour=factor(myGroups))) +
  geom_point(size=5) +
  theme(legend.position="right")

#create a 'small multiples' chart to look at impact of each variable on each pricipal component
library(reshape2)
melted <- cbind(myGroups, melt(pca.res$x[,1:9]))
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=myGroups), stat="identity") +
  facet_wrap(~Var2)


