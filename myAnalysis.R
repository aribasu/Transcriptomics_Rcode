# for this script to work, you will need to output two files from Illumima's GenomeStudio software 
# the first file should contain non-normalized, non-background subtracted probe-level data with columns for at least the following: Avg_Signal, Avg_NBeads, BEAD_STD and Detection Pval
# the second file should contain probe set data for controls

# being by loading the lumi package
# if using Mouse arrays, change 'Mouse' to 'Mouse' for the mapping and .db packages below
library(lumi)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)

###############################################################################################
# Read in the GenomeStudio output files using the 'LumiR' and 'addControlData2lumi' functions
# This will create 'LumiBatch' objects from the control and experimental probe data
###############################################################################################
rawData <- lumiR("FinalReport_samples_probes_NoNorm_NoBkrnd.txt", convertNuID = TRUE, sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = "lumiMouseIDMapping")
rawData

# Read control probe data into a separate LumiBatch and take a look at these controls
rawData <- addControlData2lumi("FinalReport_controls_probes_NoNorm_NoBkrnd.txt", rawData)
rawData
controlData <- getControlData(rawData)
getControlType(controlData)
plotStringencyGene(rawData, lib = NULL, slideIndex = NULL, addLegend = TRUE, logMode = TRUE)
plotControlData(rawData, type = 'HOUSEKEEPING', slideIndex = NULL, logMode = TRUE, new = TRUE)
plotHousekeepingGene(rawData)

###############################################################################################
# QC check of data
###############################################################################################
summary(rawData, 'QC')

# choose color scheme for graphs
cols.ALL <- topo.colors(n=16, alpha=1)

hist(rawData, xlab = "log2 expression", main = "non-normalized data - histograms")
boxplot(rawData, ylab = "non-normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)
# pairs(rawData)
# MAplot(rawData)


###############################################################################################
# carry out all preprocessing steps (lumiQ, lumiB, lumiN, lumiT) using the 'lumiExpresso' function which encapsulates all these preprocessing steps
###############################################################################################
# options for normalization in Lumi include: "loess', 'quantile', 'vsn' and the default, which is robust spline normalization ('rsn') which combines the features of quantile and loess
normData <- lumiExpresso(rawData, QC.evaluation=TRUE, normalize.param=list(method='rsn'))


###############################################################################################
# repeat the summarization of QC data and graphs on the normalized data
###############################################################################################
summary(normData, 'QC')
hist(normData, xlab = "log2 expression", main = "non-normalized data - histograms", col=cols.ALL)
boxplot(normData, ylab = "non-normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

######################################################
# filter out probes that don't change (low variance)
# remove duplicate genes based on entrez ID
# remove coding sequences with no entrezIDs
######################################################
library(genefilter)
filtered_geneList <- nsFilter(normData, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE)
head(filtered_geneList)
# extract the ExpressionSet from this filtered list
filtered.eset <- filtered_geneList$eset
head(filtered.eset)

# now convert to a datamatrix that will contain only the probes after filtering
filtered.matrix <- as.matrix(filtered.eset)
probeList <- rownames(filtered.matrix)
head(filtered.matrix)


###############################################################################################
# output this filtered data to a table that can be used in Excel or other programs
# This file contains genes after normalization and filtering
###############################################################################################
# first need to get all the gene symbols and entrez IDs so you can put everything together
library(annotate)
keytypes(lumiMouseAll.db)
myGenesAll <- getSYMBOL(probeList, "lumiMouseAll.db")
myGenesAll <- as.matrix(myGenesAll)
myEntrezAll <- getEG(probeList, "lumiMouseAll.db")
myEntrezAll <- as.matrix(myEntrezAll)
write.table(cbind(myGenesAll, myEntrezAll, filtered.matrix),"normalizedFilteredData.txt", sep="\t", quote=FALSE)

###############################################################################################
# also output unfiltered expression data.  
###############################################################################################
write.exprs(normData, file='robustSplineNorm.txt')
normData.delim <- read.delim("robustSplineNorm.txt", header=TRUE, row.names=1)
normData.matrix <- as.matrix(normData.delim)
probeList2 <- rownames(normData.matrix)
library(annotate)
symbols <- getSYMBOL(probeList2, "lumiMouseAll.db")
entrezIDs <- getEG(probeList2, "lumiMouseAll.db")
write.table(cbind(symbols, entrezIDs, normData.matrix), "unfiltered_expressionData.txt", sep="\t", quote=FALSE)

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
distance <- dist(t(filtered.matrix),method="euclidean")
clusters <- hclust(distance, method = "average") 
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
melted <- cbind(myGroups, melt(pca.res$x[,1:3]))
#look at your 'melted' data
ggplot(melted) +
  geom_bar(aes(x=Var1, y=value, fill=myGroups), stat="identity") +
  facet_wrap(~Var2)


#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#read in your data from a text file that contains genes symbols as rows, and samples as columns. 
myData <- read.delim("normalizedFilteredData.txt", header=TRUE)
head(myData)
myData <- myData[,-2]
row.names(myData) <- myData[,1]
myData <- myData[,-1]
head(myData)

#column headers a bit cumbersome, so we'll change these to something more human-readable
targets <- read.delim("SunyerMouse_studyDesign.txt", sep="\t", stringsAsFactors = FALSE)
sampleLabels <- as.character(paste(targets$cellType, targets$treatment, targets$phenotype, targets$rep, sep="."))
colnames(myData) <- sampleLabels
head(myData)
geneSymbols <- row.names(myData)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(myData,
                 B1a.Naive.AVG = (B1a.naive.naive.1 + B1a.naive.naive.2)/2,
                 B1b.Naive.AVG = (B1b.naive.naive.1 + B1b.naive.naive.2)/2,
                 Mac.Naive.AVG = (Mac.naive.naive.1 + Mac.naive.naive.2)/2,
                 B1a.Ph0.AVG = (B1a.exposed.Ph0.1 + B1a.exposed.Ph0.2)/2,
                 B1a.nonPh0.AVG = (B1a.exposed.nonPh0.1 + B1a.exposed.nonPh0.2)/2,
                 B1b.Ph0.AVG = (B1b.exposed.Ph0.1 + B1b.exposed.Ph0.2)/2,
                 B1b.nonPh0.AVG = (B1b.exposed.nonPh0.1 + B1b.exposed.nonPh0.2)/2,
                 Mac.Ph0.AVG = (Mac.exposed.Ph0.1 + Mac.exposed.Ph0.2)/2,
                 LogFC.B1b.vs.B1a.naive = (B1b.Naive.AVG - B1a.Naive.AVG),
                 LogFC.B1b.vs.Mac.naive = (B1b.Naive.AVG - Mac.Naive.AVG),
                 LogFC.B1a.Ph0.vs.B1a.nonPh0 = (B1a.Ph0.AVG - B1a.nonPh0.AVG),
                 LogFC.B1b.Ph0.vs.B1b.nonPh0 = (B1b.Ph0.AVG - B1b.nonPh0.AVG),
                 LogFC.B1a.Ph0.vs.B1b.Ph0 = (B1a.Ph0.AVG - B1b.Ph0.AVG),
                 geneSymbols = geneSymbols)

#now look at this modified data table
row.names(myData) <- geneSymbols
head(myData)

#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange)
#and then display only the columns of interest (select) to see the most differentially expressed genes
myData.sort <- myData %>%
  arrange(desc(LogFC.B1a.Ph0.vs.B1b.Ph0)) %>%
  select(geneSymbols, LogFC.B1a.Ph0.vs.B1b.Ph0)
head(myData.sort)
class(geneSymbols)
#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#ways to tweek the 'select' function
#use : between two column names to select all columns between
#use 'contains', 'starts_with' or 'ends_with' to modify how you select
#can refer to columns using exact name or numerical indicator
myData.filter <- myData %>%
  filter(geneSymbols=="Ccl3" | geneSymbols=="Ccl4") %>%
  select(geneSymbols, B1a.Naive.AVG, B1b.Naive.AVG, Mac.Naive.AVG) 
head(myData.filter)

#another example using grep to match patterns
myData.filter <- myData %>%
  filter(grepl('Cxcl', geneSymbols)) %>%
  select(geneSymbols, B1a.Naive.AVG, B1b.Naive.AVG, Mac.Naive.AVG) 
myData.filter

#first reorder and clean up the data you want to graph
row.names(myData.filter) <- myData.filter[,1]
myData.filter <- select(myData.filter, -geneSymbols)
myData.filter.transpose <- as.data.frame(t(myData.filter))
cellType <- c("B1a", "B1b", "Mac")
myData.filter.transpose <- mutate(myData.filter.transpose, cellType, 
                                  phenotype=c("naive", "naive","naive"))
myData.filter.transpose

#another way to set up for graphing...don't filter before transposing
myData.graph <- select(myData,-geneSymbols)
head(myData.graph)
myData.graph.transpose <- as.data.frame(t(myData.graph))
dim(myData.graph.transpose)

#plot a simple bar graph
ggplot(myData.filter.transpose, aes(y=Ccl3)) +
  geom_bar(aes(x=cellType, fill=phenotype), stat="identity") +
  theme(axis.text.x=element_text(angle=-45))

#create a basic scatterplot using ggplot
ggplot(myData, aes(x=B1a.Ph0.AVG, y=Mac.Ph0.AVG)) +
  geom_point(shape=1) +
  geom_point(size=4)

# ##Volcano plots
# ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
# gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)
# 
# ##Construct the volcano plot
# ggplot(data=myData, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
#   geom_point(alpha=0.4, size=1.75) +
#   opts(legend.position = "none") +
#   xlim(c(-10, 10)) + ylim(c(0, 15)) +
#   xlab("log2 fold change") + ylab("-log10 p-value")

#define a tooltip that shows gene symbol and Log2 expression data when you mouse over each data point in the plot
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbols, "</b><br>",
         "B1a.Ph0.AVG: ", data$B1a.Ph0.AVG, "<br>",
         "B1a.nonPh0.AVG: ", data$B1a.nonPh0.AVG)
}

#plot the interactive graphic
myData %>% 
  ggvis(x= ~B1a.Ph0.AVG, y= ~B1a.nonPh0.AVG, key := ~geneSymbols) %>% 
  layer_points(fill = ~LogFC.B1a.Ph0.vs.B1a.nonPh0) %>%
  add_tooltip(tooltip)



