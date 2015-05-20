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
myGroups <- factor(paste(targets$cellType, targets$phenotype, sep="."))
myGroups
design <- model.matrix(~0+myGroups)
colnames(design) <- levels(myGroups)
design


###############################################################################################
#carry out hierarchical clustering on filtered data
###############################################################################################
#make some sample labels 
sampleLabels <- paste(targets$cellType, targets$phenotype, targets$rep, sep=".")
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
sampleLabels <- as.character(paste(targets$cellType, targets$phenotype, targets$rep, sep="."))
colnames(myData) <- sampleLabels
head(myData)
geneSymbols <- row.names(myData)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(myData,
                 B1a.Naive.AVG = (B1a.naive.1 + B1a.naive.2)/2,
                 B1b.Naive.AVG = (B1b.naive.1 + B1b.naive.2)/2,
                 Mac.Naive.AVG = (Mac.naive.1 + Mac.naive.2)/2,
                 B1a.Ph0.AVG = (B1a.Ph0.1 + B1a.Ph0.2)/2,
                 B1a.nonPh0.AVG = (B1a.nonPh0.1 + B1a.nonPh0.2)/2,
                 B1b.Ph0.AVG = (B1b.Ph0.1 + B1b.Ph0.2)/2,
                 B1b.nonPh0.AVG = (B1b.nonPh0.1 + B1b.nonPh0.2)/2,
                 Mac.Ph0.AVG = (Mac.Ph0.1 + Mac.Ph0.2)/2,
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

#plot a simple bar graph
ggplot(myData.filter.transpose, aes(y=Cxcl10)) +
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
library(ggvis)
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

#the goal of this script is to identify differentially expressed genes (DEG)
#you should already know what pairwise comparisons are most important to you

#I prefer to use a table of data with a non-redundant set of gene identifiers as input for this script
#this keeps my heatmaps and lists of DEGS as simple of as possible
#if you are coming to this script with array data, you've already reduced your list in this way
#RNAseq data is a bit more complicated (many transcripts, but the same gene), 
#but you achieve the same level of reduction by filtering based on annotation data
#################################################################
#additional filtering based on annotation
#if you have array data, you can skip this section
#################################################################
#our initial filtering of RNAseq data was based on solely cpm, which only removes relatively lowly expressed genes
#now that you have annotation info, you can further filter to reduce to single line of data per gene symbo using one of two options
#first, check out how many unqiue genes are represented in your data
dupFiltered <- unique(resultTable.filtered$symbol)
#use collapseRows function from WGCNA package to collapse your dataset
library(WGCNA)
#pull your rownames and unique identifiers
myIDs <- rownames(resultTable.filtered)
#retrieve your gene symbols from your data
mySymbols <- resultTable.filtered[,5]
#remove all annotation columns so you're left with only numeric data
resultTable <- resultTable.filtered[,-1:-6]
myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")
myCollapsed <- myCollapsed$datETcollapsed
#now that the matrix is collapsed to give non-redudant list of genes,
#you could set the symbols to be the row names and move on

#################################################################################################################
#if you have no biological replicates, you will not be able to leverage statistical tools to identify DE genes
#Instead, you will ONLY rely on fold changes
#if you DO have replicates, skip this section and proceed to the next part 
#################################################################################################################
#use the dplyr 'filter' command to capture all the genes that are up/down regulated x-fold in n conditions
#in this case, 'myData' is a dataframe that you generated with Log2 expression and annotation
myData.filter <- myData %>%
  filter((abs(Ecdysone.vs.PBS_18hr_gut) >= 1) | (abs(Ecdysone.vs.PBS_5hr_gut) >= 1)) %>%
  select(geneID, Ecdysone.vs.PBS_5hr_carcass, Ecdysone.vs.PBS_18hr_carcass, Ecdysone.vs.PBS_5hr_gut, Ecdysone.vs.PBS_18hr_gut)
head(myData.filter)


###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to your filtered expression data
library(limma)
fit <- lmFit(filtered.matrix, design)
#add annotation into your linear model fit
#don't really need to do this if you have RNAseq data
library(annotate)
fit$genes$Symbol <- getSYMBOL(probeList, "lumiMouseAll.db")
fit$genes$Entrez <- getEG(probeList, "lumiMouseAll.db")

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix.naive <- makeContrasts(B1a_vs_B1b = B1a.naive - B1b.naive, Mac_vs_B1a = Mac.naive - B1a.naive, Mac_vs_B1b = Mac.naive - B1b.naive, levels=design)
contrast.matrix.Mac <- makeContrasts(Exposed_vs_naive = Mac.Ph0 - Mac.naive, Exposed_vs_B1aPh0 = Mac.Ph0 - B1a.Ph0, Exposed_vs_B1bPh0 = Mac.Ph0 - B1b.Ph0, levels=design)
contrast.matrix.B1a <- makeContrasts(Ph0_vs_naive = B1a.Ph0 - B1a.naive, nonPh0_vs_naive = B1a.nonPh0 - B1a.naive, Ph0_vs_nonPh0 = B1a.Ph0 - B1a.nonPh0, levels=design)
contrast.matrix.B1b <- makeContrasts(Ph0_vs_naive = B1b.Ph0 - B1b.naive, nonPh0_vs_naive = B1b.nonPh0 - B1b.naive, Ph0_vs_nonPh0 = B1b.Ph0 - B1b.nonPh0, levels=design)
contrast.matrix.B1a_and_B1b <- makeContrasts(Ph0_vs_naive = B1a.Ph0 - B1a.naive, nonPh0_vs_naive = B1a.nonPh0 - B1a.naive, Ph0_vs_nonPh0 = B1a.Ph0 - B1a.nonPh0, Ph0_vs_naive = B1b.Ph0 - B1b.naive, nonPh0_vs_naive = B1b.nonPh0 - B1b.naive, Ph0_vs_nonPh0 = B1b.Ph0 - B1b.nonPh0, levels=design)
contrast.matrix.Ph0 <- makeContrasts(B1a_vs_B1b = B1a.Ph0 - B1b.Ph0, levels=design)
contrast.matrix.nonPh0 <- makeContrasts(B1a_vs_B1b = B1a.nonPh0 - B1b.nonPh0, levels=design)

# check each contrast matrix
contrast.matrix.naive
contrast.matrix.Mac
contrast.matrix.B1a
contrast.matrix.B1b
contrast.matrix.B1a_and_B1b
contrast.matrix.Ph0
contrast.matrix.nonPh0

# extract the linear model fit for the contrast matrix that you just defined above
fits.naive <- contrasts.fit(fit, contrast.matrix.naive)
fits.Mac <- contrasts.fit(fit, contrast.matrix.Mac)
fits.B1a <- contrasts.fit(fit, contrast.matrix.B1a)
fits.B1b <- contrasts.fit(fit, contrast.matrix.B1b)
fits.B1a_and_B1b <- contrasts.fit(fit, contrast.matrix.B1a_and_B1b)
fits.Ph0 <- contrasts.fit(fit, contrast.matrix.Ph0)
fits.nonPh0 <- contrasts.fit(fit, contrast.matrix.nonPh0)

ebFit.naive <- eBayes(fits.naive)
ebFit.Mac <- eBayes(fits.Mac)
ebFit.B1a <- eBayes(fits.B1a)
ebFit.B1b <- eBayes(fits.B1b)
ebFit.B1a_and_B1b <- eBayes(fits.B1a_and_B1b)
ebFit.Ph0 <- eBayes(fits.Ph0)
ebFit.nonPh0 <- eBayes(fits.nonPh0)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
probeList.naive <- topTable(ebFit.naive, adjust ="BH", coef=1, number=50, sort.by="logFC")
probeList.naive
probeList.Mac <- topTable(ebFit.Mac, adjust ="BH", coef=2, number=25, sort.by="logFC")
probeList.Mac
probeList.B1a <- topTable(ebFit.B1a, adjust ="BH", coef=3, number=100, sort.by="logFC")
probeList.B1a
probeList.B1b <- topTable(ebFit.B1b, adjust ="BH", coef=2, number=100, sort.by="logFC")
probeList.B1b
probeList.Ph0 <- topTable(ebFit.Ph0, adjust ="BH", coef=1, number=25, sort.by="logFC")
probeList.Ph0
probeList.nonPh0 <- topTable(ebFit.nonPh0, adjust ="BH", coef=1, number=25, sort.by="logFC")
probeList.nonPh0

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
results.naive <- decideTests(ebFit.naive, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.Mac <- decideTests(ebFit.Mac, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.B1a <- decideTests(ebFit.B1a, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.B1b <- decideTests(ebFit.B1b, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.B1a_and_B1b <- decideTests(ebFit.B1a_and_B1b, method="global", adjust.method="BH", p.value=0.01, lfc=1)
results.Ph0 <- decideTests(ebFit.Ph0, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
results.nonPh0 <- decideTests(ebFit.nonPh0, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)

vennDiagram(results.naive, include="up") #all pairwise comparisons on a B6 background

# retrieve gene symbols for the probes from above
diffSymbols.naive <- fit$genes$Symbol[results.naive[,1] !=0 | results.naive[,2] !=0 | results.naive[,3] !=0]

# retrieve entrez IDs for the probes from above
diffEntrez.naive <- fit$genes$Entrez[results.naive[,1] !=0 | results.naive[,2] !=0 | results.naive[,3] !=0]

# retrieve expression data for the probes from above
filtered.matrix.eset <- new("ExpressionSet", exprs = filtered.matrix)
diffData.naive <- filtered.matrix.eset[results.naive[,1] !=0 | results.naive[,2] !=0 | results.naive[,3] !=0]
diffData.naive <- exprs(diffData.naive)
#write out your differentially expressed genes
write.table(cbind(diffSymbols.naive, diffEntrez.naive, diffData.naive), "myDEGs_naive.xls", sep="\t", quote=FALSE)







