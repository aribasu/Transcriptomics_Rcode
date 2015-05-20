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
contrast.matrix <- makeContrasts(BM = boneMarrow.Treg - boneMarrow.Tcon, SP = spleen.Treg - spleen.Tcon, BMvsSP = boneMarrow.Treg - spleen.Treg, levels=design)

# check each contrast matrix
contrast.matrix

# extract the linear model fit for the contrast matrix that you just defined above
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
probeList <- topTable(ebFit, adjust ="BH", coef=3, number=20, sort.by="logFC")
probeList

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)
#stats <- write.fit(ebFit)
vennDiagram(results, include="up") #all pairwise comparisons on a B6 background


# take a look at what the results of decideTests looks like
results

# now pull out probeIDs from selected regions of the Venn diagram.  In this case, I want all genes in the venn.
diffProbes <- which(results[,1] !=0 | results[,2] !=0 | results[,3] !=0)
diffSymbols <- fit$genes$Symbol[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]
diffEntrez <- fit$genes$Entrez[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]

#before pulling out expression data for differentially expressed genes, convert matrix to eset with annotation
library(Biobase)
myEset.ALL <- new("ExpressionSet", exprs = filtered.matrix)
annotation(myEset.ALL) <- "lumiMouseAll.db"

# retrieve expression data for the probes from above
diffData <- filtered.eset[results[,1] !=0 | results[,2] !=0 | results[,3] !=0]


#pull the expression data back out of the eset object
diffData <- exprs(diffData)

#combine probeIDs, gene symbols and expression data for differentially expressed genes into one file
write.table(cbind(diffProbes, diffSymbols, diffEntrez, diffData),"DiffGenes.xls", sep="\t", quote=FALSE)

# take a look at each expression matrix
dim(diffData)



