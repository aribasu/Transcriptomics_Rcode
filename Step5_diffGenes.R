###############################################################################################
# use Limma to find differentially expressed genes between two or more conditions
###############################################################################################
# fit the linear model to the filtered expression set
library(limma)
fit <- lmFit(filtered.matrix, design)
library(annotate)
fit$genes$Symbol <- getSYMBOL(probeList, "lumiMouseAll.db")
fit$genes$Entrez <- getEG(probeList, "lumiMouseAll.db")

# set up a contrast matrix based on the pairwise comparisons of interest
contrast.matrix <- makeContrasts(BM = boneMarrow.Treg - boneMarrow.Tcon, SP = spleen.Treg - spleen.Tcon, BMvsSP = boneMarrow.Treg - spleen.Treg, levels=design)

# check each contrast matrix
contrast.matrix

# extract the linear model fit for the contrast matrix that you just defined above
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)


###############################################################################################
# use the topTable and decideTests functions to see the differentially expressed genes
###############################################################################################

# use topTable function to take a look at the hits
probeList <- topTable(ebFit, adjust ="BH", coef=3, number=20, sort.by="logFC")
probeList

# use the 'decideTests' function to show Venn diagram for all diffexp genes for up to three comparisons
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.59)

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

