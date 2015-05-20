#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

######
#use the collapseRows function from the WGCNA package to get one gene symbol per line
library(WGCNA)
myIDs <- rownames(resultTable.filtered)
resultTable <- resultTable.filtered[,-1:-5]
mySymbols <- resultTable.filtered[,6]
myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")

dupFiltered <- unique(myData$Ripk3_Casp8.untreated.1)
myData.dupFiltered <- myData[dupFiltered,]
#read in your data from a text file that contains genes symbols as rows, and samples as columns. 
myData <- read.delim("normalizedFiltered.txt", header=TRUE)
head(myData)
myData <- myData[,-2]
row.names(myData) <- myData[,2]
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


