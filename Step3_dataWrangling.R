#this script walks thorough some basic data wrangling for organizing expression data spreadsheets and ends with
#how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#read in your data from a text file that contains genes (or probe IDs) as rows, and samples as columns. 
myData <- read.delim("normalizedFilteredData.txt", header=TRUE, row.names=1)
myData
#If you prefer, this data could be converted to a 'local' dataframe using the dplyr "tbl_df" function
#this simply makes reading the table a bit easier
#myData.local <- tbl_df(myData)
#myData.local

#column headers a bit cumbersome, so we'll change these to something more human-readable
sampleLabels <- read.delim("studyDesign.txt", sep="\t", stringsAsFactors = FALSE)
sampleLabels.ALL <- as.character(paste(sampleLabels$genotype, sampleLabels$treatment, sampleLabels$replicate, sep="."))
colnames(myData) <- sampleLabels.ALL
myData
geneSymbols <- row.names(myData)

#use the dplyr 'mutate' command to get averages and fold changes for all your replicates
myData <- mutate(myData,
       my.AVG.cond1 = (my.sample.rep1 + my.sample.rep2)/2,
       my.LogFC = (my.AVG.cond1 - my.AVG.cond2),
       geneSymbols = geneSymbols)

#use dplyr "arrange" and "select" functions to sort by LogFC column of interest (arrange) 
#and then display only the columns of interest (select) to see the most differentially expressed genes
myData.sort <- myData %>%
  arrange(desc(my.LogFC)) %>%
  select(geneSymbols, my.LogFC)
head(myData.sort)

#use dplyr "filter" and "select" functions to pick out genes of interest (filter)
#and again display only columns of interest (select)
myData.filter <- myData %>%
  filter(geneSymbols=="Il1a" | geneSymbols=="Casp8") %>%
  select(geneSymbols, my.AVG.cond1, my.AVG.cond2, my.AVG.cond)
head(myData.filter)

#make simple line or bar graph
#first reorder and clean up the data you want to graph
library(dplyr)
row.names(myData.filter) <- myData.filter[,1]
myData.filter <- select(myData.filter, -geneSymbols)
myData.filter.transpose <- as.data.frame(t(myData.filter))
treatments <- row.names(myData.filter.transpose)
myData.filter.transpose <- mutate(myData.filter.transpose, treatments, 
                                  genotype=c("WT", "het","KO"))
myData.filter.transpose


