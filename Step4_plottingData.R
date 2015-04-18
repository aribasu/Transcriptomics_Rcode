#this script walks thorough how to create publication quality graphics from transcriptomic data generated (regardless of platform used)
#to start this script you need a file with all your expression data and some non-redundant identifiers as row names (usually gene symbols)
#you also need a study design file

#for data wrangling, we'll using the popular dplyr package
#for creating graphics, we'll use the ggplot2 and ggvis packages which employ a 'grammar of graphics' approach
#load the packages
library(ggplot2)
library(dplyr)
library(ggvis)

#simple bar graph
ggplot(myData.filter.transpose, aes(y=Il1a)) +
  geom_bar(aes(x=treatments, fill=genotype), stat="identity") +
  theme(axis.text.x=element_text(angle=-45))

#create a basic scatterplot using ggplot
ggplot(myData, aes(x=B6.untreated.AVG, y=B6.LPS_6hr.AVG)) +
  geom_point(shape=1) +
  geom_point(size=4)

##Volcano plots
##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)

##Construct the volcano plot
ggplot(data=myData, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  opts(legend.position = "none") +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")

#define a tooltip that shows gene symbol and Log2 expression data when you mouse over each data point in the plot
tooltip <- function(data, ...) {
  paste0("<b>","Symbol: ", data$geneSymbols, "</b><br>",
         "Ripk3.LPS.6hr: ", data$Ripk3.LPS_6hr.AVG, "<br>",
         "dblKO.LPS.6hr: ", data$Ripk3_Casp8.LPS_6hr.AVG)
}

#plot the interactive graphic
myData %>% 
  ggvis(x= ~Ripk3.LPS_6hr.AVG, y= ~Ripk3_Casp8.LPS_6hr.AVG, key := ~geneSymbols) %>% 
  layer_points(fill = ~LogFC_dblKOvsRipk3_2hr) %>%
  add_tooltip(tooltip)

