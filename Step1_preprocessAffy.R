library(oligo)
#the pd.mogene.2.0.st maps probes to probnormDatas during the normalization/summarization step
library(pd.mogene.2.0.st)
#if you don't have the pd.mogene.2.0.st package, you will need to generate it and load from source as follows:
#install.packages("pd.mogene.2.0.st", repos=NULL,type="source")
#the mogene20sttranscriptcluster.db package maps probnormDatas to genes in the annotation step
library(mogene20sttranscriptcluster.db)
#if you don't have the mogene20sttranscriptcluster.db package, you will need to generate it and load from source as follows:
#install.packages("mogene20sttranscriptcluster.db", repos=NULL,type="source")

#take a look at these annotation packages in more detail
pd.mogene.2.0.st
mogene20sttranscriptcluster.db
cols(mogene20sttranscriptcluster.db)
keytypes(mogene20sttranscriptcluster.db)

#read in the raw .CEL files (will take ~5-10min)
rawData <- read.celfiles(list.celfiles())
rawData

#normalize using RMA (will take <5min)
normData <- rma(rawData)
normData
