##############################################################################################################################
#This script carries out the steps involved in analysis of RNAseq data.  
#depending on your data and interests, different parts of this script may not apply to you
##############################################################################################################################
#begin by loading the packages required for RNAseq data
library(Rsubread)
library(limma)
library(edgeR)
library(ShortRead)
options(digits=2)
#library(Biostrings)

#read in your study design file
targets <- readTargets("Igor_Ripk3_BMMC_studyDesign_RNAseq.txt", row.names=NULL)
targets
groups <- factor(paste(targets$genotype, targets$treatment, sep="."))
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$genotype, targets$treatment, targets$rep, sep=".")

#set-up your experimental design
design <- model.matrix(~0+groups)
colnames(design) <- levels(groups)
design

# ##############################################################################################################################
# #you can check read quality using shortRead package
# #but I usually find it is better to do this on the sequencer or Illumina's BaseSpace website
# ##############################################################################################################################
# myFastq <- targets$fastq
# #collecting statistics over the files
# qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
# #create and view a report
# browseURL(report(qaSummary))

##############################################################################################################################
#build index from your reference genome (expect this to take about 20 min on 8G RAM for mouse genome)
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
#this only needs to be done once, then index can be reused for future alignments
##############################################################################################################################
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa")

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads1 <- targets$fastq[9]
reads2 <- targets$fastq[21] 
align(index="mouse", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE_sample9.BAM", tieBreakHamming=TRUE,unique=TRUE,indels=5, nthreads=8)

##############################################################################################################################
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#will take about 1-2min per .bam file.
#for total transcriptome data summarized to mouse/human .gtf, expect about 50-60% of reads summarize to genes (rest is non-coding)
##############################################################################################################################
#read in text file with bam file names
MyBAM <- read.delim("BAMfiles_names.txt", header=T)
MyBAM <- as.character(MyBAM[,1])
#summarize aligned reads to genomic features (i.e. exons)
fc <- featureCounts(files=MyBAM, annot.ext="Mus_musculus.GRCm38.79.gtf", isGTFAnnotationFile=TRUE, GTF.featureType = "exon",
                    GTF.attrType="gene_id", useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=2, nthreads=8)
#use the 'DGEList' function from EdgeR to make a 'digital gene expression list' object
DGEList <- DGEList(counts=fc$counts, genes=fc$annotation)
load("DGEList")
#retrieve all your gene/transcript identifiers from this DGEList object
myEnsemblIDs <- DGEList$genes$GeneID
dim(fc$counts)
tail(fc$counts)

##############################################################################################################################
#Normalize unfiltered data using 'voom' function in Limma package
#This will normalize based on the mean-variance relationship
#will also generate the log2 of counts per million based on the size of each library (also a form of normalization)
##############################################################################################################################
normData.unfiltered <- voom(DGEList, design, plot=TRUE)
exprs.unfiltered <- normData.unfiltered$E
exprs.matrix.unfiltered <- as.matrix(exprs.unfiltered)
#note that because you're now working with Log2 CPM, much of your data will be negative number (log2 of number smaller than 1 is negative) 
head(exprs.matrix.unfiltered)

#if you need RPKM for your unfiltered, they can generated as follows
#Although RPKM are commonly used, not really necessary since you don't care to compare two different genes within a sample
rpkm.unfiltered <- rpkm(DGEList, DGEList$genes$Length)
rpkm.unfiltered <- log2(rpkm.unfiltered + 1)

##############################################################################################################################
#Filtering your dataset
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
##############################################################################################################################
cpm.matrix.filtered <- rowSums(cpm(DGEList) > 10) >= 2
DGEList.filtered <- DGEList[cpm.matrix.filtered,]
dim(DGEList.filtered)
rpkm.filtered <- rpkm(DGEList.filtered, DGEList.filtered$genes$Length) #if you prefer, can use 'cpm' instead of 'rpkm' here
rpkm.filtered <- log2(rpkm.filtered + 1)

#Use Voom again to normalize this filtered dataset
normData <- voom(DGEList.filtered, design, plot=TRUE)
exprs <- normData$E
exprs.matrix.filtered <- as.matrix(exprs)
head(exprs.matrix.filtered)


##############################################################################################################################
#annotate your normalized data using the organism-specific database package
##############################################################################################################################
library(org.Mm.eg.db)
library(AnnotationDbi)
#If we want to know what kinds of data are retriveable via the 'select' command, look at the columns of the annotation database
columns(org.Mm.eg.db)
#If we want to know what kinds of fields we could potentially use as keys to query the database, use the 'keytypes' command
keytypes(org.Mm.eg.db)
#transform you identifiers to entrezIDs
myAnnot.unfiltered <- select(org.Mm.eg.db, keys=rownames(exprs.matrix.unfiltered), keytype="ENSEMBL", columns=c("ENTREZID", "REFSEQ", "UNIGENE", "GENENAME", "SYMBOL"))
myAnnot.filtered <- select(org.Mm.eg.db, keys=rownames(exprs.matrix.filtered), keytype="ENSEMBL", columns=c("ENTREZID", "REFSEQ", "UNIGENE", "GENENAME", "SYMBOL"))
resultTable.unfiltered <- merge(myAnnot.unfiltered, exprs.matrix.unfiltered, by.x="ENSEMBL", by.y=0)
resultTable.filtered <- merge(myAnnot.filtered, exprs.matrix.filtered, by.x="ENSEMBL", by.y=0)
head(resultTable.unfiltered)
#add more appropriate sample names as column headers
colnames(resultTable.unfiltered) <- c("Ensembl", "entrez", "refseq", "Unigene", "name", "symbol", sampleLabels)
colnames(resultTable.filtered) <- c("Ensembl", "entrez", "refseq", "Unigene", "name", "symbol", sampleLabels)
#now write these annotated datasets out
write.table(resultTable.unfiltered, "normalizedUnfiltered.txt", sep="\t", quote=FALSE)
write.table(resultTable.filtered, "normalizedFiltered.txt", sep="\t", quote=FALSE)
head(resultTable.unfiltered)

#########################################################################
#If there isn't an R database package for your organism,
#you can still annotate your data using the biomaRt package
#########################################################################
library(biomaRt)
#list all the available 'marts' to choose from
listMarts()
#select your mart of interest
vectorBase <- useMart("vb_gene_mart_1502")
#list the databases available in that mart
listDatasets(vectorBase)
#select the database you want to work with
aedesDB <- useDataset("aaegypti_eg_gene", mart=vectorBase)
#list the various filters, or keys, that can used to access this database
listFilters(aedesDB)
#list all the attributes that can be retrieved from the database
listAttributes(aedesDB)
#pull out your rownames of your dataset to use as one of the filters
myEnsemblIDs.filtered <- rownames(rpkm.filtered)
myEnsemblIDs.unfiltered <- rownames(rpkm.unfiltered)
#attach annotation information to your file
myAnnot.filtered <- getBM(attributes = c("external_gene_id", "ensembl_gene_id", "description"),
                          filters = "ensembl_gene_id",   #the kind of data you'll use as keys
                          values = myEnsemblIDs.filtered,   #the exact data you'll use as keys
                          mart=aedesDB)
myAnnot.unfiltered <- getBM(attributes = c("external_gene_id", "ensembl_gene_id", "description"),
                            filters = "ensembl_gene_id",
                            values = myEnsemblIDs.unfiltered,
                            mart=aedesDB)
#merge the annotation information with the expression data
resultTable.filtered <- merge(myAnnot.filtered, rpkm.filtered, by.x="ensembl_gene_id", by.y=0)
resultTable.unfiltered <- merge(myAnnot.unfiltered, rpkm.unfiltered, by.x="ensembl_gene_id", by.y=0)

# ##############################################################################################
# #additional filtering based on annotation
# #can skip this, since I like to do this right before differential gene expression analysis
# ##############################################################################################
# #our initial filtering above based on cpm, only removes relatively low expressed genes
# #now that you have annotation info, you can further filter to reduce to single line of data per gene symbo using one of two options
# #first, check out how many unqiue genes are represented in your data
# dupFiltered <- unique(resultTable.filtered$symbol)
# #use collapseRows function from WGCNA package to collapse your dataset
# library(WGCNA)
# #pull your rownames and unique identifiers
# myIDs <- rownames(resultTable.filtered)
# #retrieve your gene symbols from your data
# mySymbols <- resultTable.filtered[,5]
# #remove all annotation columns so you're left with only numeric data
# resultTable <- resultTable.filtered[,-1:-6]
# myCollapsed <- collapseRows(resultTable, mySymbols, myIDs, method = "MaxMean")
# myCollapsed <- myCollapsed$datETcollapsed
# #now that the matrix is collapsed to give non-redudant list of genes,
# #you could set the symbols to be the row names and move on


