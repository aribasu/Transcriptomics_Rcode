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
targets <- readTargets("Igor_Ripk3_BMMC_studyDesign.txt", row.names=NULL)
targets
genotype <- factor(targets$description)
#create some more human-readable labels for your samples using the info in this file
sampleLabels <- paste(targets$description, targets$rep, sep=".")

#set-up your experimental design
#design <- model.matrix(~genotype)
design <- model.matrix(~0+genotype)
colnames(design) <- levels(genotype)
design

##############################################################################################################################
#check read quality using shortRead package
##############################################################################################################################
myFastq <- targets$fastq
#collecting statistics over the files
qaSummary <- qa("B6-WT-untreat-rep2_S2_mergedLanes_read1.fastq", type="fastq")
#creat and view a report
browseURL(report(qaSummary))


##############################################################################################################################
#build index from your reference genome (expect this to take about 20 min on 8G RAM for mouse genome)
#you must have already downloaded the fasta file for your genome of interest and have it in the working directory
##############################################################################################################################
buildindex(basename="mouse",reference="Mus_musculus.GRCm38.dna.primary_assembly.fa", memory=64000)

##############################################################################################################################
#align your reads (in the fastq files) to your indexed reference genome that you created above
#expect this to take about 45min for a single fastq file containing 25 million reads
#the output from this is a .bam file for each of your original fastq files
##############################################################################################################################
reads1 <- targets$fastq[2] 
reads2 <- targets$fastq[14] 
align(index="mouse", readfile1=reads1, readfile2=reads2, input_format="gzFASTQ",output_format="BAM",
      output_file="alignmentResultsPE.BAM", tieBreakHamming=TRUE,unique=TRUE,indels=5. nthreads = 8)

##############################################################################################################################
#now that your reads are aligned to the reference genome
#use the 'featureCounts' function to summarize read counts to genomic features (exons, genes, etc)
#summarizing reads to features will take less than 1min per .bam file.
#for canine GTF file, I get ~50% of reads summarizing to features
##############################################################################################################################
fc <- featureCounts(files=targets$OutputFile, annot.ext="Canis_familiaris.CanFam3.1.77.gtf", isGTFAnnotationFile=TRUE, 
                    useMetaFeatures=TRUE, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, strandSpecific=1)
x <- DGEList(counts=fc$counts, genes=fc$annotation)
#if you need RPKM, they can generated as follows
x_rpkm <- rpkm(x,x$genes$Length)
head(x_rpkm)

#Filtering. 
#Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
isexpr <- rowSums(cpm(x) > 10) >= 2
x <- x[isexpr,]

#Normalization. 
#Perform voom normalization:
normData <- voom(x, design, plot=TRUE)
exprs <- normData$E
exprs.matrix <- as.matrix(exprs)
head(exprs.matrix)

#explore your data using some standard graphs
#choose color scheme for graphs
cols.ALL <- topo.colors (n=9, alpha=1)
hist(exprs, xlab = "log2 expression", main = "normalized data - histograms", col=cols.ALL)
boxplot(exprs, ylab = "normalized log2 expression", main = "non-normalized data - boxplots", col=cols.ALL)

