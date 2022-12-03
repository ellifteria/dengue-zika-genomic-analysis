library("DESeq2") # load in DESeq 2 library

gene.counts.file = "./gene_counts.csv"
gene.annotations.file = "./gene_annotations copy.csv"

cts <- as.matrix(read.csv(gene.counts.file,
			   row.names = "gene_id"
			  )) # restructure gene count data into an R matrix


coldata <- read.csv(gene.annotations.file, row.names=1) # read in column data from the gene annotation data
coldata <- coldata[, c("donor", "infection_status")] # restructure coldata into donor and infection_status
coldata$donor <- factor(coldata$donor) # factor column data by donor
coldata$infection_status <- factor(coldata$infection_status) # factor column data by infection_status
rownames(coldata) <- sub("fb", "", rownames(coldata)) # assign row names for column data

print(rownames(coldata))

cts <- cts[, rownames(coldata)] # add to gene count matrix

cts[1:5, ] # display first 5 rows of gene count data

coldata # display all column data

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~infection_status) # prepare gene count matrix for DESeq analysis

dds$infection_status <- relevel(dds$infection_status, ref = "uninfected")

dds <- DESeq(dds)  # run DESeq analysis on data matrix
res <- results(dds) # generate a results table for the DESeq analysis performed on the gene count data
res  # display the results table for the DESeq analysis

res <- results(dds, contrast=c("infection_status", "DENV_infected", "uninfected")) # generate a results table for the DESeq analysis with the contrast specified
res # display the results table

resultsNames(dds) # print the coefficients of the results table

resLFC <- lfcShrink(dds, coef="infection_status_DENV_infected_vs_uninfected", type="apeglm") # shrink the effect size for the coefficient `donor_B_vs_A` to facilitate easier visualization and ranking
resLFC # display the results table after effect size shrinkage

resOrdered <- res[order(res$pvalue),] # order the results by smallest p-value
summary(res) # display a summary of the results table

sum(res$padj < 0.1, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.1

res05 <- results(dds, alpha=0.05) # generate a new results table with alpha set to 0.05 instead of the default value of 0.1
summary(res05) # display a summary of the new results table

sum(res05$padj < 0.05, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.05

library("IHW") # load the IHW (independent hypothesis weighting) library

resIHW <- results(dds, filterFun=ihw) # generate a results table filtered by the independent hypothesis weighting function
summary(resIHW) # display a summary of the results table

sum(resIHW$padj < 0.1, na.rm=TRUE) # calculate the number of adjusted p-values less than 0.1

metadata(resIHW)$ihwResult # display the metadata for independent hypothesis weighting results of the results table 

plotMA(res, ylim=c(-2,2)) # plot the log2 fold change for the mean of the normalized counts of all the samples in the gene count DESeq analysis results table

plotMA(resLFC, ylim=c(-2, 2)) # plot the log2 fold change for the mean of the normalized counts of all the samples in the gene count DESeq analysis results table post-shrinkage

# idx <- identify(res$baseMean, res$logFoldChange) # interactively detect the row number of individual genes
rownames(res)[idx] # display the gene identifiers of the selected genes

resultsNames(dds) # display the coefficients of the gene count DESeq analysis

resNorm <- lfcShrink(dds, coef="infection_status_DENV_infected_vs_uninfected", type="normal") # perform shrinkage using the normal method; `coef=2` because we are interested in `donor_B_vs_A`
resAsh <- lfcShrink(dds, coef="infection_status_DENV_infected_vs_uninfected", type="ashr") # perform shrinkage using the adaptive shrinkage estimator: `coef=2` because we are interested in `donor_B_vs_A`

par(mfrow=c(1,3), mar=c(4,4,2,1)) # set up plotting so that it displays plots adjacent to each other
xlim <- c(1,1e5) # set limits on the x-axis to be used across plots
ylim <- c(-3,3) # set limits on the y-axis to be used across plots
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm") # plot post-shrinkage data using apeglm shrinkage estimator
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal") # plot post-shrinkage data using normal shrinkage estimator
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")# plot post-shrinkage data using adaptive shrinkage estimator, ashr

plotCounts(dds, gene=which.min(res$padj), intgroup="infection_status") # plot normalized gene counts by donor

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="infection_status", returnData=TRUE) # get data from plot of normalized gene counts by donor

library("ggplot2") # load ggplot2 library
ggplot(d, aes(x=infection_status, y=count)) + 
	geom_point(position=position_jitter(w=0.1, h=0)) +
	scale_y_log10(breaks=c(25, 100, 400)) # plot data in customized plot

mcols(res)$description # display information on the variables and tests performed

resOrdered # display the ordered results

write.csv(as.data.frame(resOrdered),
	  file="DENV_infected_results.csv") # write the ordered results as a dataframe to a csv file

resSig <- subset(resOrdered, padj < 0.15) # get the subset of the results with an adjusted p-value less than 0.1
resSig # display the significant results

write.csv(as.data.frame(resSig),
	  file="DENV_infected_results_significant.csv") # write the significant ordered results as a data frame to a csv file

colData(dds) # display the column data for the DESeq analysis data

ddsMF <- dds # create a copy of the DESeq analysis for multi-factored design analysis

levels(ddsMF$infection_status) # display the levels of the `infection_status` 

levels(ddsMF$infection_status) <- sub("-.*", "", levels(ddsMF$infection_status)) # format level names so only includes letters
levels(ddsMF$infection_status) # display levels

design(ddsMF) <- formula(~ infection_status + donor) # set up multi-factor design experiment analysis
ddsMF <- DESeq(ddsMF) # re-run DESeq analysis with multi-factor design

resMF <- results(ddsMF) # generate results table for multi-factor design
head(resMF) # display header row from results table

resMFinfection_status <- results(ddsMF,
		     contrast=c("infection_status", "DENV_infected", "uninfected")) # generate results table with contrast specified
head(resMFinfection_status) # display header list from results table with contrast specified

vsd <- vst(dds, blind=FALSE) # transform DESeq analysis dataset using vst
rld <- rlog(dds, blind=FALSE) # transform DESeq analysis dataset using rlg
head(assay(vsd), 3) # display header row of vsd transformed data

ntd <- normTransform(dds) # compute normal transform giving log2(n + 1) of DESeq analysis data
library("vsn") # load bioconductor library for variance stabilization and calibration
meanSdPlot(assay(ntd)) # create mean standard deviation plot of normal tranform data

meanSdPlot(assay(vsd)) # create mean standard deviation plot of vst transform data

meanSdPlot(assay(rld)) # create mean standard deviation plot of rld tranform data

library("pheatmap") # load library for creating clustered heat maps
select <- order(rowMeans(counts(dds, normalized=TRUE)),
		decreasing=TRUE)[1:20] # define selection criteria for data for heat map
df <- as.data.frame(colData(dds)[,c("donor", "infection_status")]) # restructure data as a data frame for heat map
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
	 cluster_cols=FALSE, annotation_col=df) # create heat map of normal transform data

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
	 cluster_cols=FALSE, annotation_col=df) # create heat map of vsd tranform data

pheatmap(assay(rld)[select, ], cluster_rows=FALSE, show_rownames=FALSE,
	 cluster_cols=FALSE, annotation_col=df) # create heat map of rld transform data

sampleDists <- dist(t(assay(vsd))) # save sample distribution from vst tranform data

library("RColorBrewer") # load RColorBrewer library to edit colors for plot
sampleDistMatrix <- as.matrix(sampleDists) # make sample distributions a matrix
rownames(sampleDistMatrix) <- paste(vsd$donor, vsd$infection_status, sep="-") # set row names of sample distribution matrix
colnames(sampleDistMatrix) <- paste(vsd$donor, vsd$infection_status, sep="-") # NULL # set column names of sample distribution matrix to NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Purples")))(255) # define colors desired for plot
pheatmap(sampleDistMatrix,
	 clustering_distance_rows=sampleDists,
	 clustering_distance_cols=sampleDists,
	 col=colors) # plot sample distributions as heat map

plotPCA(vsd, intgroup=c("infection_status", "donor")) # perform principal component analysis of vst transform data and plot

pcaData <- plotPCA(vsd, intgroup=c("infection_status", "donor"), returnData=TRUE) # extract data from principal component analysis plot
percentVar <- round(100 * attr(pcaData, "percentVar")) # scale percent variances up to out of 100% instead of out of 1
ggplot(pcaData, aes(PC1, PC2, color=infection_status, shape=donor)) +
	geom_point(size=3) +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	coord_fixed() # plot principal component analysis plot customized

