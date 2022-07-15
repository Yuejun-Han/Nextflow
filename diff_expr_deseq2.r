# install Bioconductor and DESeq2; only need to run this once
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install('apeglm')

# Load DESeq2
library("DESeq2")

# be sure to set your working directory

# load matrix and meta data
cts <- as.matrix(read.table('hwk_expr_matrix.tsv', header=TRUE, row.names = 'gene_id'))
coldata <- read.table('hwk_expr_meta.tsv', sep='\t', header=TRUE)

# create deseq data set
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ stage)

# Pre-filter low expressing data
keep <- rowSums(counts(dds)) >= 400
dds <- dds[keep,]

# Factor levels set alphabetically by default
# we want to be sure our reference point is level 1 and the rest are in desired order
dds$stage <- factor(dds$stage, levels = c("stage_I", "stage_IV"))

# Run differential expression analysis
dds <- DESeq(dds)

# get results
res <- results(dds)
res

# alternative ways to get the same thing but more explicitly
#res <- results(dds, name="condition_tumor_vs_normal")
#res <- results(dds, contrast=c("condition", "tumor", "normal"))

# summarize results
summary(res)

# how many adjusted pvalues were less than 0.05?
summary(res$padj < 0.05, na.rm=TRUE)

# plotting results
plotMA(res, ylim=c(-2, 2))

# plotting with shrunken log2FC
#resLFC <- lfcShrink(dds, coef="condition_tumor_vs_normal", type="apeglm")
#plotMA(resLFC, ylim=c(-2,2))

# Save results
resOrdered <- res[order(res$padj),]
write.table(as.data.frame(resOrdered), 
          file="stage_stage_I_vs_stage_IV_results.tsv", quote = FALSE)
resOrdered
head(resOrdered, 10)
