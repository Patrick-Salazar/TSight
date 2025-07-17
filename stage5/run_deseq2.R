#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
})

# Paths
count_file  <- "stage4/counts/gene_counts.csv"
coldata     <- read_csv("stage5/config/coldata.csv")
out_dir     <- "stage5/results"
log_file    <- "stage5/logs/deseq2.log"

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(log_file), recursive=TRUE, showWarnings=FALSE)

# Redirect output
sink(log_file, split=TRUE)

# 1. Read counts
cts <- read_csv(count_file) 
rownames(cts) <- cts[[1]]
cts <- cts[,-1] 

# 2. Prepare metadata
rownames(coldata) <- coldata$sample
coldata$sample <- NULL

# 3. Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(cts),
  colData   = coldata,
  design    = ~ condition
)

# 4. Prefilter
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 5. Run DESeq
dds <- DESeq(dds)

# 6. Extract results
res <- results(dds, contrast=c("condition","tumor","normal"))
res <- lfcShrink(dds, coef="condition_tumor_vs_normal", res=res)

# 7. Save
write.csv(as.data.frame(res),
          file.path(out_dir,"deseq2_results.csv"),
          row.names=TRUE)

# 8. MA-plot
png(file.path(out_dir,"MA_plot.png"))
plotMA(res, main="DESeq2 MA-plot", ylim=c(-5,5))
dev.off()

sink()
