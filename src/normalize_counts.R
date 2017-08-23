#!/usr/bin/env Rscript

library("tximport")
library('DESeq2')
library('limma')

args = commandArgs(trailingOnly=TRUE)

# import kallisto's results
s2c <- read.table(file.path(args[1]),
                  header=TRUE,
                  stringsAsFactors=TRUE)
kal_dirs <- sapply(s2c$sample,
                   function(id) file.path(args[2],
					  id,
					  args[3],
					  "abundance.tsv"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c$ibatch = as.numeric(factor(s2c$batch))

txi.kallisto <- tximport(s2c$path,
			 type="kallisto",
			 txOut=TRUE,
			 countsCol='est_counts')
counts <- as.matrix(txi.kallisto$counts)
colnames(counts) <- s2c$sample
for (col in colnames(counts)){
	counts[, col] <- sapply(counts[, col], as.integer)
}
rownames(s2c) <- colnames(counts)

# apply variance stabilizing transformation (vst)
dds <- DESeqDataSetFromMatrix(counts,
			      s2c,
			      ~ 1)
tvst <- vst(dds)
write.table(assay(tvst),
	    file=args[4],
	    sep="\t")

# correct for batch effects
# optionally preserve treatment conditions
design <- model.matrix(~strain, s2c)
if (as.logical(args[6]))
{
  batch <- removeBatchEffect(assay(tvst),
                             batch=s2c$batch,
			     design=design)
} else {
  batch <- removeBatchEffect(assay(tvst),
                             batch=s2c$batch)
}
write.table(batch,
	    file=args[5],
	    sep="\t")
