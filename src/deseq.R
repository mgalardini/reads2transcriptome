#!/usr/bin/env Rscript

library("tximport")
library('DESeq2')

args = commandArgs(trailingOnly=TRUE)

# rna samples
# directory
# kallisto directory

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

dds <- DESeqDataSetFromMatrix(counts,
			      s2c,
	     	              ~mstrain)
dds <- DESeq(dds)

for(i in 2:length(unique(s2c$mstrain)))
{
    res <- results(dds,
                   contrast=c('mstrain',
                   toString(unique(s2c$mstrain)[i]),
                   'reference'))
    write.csv(as.data.frame(res),
              file=paste(toString(unique(s2c$mstrain)[i]),
                         '.tsv',
              sep=''))    
} 
