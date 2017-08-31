#!/usr/bin/env Rscript

library("tximport")
library('DESeq2')
library('BiocParallel')

args = commandArgs(trailingOnly=TRUE)

# rna samples
# directory
# kallisto directory
# p-val threshold
# log(FC) threshold (absolute)
# number of cores

# prepare samples table
s2c <- read.table(file.path(args[1]),
                  header=TRUE,
                  stringsAsFactors=TRUE)
kal_dirs <- sapply(s2c$sample,
                   function(id) file.path(args[2],
					  id,
					  args[3],
					  "abundance.tsv"))
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# import kallisto's counts
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

# run DESeq2
dds <- DESeqDataSetFromMatrix(counts,
			      s2c,
	     	              ~mstrain)
multicoreParam <- MulticoreParam(workers=as.numeric(args[6]))
dds <- DESeq(dds,
	     parallel=TRUE,
	     BPPARAM=multicoreParam)

# overall results
res <- results(dds,
	       alpha=as.numeric(args[4]),
	       lfcThreshold=as.numeric(args[5]),
	       altHypothesis='greaterAbs',
	       parallel=TRUE,
	       BPPARAM=multicoreParam)
write.csv(as.data.frame(res),
	  file='overall.tsv',
          sep='')
# strain's contrasts
for(i in 2:length(unique(s2c$mstrain)))
{
    res <- results(dds,
                   contrast=c('mstrain',
                              toString(unique(s2c$mstrain)[i]),
                              'reference'),
		   alpha=as.numeric(args[4]),
		   lfcThreshold=as.numeric(args[5]),
		   altHypothesis='greaterAbs',
		   parallel=TRUE,
		   BPPARAM=multicoreParam)
    write.csv(as.data.frame(res),
              file=paste(toString(unique(s2c$mstrain)[i]),
                         '.tsv',
              sep=''))    
}
