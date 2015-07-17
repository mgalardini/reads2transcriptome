#library('edgeR')
library('DESeq2')

# Read input data (counts)
df <- read.table('deseqFile', sep='\t', header=TRUE, row.names=1)
# Convert to matrix, with gene names
countData <- as.matrix(df)

#TODO: filter empty rows

# Read experimental setup
con <- read.table('conditions.tsv', sep='\t', header=TRUE, row.names=1)

# Define the analysis (strain effect > condition)
dds <- DESeqDataSetFromMatrix(countData, con, formula(~ condition+strain))

# Run the analysis
dds <- DESeq(dds)

# Output each strain pairwise comparisons
for(i in 2:length(unique(con$strain)))
{
    res <- results(dds,
                   contrast=c('strain',
                   toString(unique(con$strain)[i]),
                   'Control'))
    resOrdered <- res[order(res$padj),]
    write.csv(as.data.frame(resOrdered),
          file=paste(toString(unique(con$strain)[i]),
                     '.tsv',
                     sep=''))    
}
