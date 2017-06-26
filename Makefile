# Input files
READ = READ.txt.gz
GENOME = genome.gff

# Directories and parameters
SOFTDIR =
SRCDIR = $(CURDIR)/src
FASTQC = fastqc 
TRIMDIR = $(SOFTDIR)/Trimmomatic-0.36
TRIMMOMATIC = $(TRIMDIR)/trimmomatic-0.36.jar
ADAPTERS = $(TRIMDIR)/adapters/TruSeq3-SE.fa
MORETRIMMING = 
QCDIR = $(CURDIR)/QC
ORTHOLOGS = 
QDIR = $(CURDIR)/kallisto
KRAKENPARAMS = --fastq-input --gzip-compressed
KRAKENDB = kraken
KCPU = 1
BOOTSTRAPS = 1000
FRAGMENT = 130
FRAGMENTSD  = 70

# Anything below this point should not be changed

$(QCDIR):
	mkdir -p $@
$(QDIR):
	mkdir -p $@

# QC
QCREAD = $(QCDIR)/$(addsuffix _fastqc.zip, $(shell basename $(notdir $(READ)) .txt.gz))

$(QCREAD): $(QCDIR) $(READ)
	$(FASTQC) --outdir $(QCDIR) $(READ)
qc: $(QCREAD)

# Kraken (another QC)
KRAKENOUT = $(QCDIR)/kraken.out
KRAKEN = $(QCDIR)/kraken.txt

$(KRAKENOUT): $(READ) $(QCDIR)
	kraken --db $(KRAKENDB) --threads $(KCPU) $(KRAKENPARAMS) $< > $@

$(KRAKEN): $(KRAKENOUT)
	kraken-translate --db $(KRAKENDB) $< | awk -F'\t' '{print $$2}' | sort | uniq -c > $@
kraken: $(KRAKEN)
	
# Index reference
TRANSCRIPTS = genes.fasta
$(TRANSCRIPTS): $(GENOME)
	$(SRCDIR)/gff2genes $< $(ORTHOLOGS) > $@

INDEX = reference
$(INDEX): $(TRANSCRIPTS)
	kallisto index -i $@ $<
index: $(INDEX)

# Quantification
QUANTS = $(QDIR)/abundance.h5
$(QUANTS): $(INDEX) $(QDIR) $(READ)
	java -jar $(TRIMMOMATIC) SE $(READ) /dev/stdout ILLUMINACLIP:$(ADAPTERS):2:30:10 $(MORETRIMMING) | kallisto quant -b 100 -i $(INDEX) -o $(QDIR) --single -l $(FRAGMENT) -s $(FRAGMENTSD) /dev/stdin
quantify: $(QUANTS)

clean:
	-rm $(INDEX)
	-rm $(TRANSCRIPTS)
	-rm -rf $(QCDIR)
	-rm -rf $(QDIR)

all: qc kraken index quantify 

.PHONY: all qc kraken index quantify clean
