# Input files
READ = READ.txt.gz
GENOME = genome.gbk

# Directories and parameters
SOFTDIR =
SRCDIR = $(CURDIR)/src
FASTQC = $(SOFTDIR)/FastQC/fastqc 
TRIMDIR = $(SOFTDIR)/Trimmomatic-0.36
TRIMMOMATIC = $(TRIMDIR)/trimmomatic-0.36.jar
ADAPTERS = $(TRIMDIR)/adapters/TruSeq3-SE.fa
QCDIR = $(CURDIR)/QC
QDIR = $(CURDIR)/kallisto
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

# Index reference
TRANSCRIPTS = genes.fasta
$(TRANSCRIPTS): $(GENOME)
	$(SRCDIR)/gbk2genes $< > $@

INDEX = reference
$(INDEX): $(TRANSCRIPTS)
	kallisto index -i $@ $<
index: $(INDEX)

# Quantification
QUANTS = $(QDIR)/abundance.h5
$(QUANTS): $(INDEX) $(QDIR)
	java -jar $(TRIMMOMATIC) SE $(READ) /dev/stdout ILLUMINACLIP:$(ADAPTERS):2:30:10 | kallisto quant -b 100 -i reference -o $(QDIR) --single -l $(FRAGMENT) -s $(FRAGMENTSD) /dev/stdin
quantify: $(QUANTS)

clean:
	-rm $(INDEX)
	-rm $(TRANSCRIPTS)
	-rm -rf $(QCDIR)
	-rm -rf $(QDIR)

all: qc index quantify 

.PHONY: all qc index quantify clean
