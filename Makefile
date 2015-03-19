# Input files
READ = READ.txt.gz
GENOME = genome.fasta

# Directories and parameters
FASTQC = FastQC/fastqc 

# Anything below this point should not be changed

# Directories
SRCDIR = $(CURDIR)/src

QCDIR = $(CURDIR)/QC
$(QCDIR):
	mkdir -p $(QCDIR)

# QC
QCREAD = $(QCDIR)/$(addsuffix _fastqc.zip, $(shell basename $(notdir $(READ)) .txt.gz))

$(QCREAD): $(QCDIR) $(READ)
	$(FASTQC) --outdir $(QCDIR) $(READ)
fastqc: $(QCREAD)

# Trim
TREAD = $(addsuffix .fq.gz, $(shell basename $(notdir $(READ)) .txt.gz))

$(TREAD): $(QCREAD) $(READ)
	cd $(QCDIR) && unzip -o $(addsuffix _fastqc.zip, $(shell basename $(notdir $(READ)) .txt.gz))
	cp $(READ) tmp.txt.gz
	for oligo in $$($(SRCDIR)/fastqc_adapters $(QCDIR)/$(addsuffix _fastqc, $(shell basename $(notdir $(READ)) .txt.gz))/fastqc_data.txt); do \
	  zcat tmp.txt.gz | trim_blast_short -l $$oligo -z > tmp.trimmed.gz; \
	  mv tmp.trimmed.gz tmp.txt.gz; \
	done
	mv tmp.txt.gz $(TREAD)
trim: $(TREAD)

# Bowtie-index
GIDX = genome.1.bt2

$(GIDX): $(GENOME)
	bowtie2-build $(GENOME) genome

# Bowtie-align
SAM = genome.sam.gz

$(SAM): $(TREAD) $(GIDX)
	bowtie2 -x genome -U $(TREAD) -t --trim5 9 --al-gz mapped.gz --un-gz unmapped.gz | gzip --stdout > $(SAM)
align: $(SAM)

all: fastqc trim align

.PHONY: all fastqc trim align
