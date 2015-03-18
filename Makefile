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
QCREAD = $(QCDIR)/$(addsuffix _fastqc.zip, $(basename $(notdir $(READ))))

$(QCREAD): $(QCDIR) $(READ)
	$(FASTQC) --outdir $(QCDIR) $(READ)
fastqc: $(QCREAD)

# Trim

# Bowtie-index
GIDX = genome.1.bt2

$(GIDX): $(GENOME)
	bowtie2-build $(GENOME) genome



all: fastqc 

.PHONY: all fastqc
