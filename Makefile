# Input files
READ = READ.txt.gz

# Directories and parameters
FASTQC = FastQC/fastqc 

# Anything below this point should be changed

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

all: fastqc 

.PHONY: all fastqc
