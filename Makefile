# Input files
READ = READ.txt.gz
GENOME = genome.fasta

# Directories and parameters
FASTQC = FastQC/fastqc 
TRIM5 = 9

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
	bowtie2 -x genome -U $(TREAD) -t --trim5 $(TRIM5) --al-gz mapped.gz --un-gz unmapped.gz | gzip --stdout > $(SAM)
align: $(SAM)

# WIG files
# (counts)
FBAM = aln.f.bam
RBAM = aln.r.bam

$(FBAM): $(SAM)
	samtools view -F 0x10 -b -o $(FBAM) -S $(SAM)
$(RBAM): $(SAM)
	samtools view -f 0x10 -b -o $(RBAM) -S $(SAM)

FSORT = aln.sorted.f.bam
RSORT = aln.sorted.r.bam
$(FSORT): $(FBAM)
	samtools sort $(FBAM) $(basename $(FSORT))
$(RSORT): $(RBAM)
	samtools sort $(RBAM) $(basename $(RSORT))

FWIG = aln.sorted.f.bigwig
RWIG = aln.sorted.r.bigwig

TOWIG = bam_to_wiggle.py
TOBIG = wigToBigWig

$(TOWIG):
	wget https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/scripts/utils/bam_to_wiggle.py && \
	sed -i 's/wigToBigWig/.\/wigToBigWig/g' $(TOWIG)

$(TOBIG):
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
	chmod 755 $(TOBIG)

$(FWIG): $(FSORT) $(TOWIG) $(TOBIG)
	python2 bam_to_wiggle.py $(FSORT) --normalize
$(RWIG): $(RSORT) $(TOWIG) $(TOBIG)
	python2 bam_to_wiggle.py $(RSORT) --normalize
counts: $(FWIG) $(RWIG)

all: fastqc trim align counts

.PHONY: all fastqc trim align counts
