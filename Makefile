# Input files
READ = READ.txt.gz
GENOME = genome.fasta
GPTT = genome.ptt
GRNT = genome.rnt

# Directories and parameters
FASTQC = FastQC/fastqc 
TRIM5 = 9

# Anything below this point should not be changed

# Directories
SRCDIR = $(CURDIR)/src
EDGEDIR = $(SOFTDIR)/EDGE_pro_v1.3.1 

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

# RPKMs
RPKM = READ.rpkm_0

$(RPKM): $(TREAD) $(GENOME) $(GPTT) $(GRNT)
	$(EDGEDIR)/edge.pl -g $(GENOME) -p $(GPTT) -r $(GRNT) -u $(TREAD) -o READ -s $(EDGEDIR)

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
	samtools view -q 10 -F 0x10 -b -o $(FBAM) -S $(SAM)
$(RBAM): $(SAM)
	samtools view -q 10 -f 0x10 -b -o $(RBAM) -S $(SAM)

FSORT = aln.sorted.f.bam
RSORT = aln.sorted.r.bam
$(FSORT): $(FBAM)
	samtools sort $(FBAM) $(basename $(FSORT))
$(RSORT): $(RBAM)
	samtools sort $(RBAM) $(basename $(RSORT))

FWIG = aln.sorted.f.bigwig
RWIG = aln.sorted.r.bigwig
FBED = aln.sorted.f.bed
RBED = aln.sorted.r.bed

TOWIG = bam_to_wiggle.py
TOBIG = wigToBigWig
TOBED = bigWigToBedGraph

$(TOWIG):
	wget https://raw.githubusercontent.com/chapmanb/bcbio-nextgen/master/scripts/utils/bam_to_wiggle.py && \
	sed -i 's/wigToBigWig/.\/wigToBigWig/g' $(TOWIG)

$(TOBIG):
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
	chmod 755 $(TOBIG)

$(TOBED):
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
	chmod 755 $(TOBED)

$(FWIG): $(FSORT) $(TOWIG) $(TOBIG)
	python2 bam_to_wiggle.py $(FSORT) --normalize
$(RWIG): $(RSORT) $(TOWIG) $(TOBIG)
	python2 bam_to_wiggle.py $(RSORT) --normalize

$(FBED): $(FWIG) $(TOBED)
	./$(TOBED) $(FWIG) $(FBED)
$(RBED): $(RWIG) $(TOBED)
	./$(TOBED) $(RWIG) $(RBED)
counts: $(FBED) $(RBED)

rpkm: $(RPKM)

# Stats
STATS = genome.stats.tsv

$(STATS): $(SAM) $(FSORT) $(RSORT)
	zcat unmapped.gz | $(SRCDIR)/get_stats genome $(FSORT) $(RSORT) > $(STATS) || rm $(STATS)
stats: $(STATS)

all: fastqc trim align counts rpkm stats

.PHONY: all fastqc trim align counts rpkm stats
