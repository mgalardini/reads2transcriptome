reads2transcriptome
===================

From second generation sequencing reads to nucleotide counts
This small pipeline works with single illumina reads (no PE), and will output the counts separated by strand. 

Usage
-----

* `make all`

Notes
-----

* You may want to run `make fastqc` first and use its results to tweak the trimming parameters
    * By default any overrepresented sequence (e.g. adapters) will be trimmed

Prerequisites
-------------

* Python + Biopython + bcbio-nextgen + pysam
* FastQC
* seq_crumbs
* bowtie2
* samtools

Copyright
---------

Copyright (C) <2015> EMBL-European Bioinformatics Institute

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
GNU General Public License for more details.

Neither the institution name nor the name reads2transcriptome
can be used to endorse or promote products derived from
this software without prior written permission.
For written permission, please contact <marco@ebi.ac.uk>.

Products derived from this software may not be called reads2transcriptome
nor may reads2transcriptome appear in their names without prior written
permission of the developers. You should have received a copy
of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.
