Allele-specific expression and allele-specific open chromatin analysis pipeline
-------------------------------------------------------------------------------

`git clone --recursive https://github.com/tobiasrausch/allis.git`

`cd allis`

`make all`


Prerequisites
-------------

Download the 1000 Genomes Reference Panel

`cd refpanel && ./download_1kGP_hg19.sh`

Build a BED file of all exonic regions

`cd R && Rscript exon.R`


Running Allis
-------------

`./src/allis.sh <rna_seq.bam> <hg19.fasta> <outprefix> <exome_variants.vcf.gz>`


Plotting Results
----------------

`Rscript R/allis.R example/HG00096.tsv.gz R/exon.hg19.bed.gz`
