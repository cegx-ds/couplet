# couplet

**couplet** will resolve paired-end reads generated from **CEGX 5-Letter Seq** sequencing kits into 4-letter genomic reads annotated with the epigenomic status of cytosines.

It takes paired-end FASTQ files as an input and generates an output FASTQ file of resolved reads.

Additionally, metrics associated with the resolution process are written to a yaml file and reads that fail to resolve are written to a discarded reads FASTQ file-pair.
