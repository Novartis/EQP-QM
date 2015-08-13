# Exon quantification pipeline - quantification module (EQP-QM)
EQP-QM is a Unix based RNA-seq quantification module which uses SAM/BAM
genome alignment files as input and creates gene, exon, and junctions counts.


## Installation of EQP-QM

For the installation of EQP-QM just copy the files and directories
obtained from unpacking the GitHub download into a directory that is contained
in your `PATH` environment variable. Please note that the relative directory
structure must be maintained as the shell scripts `eqp-setup.sh` and
`eqp-quantify.sh` assume that auxilliary scripts and Java programs are
located at certain paths relative to their location.


## Dependencies
* Python (>= version 2.6.5, imported libraries: copy, gettext, gzip,
  numpy, os, re, sets, sys, textwrap, warnings)
* Java (>= version 1.6)
* samtools (>= version 0.1.17)
* bedtools (>= version 2.24.0)


## Running EQP-QM

In order to be able to use EQP-QM to quantify genes, exons, and junctions
it is necessary to execute two preparatory steps:

1. First a number of annotation files required by EQP-QM need to be
   created. This is achieved by running the script `eqp-setup.sh` on a GTF
   file which contains the genome annotation for the genome that will be
   used to align the sample Fastq files. `eqp-setup.sh` is called as
   follows:

   > `eqp-setup.sh <GTF file> <data directory>`

   Here `<data directory>` is the directory in which all necessary
   annotation files for EQP-QM are created. It needs to be supplied in
   the quantification step. The setup will take a while to complete but it
   is only necessary to execute this step once for each genome
   annotation.

   **GTF file**:  
   Only entries with feature type `exon` which contain a `gene_id` field
   are used by `eqp-setup.sh` (and, thus, for the quantification of genes,
   exons, and junctions).

   **Caveat**: Note that the standard GTF file provided by UCSC contains the
   transcript id in the `gene_id` field (which equals the
   `transcript_id` field). In this case EQP-QM will generate counts for
   the transcripts; however, please note that these cannot be considered
   as transcript abundance estimates; in particular, one read can contribute
   to many transcripts. Ensembl GTF files work without problems.

2. Secondly the Fastq files need to be aligned against the chosen
   reference genome using a (splice-aware) aligner of your choice (e.g.
   HiSAT, STAR, or Tophat2) to create SAM/BAM files.

   **Caveat**: Please note that the chromosome names for NCBI/UCSC are
   different from the chromosome names used by Ensembl. So it is not
   possible to use a genome downloaded from NCBI/UCSC for alignment and
   an Ensembl GTF file for the genome annotation. In general, please
   make sure that the identifiers in the genome Fasta file and the GTF
   file are consistent.

Note that the above two steps are independent of each other and can be
executed in parallel.

Once the alignment and the setup of EQP-QM are finished, the quantification
step can be invoked on a SAM/BAM file via:

> `eqp-quantify.sh -d <data directory> <output directory> <SAM/BAM file>`

This will create the files `<SAM/BAM file base>`-gene.cnt,
`<SAM/BAM file base>`-gene.cnt, and `<SAM/BAM file base>`-junction.cnt in
the directory `<output directory>` if `<SAM/BAM file base>` is the basename
of `<SAM/BAM file>` without extension `.sam` or `.bam`.

The run time depends on the number of reads in the SAM/BAM file. Expect
~0.5-1h for ~10M paired-end reads. EQP-QM needs at least 10GB of main
memory.


## `eqp-quantify.sh` options

Usage: `eqp-quantify.sh <options> -d <setup dir> <output dir> <SAM/BAM file>`

where `<options>` is  
>  `[-g] [-e] [-j] [-E <exon overlap>] [-J <junction overlap>]
   [-W <min read weight>] [-s <direction>] [-o <output prefix>]
   [--nosort] [--unambig] [--unweighted] [-w <weight file>]`

`output dir`: directory used for the count files
`SAM/BAM file`: the file containing the alignments of the
   reads against the genome with the aligner.  
`-d STRING`: Use STRING as the directory that contains the auxilliary
  files (needs to be supplied)  
`-g`: compute gene counts  
`-e`: compute exon counts  
`-j`: compute junction counts  
`-E INT`: Minimal overlap of a read with an exon [5]  
`-J INT`: Minimal overlap of a read with both exons on a junction [8]  
`-W FLOAT`: Minimal weight of a read; reads with a lower weight are
          disregarded [0.01]  
`-s STRING`: process reads as strand-specific in direction STRING  
`-o STRING`: A directory STRING is created in the current working directory
    and all intermediate files are stored in a directory structure under
    this directory; furthermore, the count files are generated in the current
    working directory with the prefix STRING.  
`--nosort`: the alignment file is already sorted by names; do not sort it.  
`--unambig`: count only reads that can be assigned unambiguously to a single gene
   or exon when creating the gene or exon counts  
`--unweighted`: do not use read weights in the generation of counts  
`-w STRING`: use STRING as weight file instead of computing it from the genomic
   alignments  


