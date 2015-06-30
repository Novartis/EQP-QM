1. Installation of EQP-QM

For the installation of EQP-QM you just need to copy the files and directories
in this directory into a directory that is contained in your PATH
environment variable (or add this directory to your PATH). Please note
that the relative directory structure of this directory must be
maintained as the shell scripts eqp-setup.sh and eqp-quantify.sh assume
that auxilliary scripts and java programs are located at certain paths
relative to their location.

Caveat: The subdirectory "tools" contains two executables (bedtools and
samtools) which may not run on your architecture. Should this be the
case please download the bedtools and samtools packages and copy the
executables into the "tools" subdirectory. bedtools should be version
2.2.24 or higher and samtools version 0.1.17 or higher.


2. Running EQP-QM

In order to be able to use EQP-QM to quantify genes, exons, and junctions
it is necessary to execute two preparation steps:

 a. First the Fastq files need to be aligned against the reference
    genome using a (splice-aware) aligner of your choice (e.g. HiSAT,
    STAR, or Tophat2) to create SAM/BAM files
 b. Secondly, a number of annotation files needed by EQP-QM need to be
    created. This is achieved by running eqp-setup.sh on a GTF file
    which contains the genome annotation for the genome that was used
    for the alignment of the Fastq files in the following way:

    eqp-setup.sh <GTF file> <data directory>

    Here <data directory> is the directory in which all necessary
    annotation files for EQP-QM are created. It needs to be supplied in
    the quantification step. The setup may take a couple of minutes but it
    is only necessary to execute this step once for each genome
    annotation.

    GTF file
    The GTF file needs to contain a "gene_id" field. Note that the
    standard GTF file provided by UCSC contains the transcript id as the
    gene_id. EQP-QM will generate counts for the transcripts in this
    case, however, these cannot be considered as transcript abundance
    estimates; in particular, one read can contribute to many
    transcripts. Ensembl GTF files work without problems.

    Caveat: Please note that the chromosome names for NCBI/UCSC are
    different from the chromosome names used by Ensembl. So it is not
    possible to use a genome downloaded from NCBI/UCSC for alignment and
    an Ensembl GTF to provide the genome annotation.

    
Once eqp-setup.sh finishes, EQP-QM can be invoked on a SAM/BAM file via:

   eqp-quantify.sh -d <data directory> <output directory> <SAM/BAM file>

The run time depends on the number of reads in the SAM/BAM file. Expect
~30min for ~10M paired-end reads. EQP-QM needs at least 10GB main
memory.


3. eqp-quantify.sh options

Usage: eqp-quantify.sh <options> -d <setup dir> <output dir> <SAM/BAM file>

where <options> is
   [-g] [-e] [-j] [-E <exon overlap>] [-J <junction overlap>]
   [-W <min read weight>] [-s <direction>] [-o <output prefix>]
   [--nosort] [--unambig] [--unweighted] [-w <weight file>]

output dir: directory used for the count files
BAM file: the file containing the alignments of the
     reads against the genome with the aligner.
 -d STRING: Use STRING as the directory that contains the auxilliary
   files (needs to be supplied)
 -g: compute gene counts
 -e: compute exon counts
 -j: compute junction counts
 -E INT: Minimal overlap of a read with an exon [5]
 -J INT: Minimal overlap of a read with both exons on a junction [8]
 -W FLOAT: Minimal weight of a read; reads with a lower weight are
           disregarded [0.01]
 -s STRING: process reads as strand-specific in direction STRING
 -o STRING: A directory STRING is created in the current working directory
     and all intermediate files are stored in a directory structure under
     this directory; furthermore, the count files are generated in the current
     working directory with the prefix STRING.
 --nosort: the alignment file is already sorted by names; do not sort it.
 --unambig: count only reads that can be assigned unambiguously to a single gene
    or exon when creating the gene or exon counts
 --unweighted: do not use read weights in the generation of counts
 -w STRING: use STRING as weight file instead of computing it from the genomic
    alignments







