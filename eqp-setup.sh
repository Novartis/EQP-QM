#!/bin/sh

## Copyright 2015 Novartis Institutes for BioMedical Research
## Inc.Licensed under the Apache License, Version 2.0 (the "License"); you
## may not use this file except in compliance with the License. You may
## obtain a copy of the License at
##
## http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing,
## software distributed under the License is distributed on an "AS IS"
## BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
## implied. See the License for the specific language governing
## permissions and limitations under the License.

#$ -cwd
#$ -j y
#$ -S /bin/sh
#$ -V

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
PROG_DIR=`cd "$PROG_DIR" && pwd`
VERSION=1.0.0

## Ensure that pipes report the exit status of the first failed job
set -o pipefail

################################################################################
##
##  Read options
##
################################################################################

while [ "$1" = "-h" ]
do
  if [ "$1" = "-h" ]
  then
    shift
    PRINT_HELP="TRUE"
  fi
done


################################################################################
##
##  Print help
##
################################################################################

if [ "$PRINT_HELP" = "TRUE" -o "$2" = "" ]
then
  echo "Usage: $PROG_NAME <GTF file> <data directory>"
  echo
  exit 0
fi


################################################################################
##
##  Read arguments
##
################################################################################

INPUT_GTF_FILE_PATH=$1
shift

if [ ! -f "$INPUT_GTF_FILE_PATH" ]
then
  echo "GTF file $INPUT_GTF_FILE_PATH not found ... exiting."
  exit 1
fi

INPUT_GTF_FILE=`basename $INPUT_GTF_FILE_PATH`
INPUT_GTF_FILE_BASE=`echo $INPUT_GTF_FILE | sed -e 's/.gz$//' | sed -e 's/.gtf$//'`
INPUT_GTF_FILE_EXT=`echo $INPUT_GTF_FILE | sed -e "s/$INPUT_GTF_FILE_BASE[.]//"`

if [ "$INPUT_GTF_FILE_EXT" = "gtf.gz" ]
then
  echo "Decompressing file $INPUT_GTF_FILE_PATH"
  gunzip $INPUT_GTF_FILE_PATH
  INPUT_GTF_FILE_PATH=`echo $INPUT_GTF_FILE_PATH | sed -e 's/[.]gz$//'`
  INPUT_GTF_FILE=`basename $INPUT_GTF_FILE_PATH`
elif [ "$INPUT_GTF_FILE_EXT" != "gtf" ]
then
  echo "Unknown extension for file $INPUT_GTF_FILE ... exiting."
  exit 1
fi

PROJECT_DATA_DIR=$1
shift

if [ ! -d $PROJECT_DATA_DIR ]
then
  echo "Creating directory $PROJECT_DATA_DIR"
  mkdir -p $PROJECT_DATA_DIR
fi

if [ ! -d $PROJECT_DATA_DIR/bin ]
then
  mkdir $PROJECT_DATA_DIR/bin
fi

GENE_MODEL_PREFIX="$INPUT_GTF_FILE_BASE"

echo "export FILE_BASE=$INPUT_GTF_FILE_BASE"  > $PROJECT_DATA_DIR/bin/setup.sh
echo "export READ_LENGTH=0"                  >> $PROJECT_DATA_DIR/bin/setup.sh
echo "export RADIUS=0"                       >> $PROJECT_DATA_DIR/bin/setup.sh
echo "export PAIRED_END=TRUE"                >> $PROJECT_DATA_DIR/bin/setup.sh
echo "export STRAND_SPECIFIC=FALSE"          >> $PROJECT_DATA_DIR/bin/setup.sh
echo "export STRAND_SPECIFIC_DIRECTION=none" >> $PROJECT_DATA_DIR/bin/setup.sh

chmod a+x $PROJECT_DATA_DIR/bin/setup.sh

R_EXE=`which R`
if [ ! -f "$R_EXE" ]
then
  echo
  echo 'R not found. Please make sure that R is available ... exiting'
  exit 1
fi

PYTHON_EXE=`which python`
if [ ! -f "$PYTHON_EXE" ]
then
  echo
  echo 'python not found. Please make sure that python is available ... exiting'
  exit 1
fi

################################################################################
##
##  Set directories
##
################################################################################

PROJECT_DIR=`dirname $PROG_DIR`

if [ ! -d $PROJECT_DATA_DIR ]
then
  echo "Creating directory $PROJECT_DATA_DIR"
  mkdir $PROJECT_DATA_DIR
fi

PROJECT_GTF_DIR=$PROJECT_DATA_DIR/gtf-files
if [ ! -d $PROJECT_GTF_DIR ]
then
  echo "Directory $PROJECT_GTF_DIR does not exist ... creating"
  mkdir $PROJECT_GTF_DIR
fi

PROJECT_BED_DIR=$PROJECT_DATA_DIR/bed-files
if [ ! -d $PROJECT_BED_DIR ]
then
  echo "Directory $PROJECT_BED_DIR does not exist ... creating"
  mkdir $PROJECT_BED_DIR
fi

PROJECT_MAP_DIR=$PROJECT_DATA_DIR/map-files
if [ ! -d $PROJECT_MAP_DIR ]
then
  echo "Directory $PROJECT_MAP_DIR does not exist ... creating"
  mkdir $PROJECT_MAP_DIR
fi


################################################################################
##
##
##  Create the reference files
##
##
################################################################################

SCRIPTS_DIR=$PROG_DIR/util-scripts

################################################################################
##
## Add exon_id field to GTF file
##
################################################################################

echo "Adding exon ids to GTF file"
GTF_FILE=$PROJECT_GTF_DIR/$GENE_MODEL_PREFIX.gtf
$SCRIPTS_DIR/addExonIdGtf.py -g $INPUT_GTF_FILE_PATH -o $GTF_FILE
if [ $? -ne 0 ]
then
  echo "Problem with addExonIdGtf.py -g $INPUT_GTF_FILE_PATH -o $GTF_FILE ... exiting"
  exit 1
fi


################################################################################
##
## Compute chromosome names
##
################################################################################

cut -f 1 $GTF_FILE | sed -e 's/^/#/' | sed -e 's/$/#/' | sort -u > $PROJECT_GTF_DIR/$GENE_MODEL_PREFIX-chromosomes.txt
if [ $? -ne 0 ]
then
  echo "Problem with cut -f 1 $GTF_FILE | sort -u > $PROJECT_GTF_DIR/$GENE_MODEL_PREFIX-chromosomes.txt ... exiting."
  exit 1
fi


################################################################################
##
## Create file with mapping of junction ids to exons (junction ids are realized as the combination of two exons)
##
################################################################################

echo "Creating junction - exon mapping file - this can take a while"
export GZIPPED_JUNCTION_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_junction_exon.map.gz

if [ ! -f $GZIPPED_JUNCTION_EXON_MAP_FILE ]
then
  JUNCTION_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_junction_exon.map
  if [ ! -f $JUNCTION_EXON_MAP_FILE ]
  then
    export EXON_PIPELINE_DIR=$SCRIPTS_DIR
    export GTF_FILE
    export JUNCTION_EXON_MAP_FILE
    $R_EXE --vanilla --slave < $SCRIPTS_DIR/R/create-junction-exon-map-file.R
    if [ $? -ne 0 ]
    then
      if [ -f $JUNCTION_EXON_MAP_FILE ]
      then
        rm $JUNCTION_EXON_MAP_FILE
      fi
      echo "Problem with create-junction-exon-map-file.R ... exiting"
      exit 1
    fi
  fi

  echo "Compressing file junction - exon mapping file $JUNCTION_EXON_MAP_FILE"
  gzip -f $JUNCTION_EXON_MAP_FILE
fi


################################################################################
##
## Create the files with the mappings of GTF exon ids to EQP exon ids, the file with the mappings
## of EQP exon ids to gene ids, and the BED file with the genomic coordinates of the exons.
##
################################################################################

EXON_EXON_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_exon.map
EXON_GENE_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_gene.map
EXON_BED_FILE=$PROJECT_BED_DIR/${GENE_MODEL_PREFIX}_genome_exons.bed
echo "Creating exon - exon map file $EXON_EXON_MAP_FILE"
$SCRIPTS_DIR/createExonExonMapFile.py -g $EXON_GENE_MAP_FILE -b $EXON_BED_FILE $GTF_FILE $EXON_EXON_MAP_FILE
if [ $? -ne 0 ]
then
  echo "Problem with createExonExonMapFile.py -g $EXON_GENE_MAP_FILE -b $EXON_BED_FILE $GTF_FILE"
  echo " $EXON_EXON_MAP_FILE ... exiting"
  exit 1
fi


################################################################################
##
## Create file with mapping of EQP exon ids to junctions
##
################################################################################

GZIPPED_EXON_JUNCTION_MAP_FILE=$PROJECT_MAP_DIR/${GENE_MODEL_PREFIX}_exon_junction.map.gz
echo "Creating exon - junction map file $EXON_JUNCTION_MAP_FILE"
## The option -w switches on warning messages
$SCRIPTS_DIR/createExonJunctionMapFile.py -w -A -e $EXON_EXON_MAP_FILE -j $GZIPPED_JUNCTION_EXON_MAP_FILE -o $GZIPPED_EXON_JUNCTION_MAP_FILE
if [ $? -ne 0 ]
then
  echo "Problem with createExonJunctionMapFile.py -w -e $EXON_EXON_MAP_FILE -j $GZIPPED_JUNCTION_EXON_MAP_FILE -o $GZIPPED_EXON_JUNCTION_MAP_FILE"
  echo "   ... exiting"
  exit 1
fi

echo "EQP setup successfully completed."
exit 0
