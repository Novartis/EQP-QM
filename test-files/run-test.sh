#!/bin/bash

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
PROG_DIR=`cd "$PROG_DIR" && pwd`

cd $PROG_DIR


################################################################################
##
##  Run eqp-setup.sh
##
################################################################################

## Create reference file directory
if [ ! -d data-files ]
then
  mkdir data-files
fi

EQP_SETUP_LOCATION=`which eqp-setup.sh  2>&1 | fgrep "which:"`
if [ "$EQP_SETUP_LOCATION" != "" ]
then
  echo "$EQP_SETUP_LOCATION"
  echo "Please make sure that eqp-setup.sh is contained in a directory which is"
  echo "included in your PATH environment variable ... exiting"
  exit 1
fi

eqp-setup.sh gtf-files/Homo_sapiens.GRCh38.76-test.gtf data-files
if [ $? -ne 0 ]
then
  echo "Creation of reference files by calling"
  echo "   eqp-setup.sh test-files/gtf-files/Homo_sapiens.GRCh38.76-test.gtf test-files/data-files"
  echo "failed ... exiting."
  exit 1
fi


################################################################################
##
##  Run eqp-quantify.sh
##
################################################################################

EQP_QUANTIFY_LOCATION=`which eqp-quantify.sh 2>&1 | fgrep "which:"`
if [ "$EQP_QUANTIFY_LOCATION" != "" ]
then
  echo "$EQP_QUANTIFY_LOCATION"
  echo "Please make sure that eqp-quantify.sh is contained in a directory which is"
  echo "included in your PATH environment variable ... exiting"
  exit 1
fi

eqp-quantify.sh -d data-files count-files sam-files/test-alignment-file.bam
if [ $? -ne 0 ]
then
  echo "Quantification of file test-files/sam-files/test-alignment-file.bam by"
  echo "eqp-quantify.sh failed ... exiting."
  exit 1
fi


################################################################################
##
##  Check if the count files are correct
##
################################################################################

echo
COMPARISON_DIR="comparison-files"
for FILE in test-alignment-file-gene test-alignment-file-exon test-alignment-file-junction
do
  echo "Checking file count-files/$FILE.cnt"
  paste count-files/$FILE.cnt $COMPARISON_DIR/$FILE.cnt > count-files/test-alignment-file-gene.cmp
  DIFF=`awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (abs($2 - $4) > 0.00001) { print "Difference detected"; exit 1} }' count-files/test-alignment-file-gene.cmp`
  
  if [ "$DIFF" != "" ]
  then
    echo "The computed values by the test for file $FILE"
    echo "are not correct - please check the output."
    exit 1
  fi
  echo "... correct"
  echo
done

echo "Test successfully completed."

