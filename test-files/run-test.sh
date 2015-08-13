#!/bin/sh

PROG_NAME=`basename $0`
PROG_DIR=`dirname $0`
PROG_DIR=`cd "$PROG_DIR" && pwd`

PROJECT_DIR=`dirname $PROG_DIR`

export PATH="${PROJECT_DIR}:$PATH"

cd $PROG_DIR

## Create reference files
if [ ! -d data-files ]
then
  mkdir data-files
fi
eqp-setup.sh gtf-files/Homo_sapiens.GRCh38.76-test.gtf data-files


## Create count files
eqp-quantify.sh -d data-files count-files sam-files/test-alignment-file.bam

echo
for FILE in test-alignment-file-gene test-alignment-file-exon test-alignment-file-junction
do
  echo "Checking file count-files/$FILE.cnt"
  paste count-files/$FILE.cnt comparison-files/$FILE.cnt > count-files/test-alignment-file-gene.cmp
  DIFF=`awk 'function abs(x){return ((x < 0.0) ? -x : x)} { if (abs($2 - $4) > 0.00001) { print "Difference detected"; exit 1} }' count-files/test-alignment-file-gene.cmp`
  if [ "$DIFF" != "" ]
  then
    echo "The computed values by the test for file are not correct - please check the output."
    exit 1
  fi
  echo "... correct"
  echo
done

