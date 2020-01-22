#!/usr/bin/env python

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

################################################################################
##
##  Script to add exon_id fields to a GTF file
##
################################################################################

import sys
import argparse
import os

## execfile(os.path.join("/dlab", "ldrive", "PHCHBS-I21605", "ref-files", "Laurent-chr22", "bin", "addExonIdGtf.py"))


################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Add an exon_id field to the exon annotation in a GTF file')
parser.add_argument('-d', default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-g', dest='gtfFile', default='-', metavar="<GTF file>", help='GTF file for the exons and transcripts.')
parser.add_argument('-o', dest='outputFile', default='-', metavar="<Output file>", help='GTF output file with the exon numbers')
parser.add_argument('-w', dest="warningsOn", action='store_true', help='Flag to set warnings on')

if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  pipelineDir = os.path.join("/dlab", "ldrive", "PHCHBS-I21605", "ref-files", "Laurent-chr22")
  outputDir   = pipelineDir
  
  inputArgs = []
  inputArgs.append('-g')
  inputArgs.append(os.path.join(pipelineDir, "chr22.gtf"))
      
  inputArgs.append('-o')
  inputArgs.append(os.path.join(pipelineDir, "chr22-test.gtf"))
  
  inputArgs.append('-w')
  
  args = parser.parse_args(inputArgs)

debugLevel     = args.debugLevel
gtfFilename    = args.gtfFile
outputFilename = args.outputFile
warningsOn     = args.warningsOn

if len(sys.argv) <= 1:
  print("Using default arguments: ", file=sys.stderr)
  print("gtfFilename: " + gtfFilename, file=sys.stderr)
  print("outputFilename: " + outputFilename, file=sys.stderr)

lineNumChunkSize = 50000

################################################################################
##
## Main program
##
################################################################################

try:
  gtfFile = open(gtfFilename)
  print("Reading file: " + gtfFilename, file=sys.stderr)
except IOError as e:
  raise Exception ("File " + gtfFilename + " not found\n" + "Unix Error: " + str(e[0]) + " " + str(e[1]))

try:
  outputFile = open(outputFilename, 'w')
  print("Writing to file: " + outputFilename, file=sys.stderr)
except IOError as e:
  raise Exception ("File " + outputFilename + " not found\n" + "Unix error: " + str(e[0]) + " " + str(e[1]))

lineNum  = 0
numExons = 0
for line in gtfFile:
  if line.startswith("#"):
    continue
  
  line = line.rstrip ()
  # print "line: " + line
  gtfEntries = line.split("\t")

  chromosome  = gtfEntries[0]
  source      = gtfEntries[1]
  featureType = gtfEntries[2]
  start       = int(gtfEntries[3])
  end         = int(gtfEntries[4])
  interval    = [start, end]
  score       = float(0 if gtfEntries[5] == "." else gtfEntries[5])
  strand      = gtfEntries[6]
  frame       = gtfEntries[7]

  if featureType == "exon":
    numExons += 1
    if not "exon_id" in gtfEntries[8]:
      externalExonId = ":".join(["exon", chromosome, "-".join(map(str, interval)), strand])
      gtfEntries[8] = gtfEntries[8] + ' exon_id "' + externalExonId + '";'

    print("\t".join(gtfEntries), file=outputFile)
  
  lineNum = lineNum + 1
  if lineNum % lineNumChunkSize == 0:
    sys.stderr.write(".")
    sys.stderr.flush()

## End the line of dots
if lineNum > lineNumChunkSize:
  sys.stderr.write("\n")

## Close the files
gtfFile.close()
outputFile.close ()
print(str(numExons) + " exon entries found.", file=sys.stderr)


