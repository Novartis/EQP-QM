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
##  Script to compute a mapping between internal and external exon ids. The
##  previous approach was to take the files just from the exon_gene_rank.map
##  file; however, this file contains only one external exon id for each
##  internal exon id. But for Ensembl there may be more than one.
##
################################################################################

import sys
import argparse
import re
import os

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "pipeline-setup-scripts", "util-scripts", "createExonExonMapFile.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Convert a junction BED file to a junction count file')
parser.add_argument('-d', default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output] ')
parser.add_argument("-b", dest='genomeExonBedFile', default="", metavar="STRING", help='BED file with genomic coordinates of exons')
parser.add_argument("-g", dest='geneExonMapFile', default="", metavar="STRING", help='Mapping of internal exon id to gene ids')
parser.add_argument('gtfFile', metavar="STRING", help='GTF file with possibly additional exon identifiers')
parser.add_argument('mapFile', metavar="STRING", help='File with internal exon id to external exon id mapping')

if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  inputArgs = []
  inputArgs.append("-b")
  inputArgs.append(os.path.join(projectDir, "chr22-genome-exon.bed"))
  inputArgs.append("-g")
  inputArgs.append(os.path.join(projectDir, "chr22_exon_gene.map"))
  
  inputArgs.append(os.path.join(projectDir, "chr22-exon-id.gtf"))
  inputArgs.append(os.path.join(projectDir, "chr22_exon_exon.map"))
  args = parser.parse_args(inputArgs)

debugLevel  = args.debugLevel
gtfFilename = args.gtfFile
mapFilename = args.mapFile

geneExonMapFilename   = args.geneExonMapFile
genomeExonBedFilename = args.genomeExonBedFile

################################################################################
##
## readGtfFile
##
################################################################################

def readGtfFile (gtfFilename):

  exonIdExonLocationMap = {}
  exonLocationExonIdMap = {}
  lineNumChunkSize = 20000
 
  try:
    gtfFile = open(gtfFilename)
    print("Reading file: " + gtfFilename, file=sys.stderr)
    lineNum = 0
    for line in gtfFile:
      line = line.rstrip ()
      gtfEntries = line.split("\t")
      
      chromosome  = gtfEntries[0]
      source      = gtfEntries[1]
      featureType = gtfEntries[2]
      interval    = [int(gtfEntries[3]), int(gtfEntries[4])]
      strand      = gtfEntries[6]
      frame       = gtfEntries[7]
      
      score = 0
      if gtfEntries[5] != ".":
        score = float(gtfEntries[5])
      
      if featureType == "exon":
        exon = "/".join(map(str, [chromosome] + interval + [strand]))
        annotationEntries = gtfEntries[8].split(";")
        chromosomeStrand = chromosome + strand
        geneId = ""
        exonId = ""
        transcriptId = ""
        exonNumber = -1
        for annotationEntry in annotationEntries:
          if annotationEntry != "":
            annotationType, annotationValue, annotationEmpty = annotationEntry.strip().split('"')
            if annotationEmpty != "":
              raise Exception ("Unknown format of annotation entry: " + annotationEntry)

            annotationType = annotationType.strip()
            if annotationType == "gene_id":
              geneId = annotationValue
            elif annotationType == "transcript_id":
              transcriptId = annotationValue
            if annotationType == "exon_id":
              exonId = annotationValue
            elif annotationType == "exon_number":
              exonNumber = int(annotationValue)

        if geneId == "":
          print("WARNING: no gene id for exon " + exon + " of transcript " + transcriptId, file=sys.stderr)
        if exonId == "":
          print("WARNING: no exon id for exon " + exon + " of gene " + geneId, file=sys.stderr)
        else:
          geneExon = geneId + "/" + exon
          exonIdExonLocationMap[exonId] = geneExon
          if not geneExon in exonLocationExonIdMap:
            exonLocationExonIdMap[geneExon] = []
          if not exonId in exonLocationExonIdMap[geneExon]:
            exonLocationExonIdMap[geneExon].append(exonId)
      
      lineNum = lineNum + 1
      if lineNum % lineNumChunkSize == 0:
        sys.stderr.write(".")
        sys.stderr.flush()
    
    if lineNum > lineNumChunkSize:
      sys.stderr.write("\n")
    
    gtfFile.close()
  
  except IOError as e:
    raise Exception("File " + gtfFilename + " not found ... exiting\n" + "Unix error message: " + str(e[0]) + " " + str(e[1]))
  
  return exonIdExonLocationMap, exonLocationExonIdMap


################################################################################
##
## Main program
##
################################################################################

exonIdExonLocationMap, exonLocationExonIdMap = readGtfFile(gtfFilename)
 
try:
  mapFile = open(mapFilename, 'w')
except IOError as e:
  raise Exception("File " + mapFilename + " cannot be opened for writing.\n" +"Unix error message: " + str(e[0]) + " " + str(e[1]))

print("Writing exon - exon pairs to file: " + mapFilename, file=sys.stderr)
for exonLocationId in sorted(exonLocationExonIdMap.keys()):
  for exonId in sorted(exonLocationExonIdMap[exonLocationId]):
    print("\t".join([exonLocationId, exonId]), file=mapFile)

mapFile.close ()

geneExonMapFilename   = args.geneExonMapFile
genomeExonBedFilename = args.genomeExonBedFile

## Output gene exon map file
if geneExonMapFilename != "":
  try:
    geneExonMapFile = open(geneExonMapFilename, 'w')
  except IOError as e:
    raise Exception("File " + geneExonMapFilename + " cannot be opened for writing.\n" +"Unix error message: " + str(e[0]) + " " + str(e[1]))
  print("Writing exon gene pairs to file: " + geneExonMapFilename, file=sys.stderr)
  
  for exonLocationId in sorted(exonLocationExonIdMap.keys()):
    geneId, chromosome, start, end, strand = exonLocationId.split("/")
    exonId = "/".join([chromosome, start, end, strand])
    print("\t".join([exonId, geneId]), file=geneExonMapFile)
    
  geneExonMapFile.close ()


## Output genome exon BED file
if genomeExonBedFilename != "":
  try:
    genomeExonBedFile = open(genomeExonBedFilename, 'w')
  except IOError as e:
    raise Exception("File " + genomeExonBedFilename + " cannot be opened for writing.\n" +"Unix error message: " + str(e[0]) + " " + str(e[1]))
  print("Writing exon genome BED entries to file: " + genomeExonBedFilename, file=sys.stderr)
  
  for exonLocationId in sorted(exonLocationExonIdMap.keys()):
    geneId, chromosome, start, end, strand = exonLocationId.split("/")
    exonId = "/".join([chromosome, start, end, strand])
    
    bedStart = str(int(start)-1)
    bedEnd   = end
    print("\t".join([chromosome, bedStart, bedEnd, exonId, "0", strand]), file=genomeExonBedFile)
    
  genomeExonBedFile.close ()
 

