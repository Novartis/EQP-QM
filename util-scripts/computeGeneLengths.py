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
##  Script to compute the length of genes and transcripts based on GTF file
##
##  Assumptions: GTF file entries of the same transcript are sorted by ascending
##               coordinates for transcripts on the plus strand and descending
##               coordinates for transcripts on the minus strand; GTF filr
##               contains exonNumber entries for each exon.
##
##  Issues:
##
##  Genes and transcripts can map to different chromosomes and strands or
##  to multiple locations on the same strand of a chromosome. We use
##  the exonNumber of the GTF file to distinguish between different locations.
##
################################################################################

import sys
import argparse
import re
import os.path
import numpy

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "pipeline-setup-scripts", "util-scripts", "computeGeneLengths.py"))

################################################################################
##
## Define command line options
##
################################################################################

parser = argparse.ArgumentParser(description='Combine the count files from different')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-g', dest="gtfFilename", metavar="<GTF file>",  help='GTF file with the exon coordinates')
parser.add_argument('-o', dest="geneLengthFile", default="", metavar="<Gene len. file>", help='File with the lengths of genes')
parser.add_argument('-w', dest="warningsOn", action='store_true', help='Flag to set warnings on')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  organism = "human"
  organism_short="hs"
  
  # organism = "mouse"
  # organism_short="mm"
  
  # gene_model="ensembl"
  gene_model="refseq"

  
  inputArgs = []
  inputArgs.append("-g")
  inputArgs.append(os.path.join(projectDir, "gtf-files", gene_model + "_rna_" + organism_short + ".gtf"))
  inputArgs.append("-o")
  inputArgs.append(os.path.join(projectDir, "gtf-files", gene_model + "_rna_" + organism_short + "-gene-lengths.txt"))
  args = parser.parse_args(inputArgs)

debugLevel     = args.debugLevel
gtfFilename    = args.gtfFilename
geneLengthFile = args.geneLengthFile
warningsOn     = args.warningsOn


################################################################################
##
## head
##
################################################################################

def head (x, n=10):
  
  if len(x) == 1 or isinstance(x, str):
    return x
  
  if isinstance(x, dict):
    L = {}
  else:
    L = []
    
  for y in sorted(x)[0:n]:
    if isinstance(x, dict):
      L[y] = x[y]
    else:
      L.append(y)
      
  return L


################################################################################
##
## Check for overlap of an interval with a list of intervals
##
################################################################################
      
def overlaps (interval1, interval2):
  
  if (interval1[0] <= interval2[0] and interval2[0] <= interval1[1]) or (interval2[0] <= interval1[0] and interval1[0] <= interval2[1]):
    return True
  
  return False




################################################################################
##
## union
##
################################################################################

def union (interval1, interval2):

  return [min (interval1[0], interval2[0]), max (interval1[1], interval2[1])]


################################################################################
##
## contains
##
################################################################################

def contains (interval1, interval2):

  return (interval1[0] <= interval2[0]) and (interval1[1] >= interval2[1])


################################################################################
##
## Intersection
##
################################################################################
      
def intersection (interval1, interval2):
  
  leftEndPoint  = max (interval1[0], interval2[0])
  rightEndPoint = min (interval1[1], interval2[1])

  if leftEndPoint <= rightEndPoint:
    return rightEndPoint - leftEndPoint + 1

  return 0


################################################################################
##
## unique
##
################################################################################

def unique (seq, idfun=None):
  
   # order preserving
   if idfun is None:
       def idfun(x): return x
       
   seen = {}
   result = []
   for item in seq:
     marker = idfun (item)
     # in old Python versions:
     if marker in seen: continue
     seen[marker] = 1
     result.append(item)
   return result


################################################################################
##
## getOverlapLength
##
################################################################################
      
def getOverlapLength (interval1, interval2):
  
  leftEndPoint  = max([min(interval1), min(interval2)])
  rightEndPoint = min([max(interval1), max(interval2)])
       
  length1 = max(interval1) - min(interval1) + 1
  length2 = max(interval2) - min(interval2) + 1
  maxDist = max (length1/2, length2/2)

  ## We consider intervals that are less than half their lengths
  ## apart also overlapping
  return rightEndPoint - leftEndPoint + 1 + maxDist


################################################################################
##
## intervalListOverlap
##
################################################################################

def intervalListsOverlap (intervalList1, intervalList2):

  interval1 = [min(min(intervalList1)), max(max(intervalList1))]
  interval2 = [min(min(intervalList2)), max(max(intervalList2))]

  return getOverlapLength (interval1, interval2) >= 0
  

################################################################################
##
## readGtfFile
##
## Use the exonNumber fields of each transcript to infer how many alignments a
## transcript has. This is captured by the variable alignmentNumber.
##
## Each of the three transcript output variables has an entry for each each transcript
## and each alignmentNumber of each transcript.
##
################################################################################

def readGtfFile (gtfFilename):
  
  transcriptExonList         = {}
  transcriptChromosomeStrand = {}

  geneTranscriptMap     = {}
  
  transcriptAlignmentExonNumbers = {}

  lineNumChunkSize = 50000
 
  try:
    gtfFile = open(gtfFilename)
  except IOError, e:
    raise Exception(gtfFilename + " cannot be opened ... skipping\n" + 
                    "Unix error code and message: " + str(e[0]) + " " + str(e[1]))

  print >> sys.stderr, "Opening file: " + gtfFilename
  lineNum = 0
  oldTranscriptId = ""
  alignmentNumber = 1

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
    
    exonId      = "/".join(map(str, [chromosome] + interval + [strand]))
    
    if featureType == "exon":
      annotationEntries = gtfEntries[8].split(";")
      geneId = ""
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
          elif annotationType == "exon_number":
            exonNumber = int(annotationValue)
            
      if oldTranscriptId != transcriptId:
        alignmentNumber = 1
        oldTranscriptId = transcriptId
        
      if exonNumber == -1:
        raise Exception ("WARNING: no exon number for exon " + exonId + " of transcript " + transcriptId)
               
      if not transcriptId in transcriptExonList:
        transcriptExonList            [transcriptId] = {}
        transcriptChromosomeStrand    [transcriptId] = {}
        transcriptAlignmentExonNumbers[transcriptId] = {}

      ## Increment the alignment number if the current exon number already exists in the current alignment
      ## Note that GTF entries of other transcripts may occur between two alignments of the same transcript
      while alignmentNumber in transcriptAlignmentExonNumbers[transcriptId] and \
             exonNumber in transcriptAlignmentExonNumbers[transcriptId][alignmentNumber]:
        alignmentNumber += 1

      if not alignmentNumber in transcriptAlignmentExonNumbers[transcriptId]:
        transcriptAlignmentExonNumbers[transcriptId][alignmentNumber] = [exonNumber]
      else:
        transcriptAlignmentExonNumbers[transcriptId][alignmentNumber].append(exonNumber)

      if not alignmentNumber in transcriptChromosomeStrand[transcriptId]:
        #print "transcriptChromosomeStrand of " + transcriptId + "." + str(alignmentNumber) + " set to " + chromosome + strand
        transcriptChromosomeStrand[transcriptId][alignmentNumber] = chromosome + strand
      elif transcriptChromosomeStrand[transcriptId][alignmentNumber] != chromosome + strand:
        print >> sys.stderr, "WARNING: Exon number " + str(exonNumber) + " of transcript " + transcriptId + " on chromosome/strand " + chromosome + strand + \
              " is assigned to alignment " + str(alignmentNumber) + " on chromosome/strand " + transcriptChromosomeStrand[transcriptId][alignmentNumber]
      
      if alignmentNumber in transcriptExonList[transcriptId]:
        transcriptExonList[transcriptId][alignmentNumber].append(interval)
      else:
        transcriptExonList[transcriptId][alignmentNumber] = [interval]
        
      if exonNumber in transcriptExonList[transcriptId][alignmentNumber]:
        print >> sys.stderr ("Exon number: " + str(exonNumber) + " already stored for alignment " + str(alignmentNumber) + " of " + transcriptId)
        sys.exit(1)
      
      if geneId != "" and transcriptId != "":
        if geneId in geneTranscriptMap:
          if not transcriptId in geneTranscriptMap[geneId]:
            geneTranscriptMap[geneId].append(transcriptId)
        else:
          geneTranscriptMap[geneId] = [transcriptId]
    
    lineNum = lineNum + 1
    if lineNum % lineNumChunkSize == 0:
      sys.stdout.write(".")
      sys.stdout.flush()
  
  if lineNum > lineNumChunkSize:
    sys.stdout.write("\n")
  
  gtfFile.close()
  return transcriptChromosomeStrand, transcriptExonList, geneTranscriptMap


################################################################################
##
## computeIntervalListLength
##
################################################################################

def computeIntervalListLength (intervalList):
  
  intervalList.sort(key=lambda tup: tup[0])
  oldInterval = [-1,-1]
  combinedList = []
  for interval in sorted(intervalList):
    if oldInterval[1] < interval[0]:
      combinedList.append(oldInterval)
      oldInterval = interval
    else:
      oldInterval = [oldInterval[0], max(oldInterval[1], interval[1])]
  
  combinedList.append(oldInterval)
  
  length = 0
  for interval in combinedList[1:]:
    #print length
    #print (interval[1] - interval[0] + 1)
    length = length + interval[1] - interval[0] + 1
  
  return length


################################################################################
##
## sumIntervals
##
################################################################################

def computeSumIntervalLengths (countObjects):
  sumIntervalLengthMap = {}
  
  for countObjectId in countObjects:
    
    if not countObjectId in sumIntervalLengthMap:
      sumIntervalLengthMap[countObjectId] = {}

    countObject = countObjects[countObjectId]
    for alignmentNumber in countObject:
      sumIntervalLength = 0
      for interval in countObjects[countObjectId][alignmentNumber]:
        sumIntervalLength = sumIntervalLength + interval[1] - interval[0] + 1
      sumIntervalLengthMap[countObjectId][alignmentNumber] = sumIntervalLength
          
  return sumIntervalLengthMap



################################################################################
##
## findAlignmentMatches
##
## For a given query interval I on a chromosome/strand and a list of gene intervals
## look for gene intervals with an non-negative overlap on the same chromosome/strand.
##
################################################################################

def findAlignmentMatches (interval, chromosomeStrand, geneIntervals, geneChromosomeStrand, assignedGeneAlignmentNumbers):
  
  alignmentNumbers = []
  for alignmentNumber in [alignNum for alignNum in geneIntervals if not alignNum in assignedGeneAlignmentNumbers]:
    if chromosomeStrand == geneChromosomeStrand[alignmentNumber]:
      geneInterval  = geneIntervals[alignmentNumber]
      overlapLength = getOverlapLength (interval, geneInterval)
      if overlapLength >= 0:
        alignmentNumbers.append(alignmentNumber)
          
  return alignmentNumbers


################################################################################
##
## assignTranscriptAlignments
##
## Assign the alignments of the transcripts to the alignments of the gene. This
## is done by find the best match (via findBestAlignmentMatch) of each transcript
## alignment.
##
################################################################################

def assignTranscriptAlignments (geneId, geneTranscriptList, transcriptIntervals, transcriptExonList, transcriptChromosomeStrand):
  
  geneIntervals = {}
  geneExonList  = {}
  geneChromosomeStrand = {}
  geneAlignmentTranscripts = {}
  
  maxGeneAlignmentNumber = 0
  
  for transcriptId in geneTranscriptList:
    
    assignedGeneAlignmentNumbers = []
    
    if len(geneIntervals) == 0:
      for alignmentNumber in sorted(transcriptIntervals[transcriptId]):
        geneIntervals           [alignmentNumber] = transcriptIntervals       [transcriptId][alignmentNumber]
        geneExonList            [alignmentNumber] = transcriptExonList        [transcriptId][alignmentNumber]
        geneChromosomeStrand    [alignmentNumber] = transcriptChromosomeStrand[transcriptId][alignmentNumber]
        geneAlignmentTranscripts[alignmentNumber] = [transcriptId]
        
      maxGeneAlignmentNumber = max(transcriptIntervals[transcriptId])
    else:
      for alignmentNumber in transcriptIntervals[transcriptId]:
        ## the overlap between a transcript interval and a gene interval is used to determine
        ## the quality of a match.
        geneAlignmentNumbers = findAlignmentMatches (transcriptIntervals[transcriptId][alignmentNumber], transcriptChromosomeStrand[transcriptId][alignmentNumber], \
                                                     geneIntervals, geneChromosomeStrand, [])
        
        if len(geneAlignmentNumbers) > 0:
          
          geneAlignmentNumber = min (geneAlignmentNumbers)
          
          if len(geneAlignmentNumbers) > 1:
            ## If there is more than one match, then this transcript spans two previous (almost) disjoint transcripts
            ## and we need to unite them.
            geneAlignmentNumbers.remove (geneAlignmentNumber)
            for curGeneAlignmentNumber in geneAlignmentNumbers:
              geneIntervals[geneAlignmentNumber] = union (geneIntervals[geneAlignmentNumber], geneIntervals[curGeneAlignmentNumber])
              del geneIntervals[curGeneAlignmentNumber]
              geneExonList [geneAlignmentNumber].extend(geneExonList [curGeneAlignmentNumber])
              del geneExonList [curGeneAlignmentNumber]
              geneAlignmentTranscripts[geneAlignmentNumber].extend(geneAlignmentTranscripts[curGeneAlignmentNumber])
              del geneAlignmentTranscripts[curGeneAlignmentNumber]
              if geneChromosomeStrand[geneAlignmentNumber] != geneChromosomeStrand[curGeneAlignmentNumber]:
                print >> sys.stderr, "Chromosome strands of assigned gene alignment numbers of gene " + geneId + " differ: " + \
                      geneChromosomeStrand[geneAlignmentNumber] + " vs. " + geneChromosomeStrand[curGeneAlignmentNumber]
                sys.exit(1)
              else:
                del geneChromosomeStrand[curGeneAlignmentNumber]
                
          if geneAlignmentNumber in assignedGeneAlignmentNumbers and warningsOn:
            print >> sys.stderr, "\nWARNING: Alignment " + str(alignmentNumber) + " of transcript: " + transcriptId + \
                  " is also assigned to gene alignment number " + str(geneAlignmentNumber) + " of gene " + geneId
            for geneTranscriptId in geneTranscriptList:
              print >> sys.stderr, transcriptIntervals[geneTranscriptId]
            #raise Exception ("Stop")
          
          assignedGeneAlignmentNumbers.append(geneAlignmentNumber)
          
          geneIntervals[geneAlignmentNumber] = union (geneIntervals[geneAlignmentNumber], transcriptIntervals[transcriptId][alignmentNumber])
          geneExonList [geneAlignmentNumber].extend(transcriptExonList[transcriptId][alignmentNumber])
          geneAlignmentTranscripts[geneAlignmentNumber].append(transcriptId)
          
          if len(transcriptIntervals[transcriptId]) > 1:
            if (geneChromosomeStrand[geneAlignmentNumber] != transcriptChromosomeStrand[transcriptId][alignmentNumber] or \
                not contains (geneIntervals[geneAlignmentNumber], transcriptIntervals[transcriptId][alignmentNumber])) and warningsOn:
              print >> sys.stderr, "Gene interval " + geneChromosomeStrand[geneAlignmentNumber] + ":" +  str(geneIntervals[geneAlignmentNumber]) + \
                  " does not contain transcript interval " + transcriptChromosomeStrand[transcriptId][alignmentNumber] + ":" + \
                  str(transcriptIntervals[transcriptId][alignmentNumber]) + " for transcript " + transcriptId + "." + str(alignmentNumber)
            
        elif len(geneAlignmentNumbers) == 0:
          ## There is no good match between transcript and gene alignments
          maxGeneAlignmentNumber += 1
          geneAlignmentNumber = maxGeneAlignmentNumber
          
          geneIntervals           [geneAlignmentNumber] = transcriptIntervals       [transcriptId][alignmentNumber]
          geneExonList            [geneAlignmentNumber] = transcriptExonList        [transcriptId][alignmentNumber]
          geneChromosomeStrand    [geneAlignmentNumber] = transcriptChromosomeStrand[transcriptId][alignmentNumber]
          geneAlignmentTranscripts[geneAlignmentNumber] = [transcriptId]
          
  return geneExonList, geneChromosomeStrand, geneAlignmentTranscripts


################################################################################
##
## computeGeneLength:
##
## Create a sequence of exon interval list from the transcript exon interval lists
## for loci that are in common.
##
## In order to do this we need to assign each exon list of a transcript to an
## alignment of the gene.
##
## The length of the resulting gene alignment exon list is then computed as the genomic
## footprint of the (possibly overlapping) exons
##
################################################################################

def computeGeneLength (geneId, geneTranscriptList, transcriptExonList, transcriptChromosomeStrand, transcriptIntervals, transcriptLengths):
  
  geneExonList, geneChromosomeStrand, geneAlignmentTranscripts = \
          assignTranscriptAlignments (geneId, geneTranscriptList, transcriptIntervals, transcriptExonList, transcriptChromosomeStrand)

  maxNumTranscriptsPerGeneAlignment = max([len(geneAlignmentTranscripts[alignmentNumber]) for alignmentNumber in geneAlignmentTranscripts])

  if len(set(geneChromosomeStrand.values())) > 1 and warningsOn:
    print >> sys.stderr, "WARNING: gene " + geneId + " has alignments to " + str(len(set(geneChromosomeStrand.values()))) + " chromosome/strands"
  
  geneLengths = {}
  for alignmentNumber in geneExonList:
    geneLengths[alignmentNumber] = computeIntervalListLength (geneExonList[alignmentNumber])
    if len(geneAlignmentTranscripts[alignmentNumber]) != maxNumTranscriptsPerGeneAlignment and warningsOn:
      print >> sys.stderr, "WARNING: alignment " + str(alignmentNumber) + " of gene " + geneId + " on " + \
            geneChromosomeStrand[alignmentNumber] + " is only supported by " + \
            str(len(geneAlignmentTranscripts[alignmentNumber])) + " of " + str(maxNumTranscriptsPerGeneAlignment) + " transcripts"
      
  if min(geneLengths.values()) != max(geneLengths.values()) and warningsOn:
    print >> sys.stderr, "WARNING: gene length values for " + geneId + " vary from " + str(min(geneLengths.values())) + " to " + str(max(geneLengths.values()))
    
  return max(geneLengths.values())


################################################################################
##
## Main program
##
## Read GTF file and compute transcript length lists
##
################################################################################

#raise Exception ("stop")

transcriptChromosomeStrand, transcriptExonList, geneTranscriptMap = readGtfFile(gtfFilename)

################################################################################
##
## Compute the interval spanned by the transcript. The first min returns the minimum exon interval,
## the second min the first coordinate of this interval. The corresponding statement holds for max.
##
################################################################################

transcriptIntervals = {}
for transcriptId in transcriptExonList:
  transcriptIntervals[transcriptId] = {}
  for alignmentNumber in transcriptExonList[transcriptId]:
    transcriptIntervals[transcriptId][alignmentNumber] = [min(min(transcriptExonList[transcriptId][alignmentNumber])), \
                                                          max(max(transcriptExonList[transcriptId][alignmentNumber]))]


################################################################################
##
##  Compute the maximum transcript lengths and check the length of the different
##  alignments of the transcripts
##
################################################################################

transcriptLengthList = computeSumIntervalLengths (transcriptExonList)

maxTranscriptLengths = {}
for transcriptId in transcriptLengthList:
  transcriptLengths = [transcriptLengthList[transcriptId][alignmentNumber] for alignmentNumber in transcriptLengthList[transcriptId]]
  maxTranscriptLengths[transcriptId] = max(transcriptLengths)
  if max(transcriptLengths) != min(transcriptLengths) and warningsOn:\
     print >> sys.stderr, "WARNING: Length of alignments for transcript " + transcriptId + " range from " + str(min(transcriptLengths)) + \
       " to " + str(max(transcriptLengths))


################################################################################
##
##  Invert the transcriptChromosomeStrand list
##
################################################################################

transcriptAlignmentNumbers = {}
for transcriptId in transcriptChromosomeStrand:
 if not transcriptId in transcriptAlignmentNumbers:
   transcriptAlignmentNumbers[transcriptId] = {}
 for alignmentNumber in transcriptChromosomeStrand[transcriptId]:
   chromosomeStrand = transcriptChromosomeStrand[transcriptId][alignmentNumber]
   if not chromosomeStrand in transcriptAlignmentNumbers[transcriptId]:
     transcriptAlignmentNumbers[transcriptId][chromosomeStrand] = []
   transcriptAlignmentNumbers[transcriptId][chromosomeStrand].append(alignmentNumber)


################################################################################
##
##  Count the number of transcripts on multiple chromosome/strands and
##  the ones with multiple alignments on the same chromosome/strand
##
################################################################################

transcriptsOnMultipleChromosomeStrands      = 0
transcriptsWithMultipleAlignmentsSameStrand = 0
transcriptsWithOverlappingAlignments        = []
transcriptsWithMultipleAlignments           = 0
for transcriptId in transcriptAlignmentNumbers:
  if len(transcriptChromosomeStrand[transcriptId]) > 1:
    transcriptsWithMultipleAlignments += 1
  if len(transcriptAlignmentNumbers[transcriptId]) > 1:
    transcriptsOnMultipleChromosomeStrands += 1
  for chromosomeStrand in transcriptAlignmentNumbers[transcriptId]:
    if len(transcriptAlignmentNumbers[transcriptId][chromosomeStrand]) > 1:
      transcriptsWithMultipleAlignmentsSameStrand += 1
    for alignmentNumber1 in transcriptAlignmentNumbers[transcriptId][chromosomeStrand]:
      for alignmentNumber2 in transcriptAlignmentNumbers[transcriptId][chromosomeStrand]:
        if alignmentNumber1 != alignmentNumber2:
          if overlaps(transcriptIntervals[transcriptId][alignmentNumber1], transcriptIntervals[transcriptId][alignmentNumber2]):
            if not transcriptId in transcriptsWithOverlappingAlignments:
              transcriptsWithOverlappingAlignments.append(transcriptId)

if transcriptsWithMultipleAlignments > 0:
  print "There are " + str(transcriptsWithMultipleAlignments) + " transcripts with multiple alignments:"
  print " - " + str(transcriptsOnMultipleChromosomeStrands) + " transcripts have alignments to different chromosome/strand pairs"
  print " - " + str(transcriptsWithMultipleAlignmentsSameStrand) + " transcripts have multiple alignments to the same chromosome/strand"


################################################################################
##
## Compute the gene length for one genomic location as sum of the intervals of
## the alignments of all the gene's transcripts that overlap one genomic location.
## The overall gene length is then the maximum over all genomic locations.
##
################################################################################

geneLengths = {}
for geneId in geneTranscriptMap:
  geneTranscriptList = geneTranscriptMap[geneId]
  if len(geneTranscriptList) == 1:
    transcriptId        = geneTranscriptList[0]
    geneLengths[geneId] = maxTranscriptLengths[transcriptId]
  else:
    ## Compute the number of chromosome/strands that the transcripts of a gene
    ## align to.
    geneChromosomeStrandList = []
    for transcriptId in geneTranscriptList:
      if not transcriptChromosomeStrand[transcriptId][1] in geneChromosomeStrandList:
        geneChromosomeStrandList.append(transcriptChromosomeStrand[transcriptId][1])
        
    ## For data sets (like Ensembl) which do not have transcripts with multiple alignments
    ## and for which all transcripts of a gene map to the same chromosome/strand we consider
    ## all transcripts to be part of the same locus and base the gene length computation on
    ## all the exons in this locus.
    if len(geneChromosomeStrandList) == 1 and transcriptsWithMultipleAlignments == 0:
      geneExonList = []
      for transcriptId in geneTranscriptList:
        geneExonList.extend(transcriptExonList[transcriptId][1])
      geneLengths[geneId] = computeIntervalListLength (geneExonList)
    else:
      geneLengths[geneId] = computeGeneLength (geneId, geneTranscriptList, transcriptExonList, transcriptChromosomeStrand, transcriptIntervals, maxTranscriptLengths)
    
    maxTranscriptLengthGene = max([maxTranscriptLengths[transcriptId] for transcriptId in geneTranscriptList])
    if maxTranscriptLengthGene > geneLengths[geneId] and warningsOn:
      print "WARNING: Gene " + geneId + " has transcripts that are " + str(maxTranscriptLengthGene - geneLengths[geneId]) + " longer than " + str(geneLengths[geneId])


## Check "normal" genes, that is, genes with transcripts that align only to one chromosome strand and
## only to one location on a chromosome strand
numNormalGenes = 0
for geneId in geneTranscriptMap:
  normalGene = True
  geneIntervalList = []
  geneChromosomeStrand = ""
  geneTranscriptList = geneTranscriptMap[geneId]
  for transcriptId in geneTranscriptList:
    if normalGene and len(transcriptExonList[transcriptId]) == 1:
      for alignmentNumber in transcriptExonList[transcriptId]:
        if geneChromosomeStrand == "":
          geneChromosomeStrand = transcriptChromosomeStrand[transcriptId][alignmentNumber]
        if normalGene and geneChromosomeStrand == transcriptChromosomeStrand[transcriptId][alignmentNumber]:
          if geneIntervalList != [] and transcriptsWithMultipleAlignments > 0 and not intervalListsOverlap(geneIntervalList, transcriptExonList[transcriptId][alignmentNumber]):
            normalGene = False
          else:
            geneIntervalList = geneIntervalList + transcriptExonList[transcriptId][alignmentNumber]
        else:
          normalGene = False
    else:
      normalGene = False
  
  if normalGene:
    numNormalGenes = numNormalGenes + 1
    geneLength = computeIntervalListLength (geneIntervalList)
    if geneLength != geneLengths[geneId]:
      print "Differing gene lengths for gene " + geneId + ": " + str(geneLength) + " vs. " + str(geneLengths[geneId])
    

################################################################################
##
## Output the gene lengths
##
################################################################################
  
try:
  print "Writing lengths to " + geneLengthFile
  outputFile = open(geneLengthFile, 'w')
  
  print >> outputFile, "\t".join(["Gene Id", "Length"])
  for geneId in sorted(geneLengths.iterkeys()):
    print >> outputFile, geneId + "\t" + str(geneLengths[geneId])

  outputFile.close()
  
except IOError, e:
  print geneLengthFile + " cannot be opened ... skipping"
  print "Error: " + str(e[0]) + " " + str(e[1])
    
 


