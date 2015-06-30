#!/usr/bin/env python

################################################################################
##
##  Script to compute the map file between internal exon ids (with format:
##  <gene Id>/<chromosome>/<start>/<end>/<strand>) to external junction ids
##
################################################################################

import sys
import argparse
import os
import re
import os.path
import gzip

## execfile(os.path.join(os.environ["HOME"], "ngs", "pipelines", "exon-pipeline", "bin", "pipeline-setup-scripts", "util-scripts", "createExonJunctionMapFile.py"))

################################################################################
##
## Define command line options
##
## Note that the junction - junction map file can be used to filter out
## exon - exon junctions with too low a separation from the external exon - exon
## junction map file.
##
################################################################################

parser = argparse.ArgumentParser(description='Extract pipeline meta informations')
parser.add_argument('-d', type=int, default=0, dest="debugLevel", metavar="INT", help='debug level [0, no debug output]')
parser.add_argument('-A', dest="addThirdColumn", default=False, action="store_true", help='Flag to indicate whether a third column should be added')
parser.add_argument('-e', dest="exonExonMapFile", metavar="<Exon exon map file>", help='File with mapping of internal exon ids to external exon ids')
parser.add_argument('-j', dest="junctionExonMapFile", metavar="<Junction exon map file>", help='File with mapping of external junction ids to external exon ids')
parser.add_argument('-J', dest="junctionJunctionMapFile", default="", metavar="<Junction junction map file>", \
                    help='File with mapping of internal junction ids to external junction ids')
parser.add_argument('-o', dest="outputFile", metavar="<Output file>", help='File with mapping of internal exon ids to external junction ids')
parser.add_argument('-w', dest="printWarnings", action='store_true', help='Flag to enable the printing of warnings')


if len(sys.argv) > 1:
  args = parser.parse_args(sys.argv[1:])
else:
  exonPipelinePath = os.path.join("/dlab", "ldrive", "PHCHBS-I21605", "schuisv1", "ngs", "RNA-seq-test", "EQP2.0-speed-test", "exon-pipeline-files")
  exonPipelinePath = os.path.join("/dlab", "ldrive", "PHCHBS-I21605", "EQP", "EQP2.0", "SEQC-benchmark-instance-new", "human")
  mapPath = os.path.join(exonPipelinePath, "map-files")
  
  geneModel = "ensembl"
  geneModel = "refseq"
  organismShort = "hs"
  
  inputArgs = []
  inputArgs.append("-A")
  inputArgs.append("-e")
  inputArgs.append(os.path.join(mapPath, geneModel + "_rna_" + organismShort + "_exon_exon.map"))
  inputArgs.append("-j")
  inputArgs.append(os.path.join(mapPath, geneModel + "_rna_" + organismShort + "_junction_exon.map"))
  # inputArgs.append("-J")
  # inputArgs.append(os.path.join(mapPath, geneModel + "_rna_" + organismShort + "_100_junction_junction.map.new"))
  inputArgs.append("-o")
  inputArgs.append(os.path.join(mapPath, geneModel + "_rna_" + organismShort + "_exon_junction.map.new"))
  inputArgs.append("-w")
  
  print "No arguments given - using files in " + mapPath
  args = parser.parse_args(inputArgs)

debugLevel                  = args.debugLevel
junctionExonMapFilename     = args.junctionExonMapFile
junctionJunctionMapFilename = args.junctionJunctionMapFile
exonExonMapFilename         = args.exonExonMapFile
outputFilename              = args.outputFile
printWarnings               = args.printWarnings
addThirdColumn              = args.addThirdColumn


################################################################################
##
## open a file
##
################################################################################

def openFile (filename, mode="r"):
  try:
    if filename.endswith(".gz"):
      fileLink = gzip.open(filename, mode + "b")
    else:
      fileLink = open(filename, mode)

    return fileLink
      
  except IOError, e:
    raise Exception (filename + " cannot be opened ... exiting\n" +
                     "Unix error code and message: " + str(e[0]) + " " + str(e[1]))


################################################################################
##
## readJunctionJunctionMapFile
##
################################################################################

def readJunctionJunctionMapFile (junctionJunctionMapFilename):
  
  junctionJunctionMap = {}
  
  try:
    junctionJunctionMapFile = openFile(junctionJunctionMapFilename)
  except IOError:
    raise Exception ("Could not open file " + junctionJunctionMapFilename)

  print "Reading junction junction map file: " + junctionJunctionMapFilename
  i = 0
  for line in junctionJunctionMapFile:
    internalJunctionId, externalJunctionId = line.rstrip ().split ("\t")[0:2]
    junctionJunctionMap[externalJunctionId] = internalJunctionId
            
    i = i + 1
    if i % 200000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    
  junctionJunctionMapFile.close()
  if i > 200000:
    sys.stderr.write("\n")
    
  return junctionJunctionMap


################################################################################
##
## readExonExonMapFile
##
################################################################################

def readExonExonMapFile (exonExonMapFilename):
  exonExonMap = {}
  
  try:
    exonExonMapFile = openFile(exonExonMapFilename)
  except IOError:
    print "Could not open file %s" % exonExonMapFilename
    sys.exit(2)
  print "Reading exon exon map file: " + exonExonMapFilename
  i = 0
  for line in exonExonMapFile:
    internalExonId, externalExonId = line.rstrip ().split ("\t")
    geneId, chromosome, start, end, strand = internalExonId.split ("/")
    if not externalExonId in exonExonMap:
      exonExonMap[externalExonId] = {}
    if not geneId in exonExonMap[externalExonId]:
      exonExonMap[externalExonId][geneId] = set([])
    exonExonMap[externalExonId][geneId].add (internalExonId)
            
    i = i + 1
    if i % 50000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()
    
  exonExonMapFile.close()
  if i > 50000:
    sys.stderr.write("\n")
    
  return exonExonMap


################################################################################
##
## Concatenate exon sequences
##
################################################################################

if junctionJunctionMapFilename != "":
  junctionJunctionMap = readJunctionJunctionMapFile (junctionJunctionMapFilename)
exonExonMap = readExonExonMapFile (exonExonMapFilename)

try:
  junctionExonMapFile = openFile (junctionExonMapFilename)
except IOError:
  print "Could not open file %s" % junctionExonMapFilename
  sys.exit(2)

try:
  outputFile = openFile (outputFilename, 'w')
except IOError:
  print "Could not open file %s" % outputFilename
  sys.exit(2)

print "Reading file " + junctionExonMapFilename + " and"
print "writing to file " + outputFilename
excludedExons = {}
i = 0
for line in junctionExonMapFile:
  
  externalJunctionId, externalExonId1, externalExonId2 = line.rstrip (). split ("\t")

  if externalExonId1 in exonExonMap and externalExonId2 in exonExonMap:
    if junctionJunctionMapFilename == "" or externalJunctionId in junctionJunctionMap:
      for geneId in [geneId for geneId in exonExonMap[externalExonId1] if geneId in exonExonMap[externalExonId2]]:
        for internalExonId1 in exonExonMap[externalExonId1][geneId]:
          for internalExonId2 in exonExonMap[externalExonId2][geneId]:
            if not addThirdColumn:
              print >> outputFile, "\t".join([internalExonId1, externalJunctionId])
              print >> outputFile, "\t".join([internalExonId2, externalJunctionId])
            else:
              print >> outputFile, "\t".join([internalExonId1, externalJunctionId, internalExonId1 + "+" + internalExonId2])
              print >> outputFile, "\t".join([internalExonId2, externalJunctionId, internalExonId1 + "+" + internalExonId2])
            
    i = i + 1    
    if i % 200000 == 0:
      sys.stderr.write(".")
      sys.stderr.flush()

  else:
    if not externalExonId1 in exonExonMap and not externalExonId1 in excludedExons:
      if printWarnings:
        print >> sys.stderr, "WARNING: exon " + externalExonId1 + " has been excluded"
      excludedExons[externalExonId1] = 1
    if not externalExonId2 in exonExonMap and not externalExonId2 in excludedExons:
      if printWarnings:
        print >> sys.stderr, "WARNING: exon " + externalExonId2 + " has been excluded"
      excludedExons[externalExonId2] = 1
            
junctionExonMapFile.close()
outputFile.close ()
if i > 50000:
  sys.stderr.write("\n")
  
print "Mapping information for " + str(2*i) + " internal exon id/external junction id pairs written."

  

