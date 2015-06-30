#call with: source(file.path(Sys.getenv("HOME"), "ngs", "pipelines", "eqp-qm-public", "trunc", "util-scripts", "R", "create-junction-exon-map-file.R"))

pipeline.dir <- Sys.getenv("EXON_PIPELINE_DIR")

pipeline_R.dir <- file.path(pipeline.dir, "R")

source(file.path(pipeline_R.dir, "exon-util-lib.R"))
source(file.path(pipeline_R.dir, "pipeline-setup-files", "create-exon-gene-rank-util-lib.R"))
source(file.path(pipeline_R.dir, "project-setup-files", "junction-util-lib.R"))

###############################################################################
##
## noOverlap
##
###############################################################################

noOverlap <- function (exonId1, exonId2) {

  exonId1.fields <- strsplit (exonId1, "/")[[1]]
  exonId2.fields <- strsplit (exonId2, "/")[[1]]

  chromosomeId1 <- exonId1.fields[1]
  exonStart1    <- as.integer(exonId1.fields[2])
  exonEnd1      <- as.integer(exonId1.fields[3])
  exonStrand1   <- exonId1.fields[4]

  chromosomeId2 <- exonId2.fields[1]
  exonStart2    <- as.integer(exonId2.fields[2])
  exonEnd2      <- as.integer(exonId2.fields[3])
  exonStrand2   <- exonId2.fields[4]

  ## Note that if two exons belong to the same gene but to two
  ## different chromosomes or strands we do not want to consider
  ## the junction between them so we just say that the exons overlap.
  if (exonStrand1 != exonStrand2 || chromosomeId1 != chromosomeId2) {
    return (FALSE)
  }

  ## Ensure that there are at least two bases separating the two exons
  if (exonEnd1 < exonStart2 - 2 || exonEnd2 < exonStart1 - 2) {
    return (TRUE)
  }

  return (FALSE)
  
}

###############################################################################
##
## Main program
##
###############################################################################

###################### Read GTF file ##########################################

## Read gtf file
gtf.file <-  Sys.getenv("GTF_FILE")
gtf.frame.all <- readGffFile(gtf.file)
gtf.frame <- subset(gtf.frame.all, gtf.frame.all[["Feature type"]] == "exon")

## Create frame with exon information only
gtf.chromosome.strand.ids <- paste(gtf.frame[["Chromosome"]], gtf.frame[["Start"]], gtf.frame[["End"]], gtf.frame[["Strand"]], sep="/")
gtf.gene.ids <- getAttributeFieldGff(gtf.frame[["Attributes"]], "gene_id")
gtf.exon.ids <- getAttributeFieldGff(gtf.frame[["Attributes"]], "exon_id")

if (sum(is.na(gtf.exon.ids)) > 0) {
  print(paste("No exon ids available for", sum(is.na(gtf.exon.ids)), "of", length(gtf.exon.ids), "exons."))
  print("Creating artificial exons.")
  gtf.exon.ids[is.na(gtf.exon.ids)] <-
    paste("exon", gtf.frame[["Chromosome"]], paste(gtf.frame[["Start"]], gtf.frame[["End"]], sep="-"), sep=":")[is.na(gtf.exon.ids)]
}


###############################################################################
##
## Create exon gene rank frame
##
###############################################################################

exon.gene.rank.frame <- createExonGeneRankMapFrame (gtf.chromosome.strand.ids, gtf.gene.ids, gtf.exon.ids)
print(paste(nrow(exon.gene.rank.frame), "exon gene rank mappings found for", length(unique(exon.gene.rank.frame[["Exon Id"]])), "exons."))

radius      <- 30
maxExon.num <- 1

#stopifnot(FALSE)

###################################################################################
##
## Compute junctions
##
###################################################################################

exon.chr.ids        <- as.character(exon.gene.rank.frame[["Exon Chromosome Id"]])
names(exon.chr.ids) <- as.character(exon.gene.rank.frame[["Gene Id"]])

exon.gene.rank.ind.list  <- tapply(exon.chr.ids, names(exon.chr.ids), getSortedExonsInd)
exon.id.list             <- tapply(as.character(exon.gene.rank.frame[["Exon Id"]]), names(exon.chr.ids), c)
exon.chr.id.list         <- tapply(as.character(exon.gene.rank.frame[["Exon Chromosome Id"]]), names(exon.chr.ids), c)

print(paste("Finding exon indices"))
exon.sep <- 3
exon.gene.rank.index.list <- tapply(exon.chr.ids, names(exon.chr.ids), findExonIndices, maxExon.num, radius, exon.sep)

## exon.gene.rank.index.list list has the following structure
##
## for each gene i (with n exons):
##
## exon.gene.rank.index.list[[i]]
##  -> left                     -> right
##     (n-1 element list,          (n-1 element list,
##      one for each exon           one for each exon
##      except the last)            except the first)
##     Element i (corresponds      Element i (corresponds
##     to exon i):                 to exon i+1)
##       -> list of vectors         -> list of vectors
##          each ends with i           each starts with i + 1

###################################################################################
##
## Reduce the number of allowed exons on each side of the junction to reduce the
## number of junctions for genes with many (short) exons
##
###################################################################################

print(paste("Computing combinations"))
numCombinations <- sum(sapply(exon.gene.rank.index.list, getGeneCombinationNum))

junction.id  <- array("", dim=numCombinations)
exon.id1     <- array("", dim=numCombinations)
exon.id2     <- array("", dim=numCombinations)
junction.ind <- 1
for (i in 1:length(exon.gene.rank.index.list)) {
  old.left.exon.ind  <- -1
  old.right.exon.ind <- -1
  
  external.exon.ids.sorted <- exon.id.list[[i]][exon.gene.rank.ind.list[[i]]]
  internal.exon.ids.sorted <- exon.chr.id.list[[i]][exon.gene.rank.ind.list[[i]]]
  
  gene.left.exon.index.list  <- exon.gene.rank.index.list[[i]][["left"]]
  gene.right.exon.index.list <- exon.gene.rank.index.list[[i]][["right"]]

  if (length(gene.left.exon.index.list) > 0 && length(gene.right.exon.index.list) > 0) {
    for (j in 1:length(gene.left.exon.index.list)) {   
      
      left.exon.ind <- gene.left.exon.index.list[[j]][[1]][length(gene.left.exon.index.list[[j]][[1]])]

      if (left.exon.ind != old.left.exon.ind) {

        old.left.exon.ind <- left.exon.ind
        
        for (k in 1:length(gene.right.exon.index.list)) {
          
          right.exon.ind <- gene.right.exon.index.list[[k]][[1]][1]

          if (right.exon.ind != old.right.exon.ind) {
            
            old.right.exon.ind <- right.exon.ind

            if (noOverlap(internal.exon.ids.sorted[left.exon.ind], internal.exon.ids.sorted[right.exon.ind])) {
              
              if (left.exon.ind < right.exon.ind) {
                if (length(grep("/-$", internal.exon.ids.sorted[left.exon.ind])) > 0 && length(grep("/-$", internal.exon.ids.sorted[right.exon.ind])) > 0) {
                  junction.id[junction.ind] <- paste("junction", external.exon.ids.sorted[right.exon.ind], external.exon.ids.sorted[left.exon.ind], sep="|")
                  exon.id1   [junction.ind] <- external.exon.ids.sorted[right.exon.ind]
                  exon.id2   [junction.ind] <- external.exon.ids.sorted[left.exon.ind]
                } else {
                  junction.id[junction.ind] <- paste("junction", external.exon.ids.sorted[left.exon.ind], external.exon.ids.sorted[right.exon.ind], sep="|")
                  exon.id1   [junction.ind] <- external.exon.ids.sorted[left.exon.ind]
                  exon.id2   [junction.ind] <- external.exon.ids.sorted[right.exon.ind]
                }
                
                junction.ind <- junction.ind + 1
              }
            }
          }
        }
      }
    }
  }
}

###################################################################################
##
## Write junction exon map file
##
###################################################################################

if (! all((junction.id != "") == (exon.id1 != ""))) {
  print(paste("Non-empty junction ids and exon ids discrepancy."))
}

junctionExonMap.frame <- data.frame(junction.id[junction.id != ""], exon.id1[exon.id1 != ""], exon.id2[exon.id2 != ""])
junctionExonMap.file  <- Sys.getenv("JUNCTION_EXON_MAP_FILE")

print(paste("Writing", nrow(junctionExonMap.frame), "junction exon mapping records to file", junctionExonMap.file))
write.table(junctionExonMap.frame, junctionExonMap.file, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, na="NA")


