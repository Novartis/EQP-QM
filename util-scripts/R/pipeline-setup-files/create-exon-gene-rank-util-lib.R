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

## call with: source(file.path(Sys.getenv("HOME"), "ngs", "pipelines", "exon-pipeline", "R", "pipeline-setup-files", "create-exon-gene-rank-util-lib.R"))

###############################################################################
##
## Unify gene exon whose start or end position differs by at most t bases
##
###############################################################################

unifyExons <- function (exons, t, check=FALSE) {

  exon.num  <- length(exons)

  if (exon.num == 1) {
    return(list(unified.exons=exons, exon.ind=1))
  }
  
  exon.mat   <- convertToMatrix(exons, "/")
  exon.order <- order(exon.mat[,1], as.integer(exon.mat[,2]), as.integer(exon.mat[,3]))

  exon.chr    <- exon.mat[exon.order,1]
  exon.start  <- as.integer(exon.mat[exon.order,2])
  exon.end    <- as.integer(exon.mat[exon.order,3])
  exon.strand <- exon.mat[exon.order,4]

  exon.sort <- paste(exon.chr, exon.start, exon.end, exon.strand, sep="/")

  chromosome.equal.bool    <- exon.chr[1:(exon.num-1)] == exon.chr[2:exon.num]
  start.end.incorrect.bool <- exon.end[1:(exon.num-1)] > exon.start[2:exon.num] & chromosome.equal.bool

  start.diff <- abs(exon.start[1:(exon.num-1)] - exon.start[2:exon.num])
  end.diff   <- abs(exon.end  [1:(exon.num-1)] - exon.end  [2:exon.num])

  ## Check distribution of overlaps
  if (check) {
    overlap.exon.num <- array(0, dim=30)
    for (t in 1:30) {
      start.diff.bool <- start.diff <= t
      end.diff.bool   <- end.diff   <= t
      
      overlap.exon.num[t] <- sum(start.diff.bool & end.diff.bool & start.end.incorrect.bool)
      print(paste(t, sum(start.diff.bool & end.diff.bool & start.end.incorrect.bool)))
    }
    plot(overlap.exon.num)
  }

  diff.bool <- start.diff <= t & end.diff <= t
  ## different possibility: diff.bool  <- start.diff + end.diff <= t
  
  overlap.exon.num <- sum(diff.bool & start.end.incorrect.bool)
  print ("Number of overlapping exons:", overlap.exon.num)

  overlap.exon.ind     <- which(c(FALSE, diff.bool & start.end.incorrect.bool))
  overlap.all.exon.ind <- sort(c(overlap.exon.ind, overlap.exon.ind-1))

  ## We only unify exons on the same chromosome
  unified.exon.chr    <- exon.chr

  exon.ind.start      <- array(0, dim=exon.num)
  exon.ind.end        <- array(exon.num, dim=exon.num)
  unified.exon.start  <- array(0, dim=exon.num)
  unified.exon.end    <- array(0, dim=exon.num)
  k <- 1
  l <- 1
  ## For each sequence of consecutive indices in overlap.all.exon.ind set the end
  ## of the (sorted) exons to the end of the last exon indexed by the consecutive
  ## sequence in overlap.all.exon.ind. (Note that the start is already the minimum
  ## due to the sorting of the exons).
  for (i in 1:exon.num) {
    if (k <= length(overlap.exon.ind) && i == overlap.exon.ind[k]) {
      exon.ind.start[i]   <- l
      unified.exon.end[l] <- max(unified.exon.end[l], exon.end[i])
      k <- k + 1
    } else {
      if (i - 1 > l) {
        exon.ind.end[l] <- i - 1
        for (j in (l+1):(i-1)) {
          exon.ind.start[j]     <- exon.ind.start[l]
          exon.ind.end[j]       <- exon.ind.end[l]
          unified.exon.start[j] <- unified.exon.start[l]
          unified.exon.end[j]   <- unified.exon.end[l]
        }
      } else {
        exon.ind.end[i-1] <- i - 1
      }
      exon.ind.start[i] <- i
      exon.ind.end[i]   <- i
      unified.exon.start [i] <- exon.start[i]
      unified.exon.end   [i] <- exon.end[i]
      l <- i
    }
  }

  if (exon.num > l) {
    exon.ind.end[l] <- exon.num
    for (j in (l+1):(exon.num)) {
      exon.ind.start[j]     <- exon.ind.start[l]
      exon.ind.end[j]       <- exon.ind.end[l]
      unified.exon.start[j] <- unified.exon.start[l]
      unified.exon.end[j]   <- unified.exon.end[l]
    }
  } else {
    exon.ind.end[exon.num] <- exon.num
  }


  ## At the end of this procedure we should not have any unified exons within t bases anymore
  unified.exon.num <- length(unified.exon.chr)
  chromosome.unified.equal.bool    <- unified.exon.chr[1:(unified.exon.num-1)] == unified.exon.chr[2:unified.exon.num]
  start.end.unified.incorrect.bool <-
    unified.exon.end[1:(unified.exon.num-1)] > unified.exon.start[2:unified.exon.num] & chromosome.unified.equal.bool

  start.unified.diff <- abs(unified.exon.start[1:(unified.exon.num-1)] - unified.exon.start[2:unified.exon.num])
  end.unified.diff   <- abs(unified.exon.end  [1:(unified.exon.num-1)] - unified.exon.end  [2:unified.exon.num])

  unified.exon.equal.bool <- chromosome.unified.equal.bool & start.unified.diff == 0 & end.unified.diff == 0

  unified.diff.bool <- start.unified.diff <= t & end.unified.diff <= t
  overlap.unified.exon.num <- sum(unified.diff.bool & start.end.unified.incorrect.bool & ! unified.exon.equal.bool)

  if (overlap.unified.exon.num > 0) {
    print(paste("Unification for", t, "bases failed:", overlap.unified.exon.num, "unified exons still overlap."))
  }

  unified.exons <- paste(unified.exon.chr, unified.exon.start, unified.exon.end, sep="/")

  exon.inv.order <- 1:exon.num
  exon.inv.order[exon.order] <- 1:exon.num

  exon.diff.bool <- gsub("/[+-]", "", exons) != unified.exons[exon.inv.order]

  ## length(unique(unified.exons[exon.inv.order][exon.diff.bool])) ==
  ##  length(unique(gsub("/[+-]", "", exons))) - length(unique(unified.exons[exon.inv.order]))

  return(list(unified.exons=unified.exons[exon.inv.order], exon.ind.start=exon.ind.start[exon.inv.order],
              exon.ind.end=exon.ind.end[exon.inv.order]))

}

###############################################################################
##
## formatExonRank returns zero-filled rank (e.g. 003) or (for overlapping exons) 
##
###############################################################################

formatExonRank  <- function (rank, subrank, maxRank) {

  if (subrank == 0) {
    return (gsub(" ", "0", format(rank, width=ceiling(log10(maxRank + 1)))))
  }

  #return (paste(gsub(" ", "0", format(rank, width=ceiling(log10(maxRank + 1)))), gsub(" ", "0", format(subrank, width=ceiling(log10(min(maxRank / 3, 10) + 1)))), sep="/"))
  return (paste(gsub(" ", "0", format(rank, width=ceiling(log10(maxRank + 1)))), subrank, sep="/"))
  
}

###############################################################################
##
## formatExonSubrank returns zero-filled subrank (e.g. 05/09) if there are
## more than nine subranks for an exon
##
###############################################################################

formatExonSubrank  <- function (exonRanks) {

  majorRanks <- unlist(lapply(strsplit(exonRanks, "/"), function (x) { return (x[1]) })) ## Ugly, but avoids paranthesis imbalance of `[`, 1
  subRanks   <- unlist(lapply(strsplit(exonRanks, "/"), function (x) { return (x[2]) })) ## Ugly, but avoids paranthesis imbalance of `[`, 2
  
  subRank.len <- length (subRanks)
  if (subRank.len <= 9) {
    if (subRank.len == 1 && is.na(subRanks[1])) {
      returnRanks <- majorRanks
    } else {
      returnRanks <- paste(majorRanks, subRanks, sep="/")
    }
  } else {
    returnRanks <- paste(majorRanks, gsub(" ", "0", format(as.integer(subRanks), width=ceiling(log10(subRank.len + 1)))), sep="/")
  }

  return (returnRanks)

}


###############################################################################
##
## reformatExonRank returns a zero-filled subrank (e.g. 05/09) if there are
## more than nine subranks for an exon rank
##
###############################################################################

reformatExonRanks  <- function (exonRanks) {

  majorRanks <- unlist(lapply(strsplit(exonRanks, "/"), function (x) { return (x[1]) })) ## Ugly, but avoids paranthesis imbalance of `[`, 1
 
  return (as.character(unlist(tapply(exonRanks, majorRanks, formatExonSubrank))))
  
}


###############################################################################
##
## computeExonRank
##
## Compute the rank of an exon in the list of exons of a gene. Ranks are
## zero-padded strings which can be easily sorted lexicographically.
##
## Note that exons may occur in duplicate
##
## Note that sorting by chromosome/start/end is consistent for
## the opposite strand:
##
## 1.
##      ____________|-----|_________
##      _________|--------|_________   
##
## the lower interval goes first from left to right and, thus, second from
## right to left (which is correct as in this order it starts at the same
## position but is longer)
##
## 2.
##      _________|-----|____________
##      _________|--------|_________   
##
## the upper interval goes first from left to right and, thus, second from
## right to left (which is correct as in this order it starts after the start
## of the lower interval)
##
## In the above cases the gene rank of the exons is not increased but instead
## letters are appended to the rank, e.g. 4a, 4b, ... instead of 4, 5, ... etc.
##
###############################################################################

computeExonRank <- function (exons) {

  exon.unique     <- unique(exons)
  exon.unique.num <- length(exon.unique)
  
  exon.unique.mat   <- convertToMatrix(exon.unique, "/")

  geneId <- unique(exon.unique.mat[,1])

  if (length(geneId) > 1) {
    print(paste("Too many gene ids:", paste(geneId, collapse=", ")))
  }

  exon.unique.chr    <- exon.unique.mat[,2]
  exon.unique.strand <- exon.unique.mat[,5]

  exon.unique.strand.ind <- array(1, dim=length(exon.unique))
  exon.unique.strand.ind[exon.unique.mat[,5] == "-"] <- 2

  exon.unique.order <- order(exon.unique.chr, exon.unique.strand.ind, as.integer(exon.unique.mat[,3]), as.integer(exon.unique.mat[,4]))
  exon.unique.sort  <- exon.unique[exon.unique.order]

  ## Compute consecutive indices for which the chromosome and strand remain the same. Note that we can have a
  ## situation such 4x chrUn_random/-, 4x chrUn_random/+, and 4x chrUn_random/- (100041253, NM_001200055) which prevents 
  ## the application of functions such as unique or grep.
  exon.unique.chromosome.sort        <- as.character(exon.unique.mat[exon.unique.order,2])
  exon.unique.chromosome.start.sort  <- as.integer(exon.unique.mat[exon.unique.order,3])
  exon.unique.chromosome.end.sort    <- as.integer(exon.unique.mat[exon.unique.order,4])
  exon.unique.strand.sort            <- as.character(exon.unique.mat[exon.unique.order,5])
  exon.unique.chromosome.strand.sort <- paste(exon.unique.chromosome.sort, exon.unique.strand.sort, sep="/")
  
  exon.unique.chromosome.ind.list <- vector("list", 1)
  exon.unique.chromosome.ind.list[[1]] <- 1
  list.index <- 1
  exon.rank  <- 1
  exon.subrank <- 0
  overlapping.exon.index <- -1
  exon.maxRank <- length(exon.unique.chromosome.strand.sort)
  exon.unique.chromosome.rank.list <- vector("list", 1)
  exon.unique.chromosome.rank.list[[1]] <- formatExonRank(exon.rank, exon.subrank, exon.maxRank)
  if (length(exon.unique.chromosome.strand.sort) >= 2) {
    for (i in 2:length(exon.unique.chromosome.strand.sort)) {
      if (exon.unique.chromosome.strand.sort[i-1] == exon.unique.chromosome.strand.sort[i]) {
        exon.unique.chromosome.ind.list[[list.index]] <- c(exon.unique.chromosome.ind.list[[list.index]], i)
        if (exon.unique.chromosome.end.sort[i-1] < exon.unique.chromosome.start.sort[i] ||
            (overlapping.exon.index >= 1 && exon.unique.chromosome.end.sort[overlapping.exon.index] < exon.unique.chromosome.start.sort[i])) {
          exon.rank <- exon.rank + 1
          exon.subrank <- 0
          overlapping.exon.index <- -1
          exon.unique.chromosome.rank.list[[list.index]] <-
            c(exon.unique.chromosome.rank.list[[list.index]], formatExonRank(exon.rank, exon.subrank, exon.maxRank))
        } else {
          if (exon.subrank == 0) {
            exon.subrank <- exon.subrank + 1
            overlapping.exon.index <- i-1
            exon.unique.chromosome.rank.list[[list.index]][length(exon.unique.chromosome.rank.list[[list.index]])] <-
              formatExonRank(exon.rank, exon.subrank, exon.maxRank)
          }
          exon.subrank <- exon.subrank + 1
          exon.unique.chromosome.rank.list[[list.index]] <-
            c(exon.unique.chromosome.rank.list[[list.index]], formatExonRank(exon.rank, exon.subrank, exon.maxRank))
        }
      } else {
        list.index   <- list.index + 1
        exon.rank    <- exon.rank + 1
        exon.subrank <- 0
        exon.unique.chromosome.rank.list[[list.index]] <- formatExonRank(exon.rank, exon.subrank, exon.maxRank)
        exon.unique.chromosome.ind.list[[list.index]] <- i
      }

    }
  }

  ## We want to make the exon rank lexicographically comparable and, hence,
  ## need to reformat the subrank if there are more than nine subranks for a rank.
  gene.exon.rank          <- reformatExonRanks (unlist(exon.unique.chromosome.rank.list))
  names(exon.unique.sort) <- paste(geneId, "exon", gene.exon.rank, sep="-")
  exon.unique.ids.sort    <- names(exon.unique.sort)


  ## Invert the order the exons on the reverse strand of the gene 
  for (i in 1:length(exon.unique.chromosome.ind.list)) {
    exon.unique.ind <- exon.unique.chromosome.ind.list[[i]]
    exon.unique.ind.len <- length(exon.unique.ind)

    if (exon.unique.strand.sort[exon.unique.ind[1]] == "-" && exon.unique.ind.len > 1) {

      majorRanks <- unlist(lapply(strsplit(gene.exon.rank[exon.unique.ind], "/"), function (x) { return (x[1]) })) ## Ugly, but avoids paranthesis imbalance of `[`, 1
      subRanks   <- unlist(lapply(strsplit(gene.exon.rank[exon.unique.ind], "/"), function (x) { return (x[2]) })) ## Ugly, but avoids paranthesis imbalance of `[`, 2

      majorRank.multiplicity <- tapply(majorRanks, majorRanks, length)

      ## Rename the major ranks: Invert the names of the major ranks but not their multiplicities
      majorRanks.inv <- rep(names(majorRank.multiplicity)[length(majorRank.multiplicity):1], majorRank.multiplicity)

      ## Some major ranks do not have subranks
      subRanks.non_NA.bool <- ! is.na(subRanks)

      ## Paste the sub ranks to their new names
      exonRank.unsorted <- majorRanks.inv
      exonRank.unsorted[subRanks.non_NA.bool] <- paste(majorRanks.inv[subRanks.non_NA.bool], subRanks[subRanks.non_NA.bool], sep="/")

      ## Sort the subranks into the correct order (the major ranks are already in the correct order)
      gene.exon.rank[exon.unique.ind] <- sort (exonRank.unsorted, decreasing=TRUE)
      
    }    
  }

  ## exon.ind <- match(exons, exon.unique[exon.unique.order], nomatch=0)
  ##if (sum(exon.ind == 0) > 0) {
  ##  print(paste("Exons do not match after unique call for gene", geneId))
  ##}

  exon.unique.sort.ind <- match(names(exon.unique.sort), exon.unique.ids.sort, nomatch=0)
  stopifnot(sum(exon.unique.sort.ind == 0) == 0)
  
  names(exon.unique.sort) <- exon.unique.ids.sort

  exon.chromosome.ids.sort <- paste(exon.unique.mat[,2], exon.unique.mat[,3], exon.unique.mat[,4], exon.unique.mat[,5], sep="/")[exon.unique.order]

  return(list(gene.id=rep(geneId, length(exon.unique.ids.sort)),
              gene.exon.rank=gene.exon.rank,
              gene.exon.rank.id=paste(geneId, "exon", gene.exon.rank, sep="-"),
              exon.chr.id=exon.chromosome.ids.sort))
  
}


###############################################################################
##
## create exon gene rank map frame
##
###############################################################################

createExonGeneRankMapFrame <- function(gtf.chromosome.strand.ids, gtf.gene.ids, gtf.exon.ids) {

  print(paste("Generating a list with unique exons for each gene"))
  ## Compute a set of unique exons for each gene
  gene.exons.list   <- tapply (gtf.chromosome.strand.ids, gtf.gene.ids, unique)

  gene.id.levels         <- levels(as.factor(gtf.gene.ids))
  names(gene.exons.list) <- gene.id.levels
  gene.exon.num          <- sapply(gene.exons.list, length)

  gene.exons        <- unlist(gene.exons.list)
  names(gene.exons) <- rep(gene.id.levels, gene.exon.num)

  ## Assign ranks to the exons
  print(paste("Computing exon ranks"))
  gene.exon.rank.list     <- tapply(paste(gtf.gene.ids, gtf.chromosome.strand.ids, sep="/"), gtf.gene.ids, computeExonRank)
  gene.exon.rank.list.len <- sapply(gene.exon.rank.list, function (x) { return(length(x[[1]])) })
  
  ##sum(gene.exon.rank.list.len != gene.exon.num)

  gene.exon.rank     <- unlist(lapply(gene.exon.rank.list, function (x) { return (x[["gene.exon.rank"]]) }))
  exon.gene.id       <- unlist(lapply(gene.exon.rank.list, function (x) { return (x[["gene.id"]]) }))
  exon.chromosome.id <- unlist(lapply(gene.exon.rank.list, function (x) { return (x[["exon.chr.id"]]) }))
  gene.exon.rank.id  <- unlist(lapply(gene.exon.rank.list, function (x) { return (x[["gene.exon.rank.id"]]) }))

  exon.chromosome.id.ind <- match(paste(exon.gene.id, exon.chromosome.id, sep="/"), paste(gtf.gene.ids, gtf.chromosome.strand.ids, sep="/"), nomatch=0)

  exon.map.frame <- data.frame(exon.chromosome.id, gtf.exon.ids[exon.chromosome.id.ind], exon.gene.id, gene.exon.rank)
  names(exon.map.frame) <- c("Exon Chromosome Id", "Exon Id", "Gene Id", "Exon gene rank")

  return (exon.map.frame)
  
}
