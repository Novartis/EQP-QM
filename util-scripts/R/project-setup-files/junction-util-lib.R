#call with: source(file.path(Sys.getenv("HOME"), "ngs", "pipelines", "exon-pipeline", "R", "project-setup-files", "junction-util-lib.R"))

###############################################################################
##
## transposeElement
##
###############################################################################

transposeElement <- function (x) {

  if (length(x) == 0) {
    return(NULL)
  }
  
  if (is.matrix(x)) {
    return(t(x))
  }

  if (is.vector(x)) {
    x.mat <- matrix(x, ncol=1)
    return (t(x.mat))
  }
  
  return(x)

}

###############################################################################
##
## matrixDim
##
###############################################################################

matrixDim <- function (x, d=1) {

  if (length(x) == 0) {
    return(0)
  }
  
  if (is.matrix(x)) {
    return(dim(x)[d])
  }

  return(1)

}

###############################################################################
##
## maxElementLength
##
###############################################################################

maxElementLength <- function (x) {

  if (length(x) == 0) {
    return (0)
  }

  if (is.list(x)) {
    return(max(unlist(lapply(x, maxElementLength))))
  }
  
  return(length(x))

}

###############################################################################
##
## return the sum of the elements indexed by the first argument
##
###############################################################################

indexSum <- function(indices, x, restrict=0) {
  
  if (is.na(restrict)) {
    return(sum(x[indices][-length(x[indices])]))
  }
  
  if (restrict == 0) {
    return(sum(x[indices]))
  }

  return(sum(x[indices][restrict]))
  
}


###############################################################################
##
## geneExonLengths
##
###############################################################################

geneExonLengths <- function(exons) {

  exon.mat   <- convertToMatrix(exons, "/")

  strand.ind <- grep("^[+-]$", exon.mat[1,])[1]
  chr.ind    <- strand.ind - 3
  start.ind  <- strand.ind - 2
  end.ind    <- strand.ind - 1

  return(as.integer(exon.mat[,end.ind]) - as.integer(exon.mat[,start.ind]) + 1)
  
}


###############################################################################
##
## exonOverlap.bool
##
## Note that we define two exons to overlap even if they abut (since we are
## interested in true junctions):
##
##           exon.end[i] (>)= exon.start[i] - exon.sep
##
## that is, the separation has to be at least exon.sep bp.
##
###############################################################################

exonOverlap.bool <- function(exon.start, exon.end, i, j, exon.sep) {

  if (exon.start[i] <= exon.start[j]) {
    return (exon.end[i] >= exon.start[j] - exon.sep)
  }

  return (exon.end[j] >= exon.start[i] - exon.sep)
  
}


###############################################################################
##
## geneExonLengths
##
###############################################################################

numChromosomeStrandGene <- function(exons) {

  exon.mat   <- convertToMatrix(exons, "/")

  strand.ind <- grep("^[+-]$", exon.mat[1,])[1]
  chr.ind    <- strand.ind - 3

  return(length(unique(paste(exon.mat[,chr.ind], exon.mat[,strand.ind], sep="/"))))
}


###############################################################################
##
## getBedEntryFieldColumn
##
###############################################################################

getBedEntryFieldColumn <- function(bedEntryField, column, summarizeFunction=NULL) {

  uniqueField <- unique(bedEntryField[,column])
  if (length(uniqueField) != 1) {
    if (is.null(summarizeFunction)) {
      print(paste("Problem with column ", column))
      print(bedEntryField)
    } else {
      uniqueField <- summarizeFunction (as.numeric (uniqueField))
    }
  }
  return(uniqueField)
  
}


###############################################################################
##
## getBedEntryField
##
###############################################################################

getBedEntryFieldList <- function(bedEntry, fieldName, column) {

  return(sapply(bedEntry[[fieldName]], getBedEntryFieldColumn, column))
  
}


###############################################################################
##
## getBedEntryFieldInteger
##
###############################################################################

getBedEntryFieldListSummarize <- function(bedEntry, fieldName, column, summarizeFunction) {

  return(sapply(bedEntry[[fieldName]], getBedEntryFieldColumn, column, summarizeFunction))
  
}


###############################################################################
##
## getSortedExonsInd
##
###############################################################################

getSortedExonsInd <- function (exons) {

  exon.mat   <- convertToMatrix(exons, "/")
  
  strand.ind <- grep("^[+-]$", exon.mat[1,])[1]
  chr.ind    <- strand.ind - 3
  start.ind  <- strand.ind - 2
  end.ind    <- strand.ind - 1

  exon.order <- order(exon.mat[,chr.ind], as.integer(exon.mat[,start.ind]), as.integer(exon.mat[,end.ind]))

  return (exon.order)

}


###############################################################################
##
## getGeneId
##
###############################################################################

getGeneId <- function (exon.ids.mat) {
  
  exon.ids.len <- ncol(exon.ids.mat)

  if (exon.ids.len > 3) {
    gene.id.list <- apply(exon.ids.mat[,1:(exon.ids.len-2)], 2, unique)
    if (max(unlist(lapply(gene.id.list, length))) > 1) {
      print(paste("Too many gene ids for", paste(exon.ids.mat, collapse=", ")))
      stopifnot(FALSE)
    }
    gene.id <- paste(unlist(gene.id.list), collapse="-")
  } else {
    gene.id  <- unique(exon.ids.mat[,1])
  }

  return (gene.id)

}


###############################################################################
##
## findLeftExonIndices
##
## Returns a list of exon indices vectors (leftIndices) for which the chromosome and 
## strand are the same as in entry i. The entries of the vectors are at most
## than i, their length is at most maxExon.num, and the sum of the
## lengths of the corresponding exons is larger than radius (unless i is too
## close to the beginning of the list).
##
###############################################################################

findLeftExonIndices <- function(exon.len, i, radius, exon.chr, exon.strand, exon.start, exon.end, exon.sep, maxExon.num=-1, gene.id=NA) {

  if (maxExon.num == 0) {
    return (NULL)
  }
  
  potential.exon.ind <- which(exon.chr == exon.chr[i] & exon.strand == exon.strand[i])

  ## left.potential.exon.ind is initialized with the index of i in potential.exon.ind.
  ## This becomes the last entry of a list of indices into potential.exon.ind
  left.potential.exon.ind <- which(potential.exon.ind == i)

  ## Find a second non-overlapping exon (if possible)
  first.exon.ind <- potential.exon.ind[left.potential.exon.ind[1]] 
  non.overlapping.exon.ind <- left.potential.exon.ind - 1
  while (non.overlapping.exon.ind >= 1 &&
         exonOverlap.bool (exon.start, exon.end, first.exon.ind, potential.exon.ind[non.overlapping.exon.ind], exon.sep)) {
    non.overlapping.exon.ind <- non.overlapping.exon.ind - 1
  }

  leftIndices <- vector("list", 1)
  ## Check if one entry suffices
  if (exon.len[i] >= radius || left.potential.exon.ind == 1 || maxExon.num == 1 || non.overlapping.exon.ind < 1) {

    j <- non.overlapping.exon.ind
    while (sum(exon.len[potential.exon.ind[left.potential.exon.ind]]) < radius && j >= 1) {
      if (! exonOverlap.bool (exon.start, exon.end, left.potential.exon.ind[1], j, exon.sep)) {
        left.potential.exon.ind  <- c(j, left.potential.exon.ind)
      }
      j <- j - 1
    }
        
    leftIndices[[1]] <- potential.exon.ind[left.potential.exon.ind]
    return(leftIndices)
    
  }

  ## Create possible sets of exon indices to the left of the junction of size at least 2
  k <- 1
  left.potential.exon.ind <- c(non.overlapping.exon.ind, left.potential.exon.ind)
  sumExonLength <- sum(exon.len[potential.exon.ind[left.potential.exon.ind]])
  while (length(left.potential.exon.ind) > 1) {
    ##    if (sum(exon.len[potential.exon.ind[left.potential.exon.ind]]) != sumExonLength) {
    ##      print(paste("Sum of exon length incorrect in findLeftExonIndices for gene", gene.id, ":",
    ##                  sum(exon.len[potential.exon.ind[left.potential.exon.ind]]), " != ", sumExonLength))
    ##      stopifnot(sum(exon.len[potential.exon.ind[left.potential.exon.ind]]) == sumExonLength)
    ##    }
    
    ## Note that this condition is satisfied for the first iteration by design but maybe violated in the iterations that follow
    if (! exonOverlap.bool (exon.start, exon.end, potential.exon.ind[left.potential.exon.ind[1]], potential.exon.ind[left.potential.exon.ind[2]], exon.sep)) {
      j <- left.potential.exon.ind[1] - 1
      while (sumExonLength < radius && j >= 1 && (maxExon.num == -1 || length(left.potential.exon.ind) < maxExon.num)) {
        if (! exonOverlap.bool (exon.start, exon.end, potential.exon.ind[left.potential.exon.ind[1]], potential.exon.ind[j], exon.sep)) {
          left.potential.exon.ind <- c(j, left.potential.exon.ind)
          sumExonLength <- sumExonLength + exon.len[potential.exon.ind[left.potential.exon.ind[1]]]
        }
        j <- j - 1
      }

      ## If we have reached the maximally allowed number of elements (maxExon.num) but the size of the
      ## exon collection is smaller than radius, then we just fill the exon collection with consecutive
      ## elements without influencing the enumeration of the remaining possibilities.
      left.potential.exon.ind.extended <- left.potential.exon.ind
      sumExonLength.extended <- sumExonLength
      j <- left.potential.exon.ind.extended[1] - 1
      while (sumExonLength.extended < radius && 1 <= j) {
        if (! exonOverlap.bool (exon.start, exon.end, potential.exon.ind[left.potential.exon.ind[1]], potential.exon.ind[j], exon.sep)) {
          left.potential.exon.ind.extended <- c(j, left.potential.exon.ind.extended)
          sumExonLength.extended <- sumExonLength.extended + exon.len[potential.exon.ind[j]]
        }
        j <- j - 1        
      }
      
      leftIndices[[k]] <- potential.exon.ind[left.potential.exon.ind.extended]
      exons <- paste(exon.chr[leftIndices[[k]]], exon.strand[leftIndices[[k]]], exon.start[leftIndices[[k]]], exon.end[leftIndices[[k]]], sep="/")
      k <- k + 1
    }

    ## Go the next possibility by descreasing the first element
    if (left.potential.exon.ind[1] == 1) {
      ## Remove the first element
      sumExonLength <- sumExonLength - exon.len[potential.exon.ind[left.potential.exon.ind[1]]]
      # print(paste("sumExonLength after removing the first element of length", exon.len[left.potential.exon.ind[1]], " - index:", ":", sumExonLength))
      left.potential.exon.ind <- left.potential.exon.ind[-1]
    }

    sumExonLength <- (sumExonLength - exon.len[potential.exon.ind[left.potential.exon.ind[1]]] +
                      exon.len[potential.exon.ind[left.potential.exon.ind[1] - 1]])
    left.potential.exon.ind[1] <- left.potential.exon.ind[1] - 1
  }
  
  return (leftIndices)
  
}


###############################################################################
##
## findRightExonIndices
##
## Returns a list of exon indices vectors (rightIndices) for which the chromosome and 
## strand are the same as in entry i. The entries of the vectors are larger
## than i, their length is at most maxExon.num, and the sum of the
## lengths of the corresponding exons is larger than radius (unless i is too
## close to the end of the list).
##
###############################################################################

findRightExonIndices <- function(exon.len, i, radius, exon.chr, exon.strand, exon.start, exon.end, exon.sep, maxExon.num=-1, gene.id=NA) {

  if (maxExon.num == 0) {
    return (NULL)
  }

  potential.exon.ind <- which(exon.chr == exon.chr[i] & exon.strand == exon.strand[i])
  
  ## right.potential.exon.ind is initialized with the index of i in potential.exon.ind. 
  ## It becomes is the first entry of a list of indices into potential.exon.ind
  right.potential.exon.ind <- which(potential.exon.ind == i)
  potential.exon.num <- length(potential.exon.ind)

  ## Find a second non-overlapping exon (if possible)
  first.exon.ind <- potential.exon.ind[right.potential.exon.ind[1]] 
  non.overlapping.exon.ind <- right.potential.exon.ind + 1
  while (non.overlapping.exon.ind <= potential.exon.num &&
         exonOverlap.bool (exon.start, exon.end, first.exon.ind, potential.exon.ind[non.overlapping.exon.ind], exon.sep)) {
    non.overlapping.exon.ind <- non.overlapping.exon.ind + 1
  }

  rightIndices <- vector("list", 1)
  ## Check if one entry suffices
  if (exon.len[i] >= radius || right.potential.exon.ind == potential.exon.num  || maxExon.num == 1 || non.overlapping.exon.ind > potential.exon.num) {

    right.potential.exon.ind.num <- length(right.potential.exon.ind)
    j <- non.overlapping.exon.ind
    while (sum(exon.len[potential.exon.ind[right.potential.exon.ind]]) < radius && j <= potential.exon.num) {
      if (! exonOverlap.bool (exon.start, exon.end, right.potential.exon.ind[right.potential.exon.ind.num], j, exon.sep)) {
        right.potential.exon.ind     <- c(right.potential.exon.ind, j)
        right.potential.exon.ind.num <- right.potential.exon.ind.num + 1
      }
      j <- j + 1
    }

    rightIndices[[1]] <- potential.exon.ind[right.potential.exon.ind]
    return(rightIndices)
    
  }

  ## Create possible sets of exon indices to the right of the junction of size at least 2
  k <- 1
  right.potential.exon.ind     <- c(right.potential.exon.ind, non.overlapping.exon.ind)
  right.potential.exon.ind.num <- length(right.potential.exon.ind)
  sumExonLength = sum(exon.len[potential.exon.ind[right.potential.exon.ind]])
  while (length(right.potential.exon.ind) > 1) {
    ##    if (sum(exon.len[potential.exon.ind[right.potential.exon.ind]]) != sumExonLength) {
    ##      print(paste("Sum of exon length incorrect in findRightExonIndices for gene", gene.id, ":",
    ##                  sum(exon.len[potential.exon.ind[right.potential.exon.ind]]), " != ", sumExonLength))
    ##      stopifnot(sum(exon.len[potential.exon.ind[right.potential.exon.ind]]) == sumExonLength)
    ##    }

    ## Note that this condition is satisfied for the first iteration by design but maybe violated in the iterations that follow
    if (! exonOverlap.bool (exon.start, exon.end, potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num-1]],
                            potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]], exon.sep)) {
      j <- right.potential.exon.ind[right.potential.exon.ind.num] + 1
      while (sumExonLength < radius && j <= potential.exon.num && (maxExon.num == -1 || length(right.potential.exon.ind) < maxExon.num)) {
        if (! exonOverlap.bool (exon.start, exon.end, potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]], potential.exon.ind[j], exon.sep)) {
          right.potential.exon.ind     <- c(right.potential.exon.ind, j)
          right.potential.exon.ind.num <- right.potential.exon.ind.num + 1
          sumExonLength <- sumExonLength + exon.len[potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]]]
        }
        j <- j + 1
      }
      
      ## If we have reached the maximally allowed number of elements (maxExon.num) but the size of the
      ## exon collection is smaller than radius, then we just fill the exon collection with consecutive
      ## elements without influencing the enumeration of the remaining possibilities.
      right.potential.exon.ind.extended     <- right.potential.exon.ind
      right.potential.exon.ind.extended.num <- right.potential.exon.ind.num
      sumExonLength.extended <- sumExonLength
      j <- right.potential.exon.ind[right.potential.exon.ind.num] + 1
      while (sumExonLength.extended < radius && j <= potential.exon.num) {
        if (! exonOverlap.bool(exon.start, exon.end, potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]], 
                               potential.exon.ind[j], exon.sep)) {
          right.potential.exon.ind.extended     <- c(right.potential.exon.ind.extended, j)
          right.potential.exon.ind.extended.num <- right.potential.exon.ind.extended.num + 1
          sumExonLength.extended <- sumExonLength.extended + exon.len[potential.exon.ind[j]]
        }
        j <-  j + 1
      }
      
      rightIndices[[k]] <- potential.exon.ind[right.potential.exon.ind.extended]
      k <- k + 1
      
    }
      
    ## Go the next possibility by increasing the last element
    if (right.potential.exon.ind[right.potential.exon.ind.num] == potential.exon.num) {
      ## Remove the last element
      sumExonLength <- sumExonLength - exon.len[potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]]]
      right.potential.exon.ind     <- right.potential.exon.ind[-right.potential.exon.ind.num]
      right.potential.exon.ind.num <- right.potential.exon.ind.num - 1
    }
    
    sumExonLength <- (sumExonLength - exon.len[potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num]]] +
                      exon.len[potential.exon.ind[right.potential.exon.ind[right.potential.exon.ind.num] + 1]])
    right.potential.exon.ind[right.potential.exon.ind.num] <- right.potential.exon.ind[right.potential.exon.ind.num] + 1
      
  }

  return (rightIndices)
  
}


###############################################################################
##
## findExonIndices
##
## Returns a list with the left and the right exon indices for each element of 
## exons. So we have the following structure. Let the number of exons be n.
##
##  -> left                     -> right
##     (n-1 element list,          (n-1 element list,
##      one for each exon           one for each exon
##      except the last)            except the first)
##     Element i (corresponds      Element i (corresponds
##     to exon i):                 to exon i+1)
##       -> list of vectors         -> list of vectors
##          each ends with i           each starts with i + 1
##
###############################################################################

findExonIndices <- function(exons, maxExon.num, radius, exon.sep, printGeneId=FALSE) {

  exon.num <- length(exons)

  if (exon.num == 1) {
    return(c())
  }

  exon.mat   <- convertToMatrix(exons, "/")
  strand.ind <- grep("^[+-]$", exon.mat[1,])[1]
  chr.ind    <- strand.ind - 3
  start.ind  <- strand.ind - 2
  end.ind    <- strand.ind - 1

  exon.order   <- order(exon.mat[,chr.ind], as.integer(exon.mat[,start.ind]), as.integer(exon.mat[,end.ind]))
  exon.chr     <- exon.mat[exon.order,chr.ind]
  exon.strand  <- exon.mat[exon.order,strand.ind]

  stopifnot(length(names(exons)) > 0)
  exon.ids.mat <- convertToMatrix(names(exons)[exon.order], "-")
  gene.id      <- getGeneId(exon.ids.mat)

  if (printGeneId) {
    print (gene.id)
  }

  if (length(gene.id) != 1) {
    print(paste("Too many (or few) gene ids for", paste(exons, collapse=", ")))
  }

  exon.start <- as.integer(exon.mat[exon.order,start.ind])
  exon.end   <- as.integer(exon.mat[exon.order,end.ind])

  exon.len <- exon.end - exon.start + 1

  ## Process the exon ids
  ## exon.ids.mat <- convertToMatrix(names(exons)[exon.order], "-")

  leftIndices  <- vector("list", length(exon.len) - 1)
  for (i in 1:(length(exon.len) - 1)) {
    ## print(paste("Left:", i))
    leftIndices[[i]] <- findLeftExonIndices(exon.len, i, radius, exon.chr, exon.strand, exon.start, exon.end, exon.sep, maxExon.num, gene.id)
  }
  
  rightIndices <- vector("list", length(exon.len) - 1)
  for (i in 2:length(exon.len)) {
    ## print(paste("Right:", i))
    rightIndices[[i-1]] <- findRightExonIndices(exon.len, i, radius, exon.chr, exon.strand, exon.start, exon.end, exon.sep, maxExon.num, gene.id)
  }

  return(list(left=leftIndices, right=rightIndices))
}


###############################################################################
##
## getExonCombinationNum
## compute the number of combinations for each exon on the side <side>
##
###############################################################################

getExonCombinationNum <- function (x, side="left") {

  if (length(x) == 0) {
    return (0)
  }
  
  left.len <- sapply(x[[side]], length)
  
  return (left.len)
  
}


###############################################################################
##
## getGeneCombinationNum
##
###############################################################################

getGeneCombinationNum <- function (x) {
  
  leftNumCombinations  <- getExonCombinationNum(x, "left")
  rightNumCombinations <- getExonCombinationNum(x, "right")

  if (length(leftNumCombinations) == 0 || length(rightNumCombinations) == 0) {
    return (0)
  }

  ## We can combine the left entries of exon i with the right entries of exons i+1, i+2, ...
  ## The right entries of exons i+1, i+2, ... are contained in rightNumCombinations (r)
  ## i, i+1, i+2, ..., and, we get
  ##
  ##     sum_{j>=i} r_j  =  sum_{j} r_j - sum_{j<=i} r_j + r_i
  
  sumRightCombinations <- sum(rightNumCombinations) - cumsum(rightNumCombinations) + rightNumCombinations

  return (sum(as.numeric(leftNumCombinations * sumRightCombinations)))
     
}

###############################################################################
##
## getGeneEntryNum
##
###############################################################################

getGeneEntryNum <- function (x) {
  
  return(sum(getExonCombinationNum(x, "left")) + sum(getExonCombinationNum(x, "right")))
     
}


###############################################################################
##
## computeBedEntries
##
###############################################################################

computeBedEntries <- function (index.list, exons, radius) {

  if (length(index.list) == 0) {
    return (c())
  }

  if (0 == 1) {
    index.list <- gene.exon.index.list[[i]]
    exons <- gene.exon.list[[i]]
  }
  
  leftGeneIndices  <- index.list[["left"]]
  rightGeneIndices <- index.list[["right"]]

  exon.mat   <- convertToMatrix(exons, "/")
  strand.ind <- grep("^[+-]$", exon.mat[1,])[1]
  chr.ind    <- strand.ind - 3
  start.ind  <- strand.ind - 2
  end.ind    <- strand.ind - 1

  exon.order <- order(exon.mat[,chr.ind], as.integer(exon.mat[,start.ind]), as.integer(exon.mat[,end.ind]))
  exons.sort <- exons[exon.order]

  exon.ids.mat <- convertToMatrix(names(exons)[exon.order], "-")
  exon.ids.len <- ncol(exon.ids.mat)
  exon.ids     <- exon.ids.mat[,exon.ids.len]
  gene.id      <- getGeneId (exon.ids.mat)

  if (length(gene.id) != 1) {
    print(paste("Too many (or few) gene ids for", paste(exons, collapse=", ")))
    stopifnot(FALSE)
  }

  exon.chr    <- exon.mat[exon.order,chr.ind]
  exon.start  <- as.integer(exon.mat[exon.order,start.ind])
  exon.end    <- as.integer(exon.mat[exon.order,end.ind])
  exon.strand <- exon.mat[exon.order,strand.ind]
  exon.len    <- exon.end - exon.start + 1

  ## Take out one level in the hierarchy of the lists for the left half
  ## Compute for each exon the length of the list of index vectors
  leftGeneIndices.len  <- sum(sapply(leftGeneIndices, length))
  leftGeneIndices.list <- vector("list", leftGeneIndices.len)
  k <- 1
  for (i in 1:length(leftGeneIndices)) {
    for (j in 1:length(leftGeneIndices[[i]])) {
      leftGeneIndices.list[[k]] <- leftGeneIndices[[i]][[j]]
      k <- k+1
    }
  }

  ## Process list of exon indices vectors
  leftJunctionExon   <- array(0, dim=leftGeneIndices.len)
  leftJunctionPrefix <- array("", dim=leftGeneIndices.len)
  leftBedEntries     <- vector("list", leftGeneIndices.len)
  for (i in 1:leftGeneIndices.len) {
    leftExonIndices     <- leftGeneIndices.list[[i]]
    leftExonIndices.num <- length(leftExonIndices)
    
    leftJunctionExon[i] <- exon.ids[leftExonIndices[leftExonIndices.num]]
      
    leftJunctionPrefix[i] <- ""
    if (leftExonIndices.num > 1) {
      leftJunctionPrefix[i] <- paste(exon.ids[leftExonIndices[1:(leftExonIndices.num-1)]], collapse="-")
    }

    if (radius > 0) {
      curRadius <- radius
    } else {
      curRadius <- .Machine$integer.max
    }
    leftBedEntries[[i]] <- array("", dim=c(leftExonIndices.num, 11))
    for (j in leftExonIndices.num:1) {
      
      leftBedEntries[[i]][j, 1] <- exon.chr[leftExonIndices[j]]
      leftBedEntries[[i]][j, 2] <- max(exon.start[leftExonIndices[j]], exon.end[leftExonIndices[j]] - curRadius + 1) - 1
      leftBedEntries[[i]][j, 3] <- exon.end[leftExonIndices[j]]
      leftBedEntries[[i]][j, 5] <- "0"
      leftBedEntries[[i]][j, 6] <- exon.strand[leftExonIndices[j]]
      leftBedEntries[[i]][j, 7] <- exons.sort[leftExonIndices[j]]
      
      curRadius <- max(0, curRadius - exon.len[leftExonIndices[j]])

    }
    names(leftBedEntries)[i] <- leftJunctionExon[i]


    curJunctionCoord <- 1
    ## Since we are using bed coordinates, we do not need to add 1 in the length calculation
    exonJunction.len <- as.integer(leftBedEntries[[i]][, 3]) - as.integer(leftBedEntries[[i]][, 2])
    for (j in 1:leftExonIndices.num) {

      ## Compute the difference between the actual exon length and the length of
      ## the exon on the junction
      exonLen.diff <- exon.len[leftExonIndices[j]] - exonJunction.len[j]
      
      leftBedEntries[[i]][j, 8]  <- curJunctionCoord - exonLen.diff - 1
      leftBedEntries[[i]][j, 9]  <- curJunctionCoord + exonJunction.len[j] - 1
      leftBedEntries[[i]][j, 10] <- exonJunction.len[j]
      curJunctionCoord <- curJunctionCoord + exonJunction.len[j]

    }
  }

  
  ## Take out one level in the hierarchy of the lists for the right half as well
  rightGeneIndices.len <- sum(sapply(rightGeneIndices, length))
  rightGeneIndices.list <- vector("list", rightGeneIndices.len)
  k <- 1
  for (i in 1:length(rightGeneIndices)) {
    for (j in 1:length(rightGeneIndices[[i]])) {
      rightGeneIndices.list[[k]] <- rightGeneIndices[[i]][[j]]
      k <- k+1
    }
  }

  ## Process list of exon indices vectors
  rightJunctionExon   <- array(0, dim=rightGeneIndices.len)
  rightJunctionSuffix <- array("", dim=rightGeneIndices.len)
  rightBedEntries     <- vector("list", rightGeneIndices.len)
  for (i in 1:rightGeneIndices.len) {
    rightExonIndices     <- rightGeneIndices.list[[i]]
    rightExonIndices.num <- length(rightExonIndices)
    
    rightJunctionExon[i] <- exon.ids[rightExonIndices[1]]
      
    rightJunctionSuffix[i] <- ""
    if (rightExonIndices.num > 1) {
      rightJunctionSuffix[i] <- paste(exon.ids[rightExonIndices[2:rightExonIndices.num]], collapse="-")
    }
      
    if (radius > 0) {
      curRadius <- radius
    } else {
      curRadius <- .Machine$integer.max/2
    }
    rightBedEntries[[i]] <- array("", dim=c(rightExonIndices.num, 11))
    for (j in 1:rightExonIndices.num) {
      
      rightBedEntries[[i]][j, 1] <- exon.chr[rightExonIndices[j]]
      rightBedEntries[[i]][j, 2] <- exon.start[rightExonIndices[j]] - 1
      rightBedEntries[[i]][j, 3] <- min(exon.start[rightExonIndices[j]] + curRadius - 1, exon.end[rightExonIndices[j]])
      rightBedEntries[[i]][j, 5] <- "0"
      rightBedEntries[[i]][j, 6] <- exon.strand[rightExonIndices[j]]
      rightBedEntries[[i]][j, 7] <- exons.sort[rightExonIndices[j]]

      curRadius <- curRadius - exon.len[rightExonIndices[j]]

    }
    names(rightBedEntries)[i] <- rightJunctionExon[i]

    ## Note that junction coordinate for the right half has to be
    ## increased by the length of the left half (usually radius but
    ## there may be smaller left halves for small exons)
    curJunctionCoord <- 1
    ## Since we are using bed coordinates, we do not need to add 1 in the length calculation
    exonJunction.len <- as.integer(rightBedEntries[[i]][,3]) - as.integer(rightBedEntries[[i]][,2])
    for (j in 1:rightExonIndices.num) {
      
      ## Compute the difference between the actual exon length and the length of
      ## the exon on the junction
      exonLen.diff <- exon.len[rightExonIndices[j]] - exonJunction.len[j]
      
      rightBedEntries[[i]][j,  8] <- curJunctionCoord - 1
      rightBedEntries[[i]][j,  9] <- curJunctionCoord + exonJunction.len[j] - 1 + exonLen.diff
      rightBedEntries[[i]][j, 10] <- exonJunction.len[j]
      curJunctionCoord <- curJunctionCoord + exonJunction.len[j]

    }
  }

  return (list(leftGeneIndicesLen=leftGeneIndices.len, leftJunctionExon=leftJunctionExon, leftJunctionPrefix=leftJunctionPrefix,
               leftBedEntries=leftBedEntries, rightGeneIndicesLen=rightGeneIndices.len, rightJunctionExon=rightJunctionExon,
               rightJunctionSuffix=rightJunctionSuffix, rightBedEntries=rightBedEntries))
}


###############################################################################
##
## combineBedEntries
##
###############################################################################

combineBedEntries <- function (bedEntry.list, gene.id, exon.sep) {

  if (0 == 1) {
    bedEntry.list <- list(exon.gene.rank.bed.entry.list[[n]])
    gene.id <- gene.id.levels[n]
  }

  ## Compute chromosome, strand, start, and end for the bed entries
  ## Note there can be more than one BED entry for one exon
  leftChromosome.list  <- lapply(bedEntry.list, getBedEntryFieldList, "leftBedEntries", 1)
  leftStart.list       <- lapply(bedEntry.list, getBedEntryFieldListSummarize, "leftBedEntries", 2, min)
  leftEnd.list         <- lapply(bedEntry.list, getBedEntryFieldListSummarize, "leftBedEntries", 3, max)
  leftStrand.list      <- lapply(bedEntry.list, getBedEntryFieldList, "leftBedEntries", 6)
  rightChromosome.list <- lapply(bedEntry.list, getBedEntryFieldList, "rightBedEntries", 1)
  rightStart.list      <- lapply(bedEntry.list, getBedEntryFieldListSummarize, "rightBedEntries", 2, min)
  rightEnd.list        <- lapply(bedEntry.list, getBedEntryFieldListSummarize, "rightBedEntries", 3, max)
  rightStrand.list     <- lapply(bedEntry.list, getBedEntryFieldList, "rightBedEntries", 6)

  ## Count the number of bed records
  bedRecord.num  <- 0
  junctionId.num <- 0
  entry.num <- 0
  for (m in 1:length(bedEntry.list)) {

    if (m %% 1000 == 0) {
      print (paste(m, "genes of", length(bedEntry.list),"genes processed."))
    }

    if (length(bedEntry.list[[m]]) > 0) {
    
      leftGeneIndices.len  <- bedEntry.list[[m]][["leftGeneIndicesLen"]]
      leftJunctionExon     <- bedEntry.list[[m]][["leftJunctionExon"]]
      leftBedEntries       <- bedEntry.list[[m]][["leftBedEntries"]]
      
      rightGeneIndices.len <- bedEntry.list[[m]][["rightGeneIndicesLen"]]
      rightJunctionExon    <- bedEntry.list[[m]][["rightJunctionExon"]]
      rightBedEntries      <- bedEntry.list[[m]][["rightBedEntries"]]
      
      for (i in 1:leftGeneIndices.len) {

        chromosomeStrand.bool <-
          rightStrand.list[[m]] == leftStrand.list[[m]][i] & rightChromosome.list[[m]] == leftChromosome.list[[m]][i]        
        if (leftStrand.list[[m]][i] == "+") {
          compatibleRightJunction.bool <- chromosomeStrand.bool & rightJunctionExon > leftJunctionExon[i]
        } else {
          compatibleRightJunction.bool <- chromosomeStrand.bool & rightJunctionExon < leftJunctionExon[i] 
        }
        
        entry.num <- entry.num + sum(compatibleRightJunction.bool)
        for (j in which(compatibleRightJunction.bool)) {
          bedEntry.num   <- dim(leftBedEntries[[i]])[1] + dim(rightBedEntries[[j]])[1]
          junctionId.num <- junctionId.num + 1
          bedRecord.num  <- bedRecord.num + bedEntry.num
        }
      }
    }
  }
  
  bedCoordinates <- array("", dim=c(bedRecord.num, 11))
  junctionId     <- array("", dim=junctionId.num)
  k <- 1
  l <- 1
  count               <- 0
  skipped.entries.num     <- 0
  skipped.entries.ind <- c()
  for (m in 1:length(bedEntry.list)) {

    if (length(bedEntry.list[[m]]) > 0) {
      leftGeneIndices.len  <- bedEntry.list[[m]][["leftGeneIndicesLen"]]
      leftJunctionExon     <- bedEntry.list[[m]][["leftJunctionExon"]]
      leftJunctionPrefix   <- bedEntry.list[[m]][["leftJunctionPrefix"]]
      leftBedEntries       <- bedEntry.list[[m]][["leftBedEntries"]]
      
      rightGeneIndices.len <- bedEntry.list[[m]][["rightGeneIndicesLen"]]
      rightJunctionExon    <- bedEntry.list[[m]][["rightJunctionExon"]]
      rightJunctionSuffix  <- bedEntry.list[[m]][["rightJunctionSuffix"]]
      rightBedEntries      <- bedEntry.list[[m]][["rightBedEntries"]]
      
      for (i in 1:leftGeneIndices.len) {
        
        chromosomeStrand.bool <-
          rightStrand.list[[m]] == leftStrand.list[[m]][i] & rightChromosome.list[[m]] == leftChromosome.list[[m]][i]        
        if (leftStrand.list[[m]][i] == "+") {
          compatibleRightJunction.bool <- chromosomeStrand.bool & rightJunctionExon > leftJunctionExon[i]
        } else {
          compatibleRightJunction.bool <- chromosomeStrand.bool & rightJunctionExon < leftJunctionExon[i]
        }

        for (j in which(compatibleRightJunction.bool)) {
          bedEntry.num <- dim(leftBedEntries[[i]])[1] + dim(rightBedEntries[[j]])[1]
          if (count > 100000) {
            print(paste("Gene ", gene.id[m], ": ", l, " of ", bedRecord.num, " BED entries processed.", sep=""))
            count <- 0
          }
          
          junctionId[k] <- paste(gene.id[m], "-junction-", leftJunctionExon[i], "-", rightJunctionExon[j], ":", leftJunctionPrefix[i],
                                 ":", rightJunctionSuffix[j], sep="")

          ## Note for exon.sep == 0, two entries could still abut (but not overlap since the left boundary BED coordinate is
          ## decreased by one)
          leftBedEntries.num  <- dim(leftBedEntries[[i]])[1]
          if (as.integer(leftBedEntries[[i]][leftBedEntries.num, 3]) + exon.sep <= as.integer(rightBedEntries[[j]][1, 2])) {
            leftBedEntries[[i]][, 4] <- rep(junctionId[k], leftBedEntries.num)
            leftBedEntries[[i]][,11] <- paste(leftBedEntries[[i]][,4], "-exon-", 1:leftBedEntries.num, sep="")
            bedCoordinates[l:(l + leftBedEntries.num - 1),] <- leftBedEntries[[i]]
            l <- l + leftBedEntries.num
            count <- count + leftBedEntries.num
            
            leftBedEntry.len <- sum(as.integer(leftBedEntries[[i]][,10]))
            
            rightBedEntries.num <- dim(rightBedEntries[[j]])[1]
            rightBedEntries[[j]][, 4] <- rep(junctionId[k], rightBedEntries.num)
            rightBedEntries[[j]][,11] <- paste(rightBedEntries[[j]][,4], "-exon-", 1:rightBedEntries.num + leftBedEntries.num, sep="")
            bedCoordinates[l:(l + rightBedEntries.num - 1),] <- rightBedEntries[[j]]

            ## Correct the start coordinates of the right half exons. Note that we cannot do this on
            ## rightBedEntries[[j]] as we use rightBedEntries[[j]] several times
            bedCoordinates[l:(l + rightBedEntries.num - 1), 8] <-
              as.integer(bedCoordinates[l:(l + rightBedEntries.num - 1), 8]) + leftBedEntry.len
            bedCoordinates[l:(l + rightBedEntries.num - 1), 9] <-
              as.integer(bedCoordinates[l:(l + rightBedEntries.num - 1), 9]) + leftBedEntry.len
            l <- l + rightBedEntries.num
            count <- count + rightBedEntries.num
          
            k <- k + 1
          } else {
            skipped.entries.num <- skipped.entries.num + 1
            skipped.entries.ind <- c(skipped.entries.ind, i, j)
          }

        }
      }
    }
  }

  print(paste(skipped.entries.num, "combined bed entries skipped."))

  return(bedCoordinates[1:(l-1),])
  
}



###############################################################################
##
## computeJunctionBedLengths
##
###############################################################################

computeJunctionBedLengths <- function (exons) {

  return(sum(geneExonLengths(exons)))

}
