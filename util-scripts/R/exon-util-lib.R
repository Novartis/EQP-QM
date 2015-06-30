#call with: source(file.path(Sys.getenv("HOME"), "projects", "exons", "R", "exon-util-lib.R"))


###############################################################################
##
## Remove trailing spaces
##
###############################################################################

## returns string w/o leading whitespace
trim.leading <- function (x)  sub("^\\s+", "", x)

## returns string w/o trailing whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)

## returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)


###############################################################################
##
## Capitalize words
##
###############################################################################

## and the better, more sophisticated version:
capitalize <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
          {s <- substring(s, 2); if(strict) tolower(s) else s},
          sep = "", collapse = " ")
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


###############################################################################
##
## Convert index vector to a logical vector
##
###############################################################################

to.bool <- function (x, len) {

  x.bool <- array(FALSE, dim=len)

  if (length(! is.na(x)) > 0) {
    x.bool[x[! is.na(x)]] <- TRUE
  }

  if (length(is.na(x)) > 0) {
    x.bool[x[is.na(x)]] <- NA
  }
  
  return (x.bool)
}


###############################################################################
##
## lengthFirst: get the length of the first element
##
###############################################################################

lengthFirst <- function (x) {
  return (length(x[[1]]))
}


###############################################################################
##
## sort.df: sort a data frame by columns
##
###############################################################################

sort.df <- function (df, sortVar, ...) {

  return (df[with(df, order(sortVar, ...)), ])
  
}


###############################################################################
##
## function to parse gff attribute fields
##
###############################################################################

getAttributeFieldGff <- function (x, field, attrsep = ";") {
  s = strsplit(gsub("\"", "", x), split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(gsub("^ *", "", atts), split = " ", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      field.value <- a[[m]][2]
    }
    else {
      field.value <- as.character(NA)
    }
    return(field.value)
  })
}

###############################################################################
##
## function to read gff files
##
###############################################################################

readGffFile <- function(gff.file, nrows = -1) {

  gtf.colClasses <- c("character", "character", "character", "integer", "integer", "character", "character", "character", "character")
  
  print(paste("Reading ", gff.file, ": ", sep=""))
  gff.frame <- read.table(gff.file, sep="\t", as.is=TRUE, quote="", header=FALSE, comment.char="#", nrows = nrows, colClasses=gtf.colClasses,
                          stringsAsFactors=FALSE)
  colnames(gff.frame) <- c("Chromosome", "Source", "Feature type", "Start", "End", "Score", "Strand", "Frame", "Attributes")
  
  print(paste("Found", nrow(gff.frame), "rows."))
  
  stopifnot(!any(is.na(gff.frame[["Start"]])), !any(is.na(gff.frame[["End"]])))
  return(gff.frame)
  
}


###############################################################################
##
## convert to matrix
##
###############################################################################

convertToMatrix <- function (x, sep="/") {

  x.split <- strsplit(x[1], sep)
  ncol    <- length(x.split[[1]])

  x.mat <- array("", dim=c(length(x), ncol))
  for (i in 1:length(x)) {
    x.mat[i,] <- strsplit(x[i], sep)[[1]]
  }

  return(x.mat)
  
}

###############################################################################
##
## Get organism short
##
###############################################################################

getOrganismShort <- function (organism, sys.variable="ORGANISM_SHORT", organism.short.map=NA) {

  if (length(organism.short.map) == 1 && is.na(organism.short.map)) {
    organism.short.map <- c("hs", "mm", "rn", "susScr", "danRer", "plasCyn", "Mmul", "Mfas", "Mfas")
    names(organism.short.map) <- c("human", "mouse", "rat", "pig", "zebrafish", "plasmodium", "rhesus", "cynomolgous", "cyno")
  }
  
  organism.short <- Sys.getenv(sys.variable, ifelse (organism %in% names(organism.short.map), organism.short.map[[organism]], "unknown"))

  return(organism.short)
  
}


###############################################################################
##
## Get genome 
##
###############################################################################

getGenome <- function (organism, sys.variable="ORG_ID", genome.map=NA) {

  if (length(genome.map) == 1 && is.na(genome.map)) {
    genome.map <- c("hg19", "mm9", "rn4", "susScr3", "danRer7", "rheMac3", "fasMac4", "fasMac4")
    names(genome.map) <- c("human", "mouse", "rat", "pig", "zebrafish", "rhesus", "cynomolgous", "cyno")
  }
  
  genome <- Sys.getenv(sys.variable, ifelse (organism %in% names(genome.map), genome.map[[organism]], "unknown"))

  return(genome)
  
}



###############################################################################
##
## Get genome build
##
###############################################################################

getGenomeBuild <- function (organism, sys.variable="GENOME_BUILD", genome.build.map=NA) {

  if (length(genome.build.map) == 1 && is.na(genome.build.map)) {
    genome.build.map <- c("GRCh37", "GRCm38")
    names(genome.build.map) <- c("human", "mouse")
  }
  
  genome.build <- Sys.getenv(sys.variable, ifelse (organism %in% names(genome.build.map), genome.build.map[[organism]], "unknown"))

  return(genome.build)
  
}


#######################################################################################
##
##   Read a simple fasta file
##
#######################################################################################

readFastaFile <- function (fileName) {

  lines <- readLines(fileName)

  identifiers <- c()
  sequences <- c()
  sequence <- ""
  k <- 1
  not_found.num <- 0
  not_found.ids <- c()
  for (i in 1:length(lines)) {
    firstChar <- substr(lines[i], 1, 1)

    if (firstChar == ">") {
      identifier <- strsplit(lines[i], " ")[[1]][1]
      identifier <- substr(identifier, 2, nchar(identifier))

      identifier.ind <- match(identifier, transcript.ids)
      if (is.na(identifier.ind)) {
        not_found.num <- not_found.num + 1
        not_found.ids <- c(not_found.ids, identifier)
      }
      
      identifiers <- c(identifiers, identifier)
      if (k > 1) {
        sequences <- c(sequences, paste(lines[i_start:(i-1)], collapse=""))
      }
      i_start <- i+1
      k <- k + 1

      if (k %% 1000 == 0) {
        print(k)
      }
    } 
  }

  
  sequences <- c(sequences, paste(lines[i_start:length(lines)], collapse=""))

  if (not_found.num > 0) {
    print(paste(not_found.num, "transcript ids not found."))
  }

  return (list(identifiers, sequences))
}


###############################################################################
##
##  checkExonCoverage
##
###############################################################################

checkExonCoverage <- function (exons, separator="/") {

  exons <- unique(exons)

  if (length(exons) == 1) {
    return (TRUE)
  }
    
  exons.mat <- convertToMatrix(exons, separator)

  startIndex <- 1
  exons.chr    <- exons.mat[,startIndex]
  exons.start  <- as.integer(exons.mat[,startIndex + 1])
  exons.end    <- as.integer(exons.mat[,startIndex + 2])
  exons.strand <- exons.mat[,startIndex + 3]
  
  exons.order <- order (exons.chr, exons.start, exons.end)

  if (sum(exons.start[exons.order][2:length(exons)] == exons.end[exons.order][1:(length(exons)-1)] + 1) != length(exons)-1) {
    return (FALSE)
  }

  return (TRUE)

}


###############################################################################
##
##  createCassetteExons
##
###############################################################################

createCassetteExons <- function (exons, separator="/", printCounter=FALSE) {

  if (0 == 1) {
    test.chr <- c(rep("X", 4), rep("Y",2))
    
    test.start <- c(5,  8, 11, 14, 2, 5)
    test.end   <- c(7, 12, 13, 16, 6, 8)
    test.num   <- length(test.chr)

    test.bool <- test.chr[1:(test.num-1)] == test.chr[2:test.num] & test.start[2:test.num] <= test.end[1:(test.num-1)]

    test.chr[c(T, ! test.bool)]
    test.start[c(T, ! test.bool)]
    test.end[1:(test.num-1)][test.bool] <- test.start[2:test.num][test.bool]

  }

  exons.num <- length(exons)
  exons.mat <- convertToMatrix(exons, separator)

  startIndex <- 1
  exons.chr    <- exons.mat[,startIndex]
  exons.start  <- as.integer(exons.mat[,startIndex + 1])
  exons.end    <- as.integer(exons.mat[,startIndex + 2])
  exons.strand <- exons.mat[,startIndex + 3]
  
  exons.order <- order (exons.chr, exons.start, exons.end)

  exons.start.end.order <- order (c(exons.chr, exons.chr), c(exons.start, exons.end), c(rep("a", length(exons.start)), rep("b", length(exons.end))))

  exons.start.end.num  <- length (exons.start.end.order)
  exons.chr.sort       <- c(exons.chr, exons.chr)[exons.start.end.order]
  exons.start.end.sort <- c(exons.start, exons.end)[exons.start.end.order]
  exons.label.sort     <- c(rep("s", length(exons.start)), rep("e", length(exons.end)))[exons.start.end.order]
  exons.strand.sort    <- c(exons.strand, exons.strand)[exons.start.end.order]
  exons.ids.sort       <- c(exons, exons)[exons.start.end.order]

  exons.chr.cassette    <- array("", dim=exons.start.end.num)
  exons.start.cassette  <- array(0, dim=exons.start.end.num)
  exons.end.cassette    <- array(0, dim=exons.start.end.num)
  exons.strand.cassette <- array("", dim=exons.start.end.num)
  exons.cassette.list   <- vector ("list", exons.start.end.num)
  
  exon.count <- 0
  exon.list  <- c()
  start      <- 1
  for (i in 1:exons.start.end.num) {

    if (i %% 100000 == 0 && printCounter) {
      print (paste(as.integer(i / 2), "exons processed"))
    }

    if (i > 1 && exons.chr.sort[i-1] != exons.chr.sort[i]) {
      stopifnot(exon.count == 0)
    }
    
    if ((start < exons.start.end.sort[i] || (start == exons.start.end.sort[i] && exons.label.sort[i] == "e")) &&
        exon.count > 0) {
      exons.chr.cassette[i]    <- exons.chr.sort[i]
      exons.start.cassette[i]  <- start
      exons.cassette.list[[i]] <- exon.list
      if (exons.label.sort[i] == "s") {
        exons.end.cassette[i] <- exons.start.end.sort[i] - 1
        start <- exons.start.end.sort[i]
        exons.strand.cassette[i] <- exons.strand.sort[i-1]
      } else {
        exons.end.cassette[i] <- exons.start.end.sort[i]
        start <- exons.start.end.sort[i] + 1
        exons.strand.cassette[i] <- exons.strand.sort[i-1]
      }
    } else if (exon.count == 0) {
      start <- exons.start.end.sort[i]
      stopifnot (exons.label.sort[i] == "s")
    }
    
    if (exons.label.sort[i] == "s") {
      exon.count <- exon.count + 1
      exon.list  <- union(exon.list, exons.ids.sort[i])
    } else {
      exon.count <- exon.count - 1
      exon.list  <- setdiff (exon.list, exons.ids.sort[i])
    }

    stopifnot (exon.count == length (exon.list))

  }

  exons.cassette.list.len <- sapply(exons.cassette.list, length)
  cassette.exon.ind       <- which(exons.cassette.list.len > 0)

  cassette.exons <- paste(exons.chr.cassette[cassette.exon.ind], as.integer(exons.start.cassette[cassette.exon.ind]),
                          as.integer(exons.end.cassette[cassette.exon.ind]), exons.strand.cassette[cassette.exon.ind], sep="/")

  cassette.exons.rep <- rep(cassette.exons, exons.cassette.list.len[cassette.exon.ind])
  
  exons.cassette.orig <- unlist(exons.cassette.list)
  
  cassette.exons.rep.mat <- convertToMatrix (cassette.exons.rep, sep="/")
  
  startIndex <- 1
  cassette.exons.rep.chr    <- cassette.exons.rep.mat[,startIndex]
  cassette.exons.rep.start  <- as.integer(cassette.exons.rep.mat[,startIndex + 1])
  cassette.exons.rep.end    <- as.integer(cassette.exons.rep.mat[,startIndex + 2])
  cassette.exons.rep.strand <- cassette.exons.rep.mat[,startIndex + 3]

  exons.cassette.orig.mat <- convertToMatrix (exons.cassette.orig, sep="/")
  
  startIndex <- 1
  exons.cassette.orig.chr    <- exons.cassette.orig.mat[,startIndex]
  exons.cassette.orig.start  <- as.integer(exons.cassette.orig.mat[,startIndex + 1])
  exons.cassette.orig.end    <- as.integer(exons.cassette.orig.mat[,startIndex + 2])
  exons.cassette.orig.strand <- exons.cassette.orig.mat[,startIndex + 3]

  if (sum(cassette.exons.rep.chr != exons.cassette.orig.chr) != 0) {
    print(paste(sum(cassette.exons.rep.chr != exons.cassette.orig.chr), "cassette exons have different chromosomes than the original exons."))
  }

  if (sum(cassette.exons.rep.start < exons.cassette.orig.start) != 0) {
    print(paste(sum(cassette.exons.rep.start < exons.cassette.orig.start), "cassette exons have a start that is smaller than the start in the original exons."))
  }

  if (sum(cassette.exons.rep.end   > exons.cassette.orig.end) != 0) {
    print(paste(sum(cassette.exons.rep.end   > exons.cassette.orig.end), "cassette exons have an end that is larger than the end in the original exons."))
  }

  cassette.exons.diff.bool <- gsub("/[+-]$", "", cassette.exons.rep) != gsub("/[+-]$", "", exons.cassette.orig)
  exons.cassette.orig.bool <- tapply (cassette.exons.rep[cassette.exons.diff.bool], exons.cassette.orig[cassette.exons.diff.bool], checkExonCoverage)

  if (sum (! exons.cassette.orig.bool) != 0) {
    print(paste(sum (! exons.cassette.orig.bool), "original exons are not covered by their cassette exons."))
    exons.cassette.orig.uncovered <- levels(as.factor(exons.cassette.orig[cassette.exons.diff.bool]))[! exons.cassette.orig.bool]
    exons.cassette.uncovered <-
      cassette.exons.rep[cassette.exons.diff.bool][grep (paste(exons.cassette.orig.uncovered, collapse="|"), exons.cassette.orig[cassette.exons.diff.bool])]
    print(paste("Original exons:", paste(exons.cassette.orig.uncovered, collapse=", ")))
    print(paste("Cassette exons:", paste(exons.cassette.uncovered, collapse=", ")))
  }


  exons.cassette.orig.min.start <-
    tapply (exons.cassette.orig.start[cassette.exons.diff.bool], exons.cassette.orig[cassette.exons.diff.bool], min)
  exons.cassette.rep.min.start  <-
    tapply (cassette.exons.rep.start[cassette.exons.diff.bool], exons.cassette.orig[cassette.exons.diff.bool], min)
  exons.cassette.orig.levels <- levels(as.factor(exons.cassette.orig[cassette.exons.diff.bool]))

  if (sum(exons.cassette.orig.min.start != exons.cassette.rep.min.start) != 0) {
    print(paste(sum(exons.cassette.orig.min.start != exons.cassette.rep.min.start),
                "original exons have a different start than the cassette exons they are covered by."))
    exons.orig.diff.start <- unique(exons.cassette.orig.levels[which(exons.cassette.orig.min.start != exons.cassette.rep.min.start)])
    for (i in 1:length(exons.orig.diff.start)) {
      exons.orig.diff.start.ind <- grep (exons.orig.diff.start[i], exons.cassette.orig[cassette.exons.diff.bool])
      if (length(exons.orig.diff.start.ind) == 1) {
        print(paste("Correcting", exons.orig.diff.start[i]))
        exons.cassette.orig.start[cassette.exons.diff.bool][exons.orig.diff.start.ind] <- exons.orig.diff.start[i]
      } else {
        print(paste("Exon with differing start position covered by more than one cassette exon:",  exons.cassette.orig.start))
      }
    }   
  }

  exons.cassette.orig.max.end <-
    tapply (exons.cassette.orig.end[cassette.exons.diff.bool], exons.cassette.orig[cassette.exons.diff.bool], max)
  exons.cassette.rep.max.end  <-
    tapply (cassette.exons.rep.end[cassette.exons.diff.bool], exons.cassette.orig[cassette.exons.diff.bool], max)

  if (sum(exons.cassette.orig.max.end != exons.cassette.rep.max.end) != 0) {
    print(paste(sum(exons.cassette.orig.max.end != exons.cassette.rep.max.end),
                "original exons have a different end than the cassette exons they are covered by."))
  }

  return (list(exons=cassette.exons.rep, exons.orig=exons.cassette.orig))
  
}


###############################################################################
##
##  computeChromosomeInterExonRegions
##
###############################################################################

computeInterExonRegionsChromosom <- function (exons, separator="/") {

  exons.num <- length(exons)
  exons.mat <- convertToMatrix(exons, separator)

  exons.chr    <- exons.mat[,1]
  exons.start  <- as.integer(exons.mat[,2])
  exons.end    <- as.integer(exons.mat[,3])
  exons.strand <- exons.mat[,4]
  
  stopifnot(length(unique(exons.chr)) == 1)
  
  exons.order <- order(exons.start, exons.end)

  interExon.chr   <- c(exons.chr[1], exons.chr[exons.order])

  interExon.start <- c(1, exons.end[exons.order] + 1)

  stopifnot (exons.start[exons.order][1] != 1)
  
  interExon.end   <- c(exons.start[exons.order] - 1, as.integer(10^9 + 1))
  interExon.strand <- rep("+", length(exons.chr) + 1)

  return (paste(interExon.chr, as.integer(interExon.start), as.integer(interExon.end), interExon.strand, sep="/"))
}


###############################################################################
##
##  computeInterExonRegions
##
###############################################################################

computeInterExonRegions <- function (exons, separator="/", printCounter=FALSE) {

  if (0 == 1) {
    test.chr <- c(rep("X", 4), rep("Y",2))
    
    test.start <- c(5,  8, 11, 14, 2, 5)
    test.end   <- c(7, 12, 13, 16, 6, 8)
    test.num   <- length(test.chr)

    test.bool <- test.chr[1:(test.num-1)] == test.chr[2:test.num] & test.start[2:test.num] <= test.end[1:(test.num-1)]

    test.chr[c(T, ! test.bool)]
    test.start[c(T, ! test.bool)]
    test.end[1:(test.num-1)][test.bool] <- test.start[2:test.num][test.bool]

  }

  exons.num <- length(exons)
  exons.mat <- convertToMatrix(exons, separator)

  startIndex <- 1
  exons.chr    <- exons.mat[,startIndex]
  exons.start  <- as.integer(exons.mat[,startIndex + 1])
  exons.end    <- as.integer(exons.mat[,startIndex + 2])
  exons.strand <- exons.mat[,startIndex + 3]
  
  exons.order <- order (exons.chr, exons.start, exons.end)

  exons.start.end.order <- order (c(exons.chr, exons.chr), c(exons.start, exons.end))

  exons.start.end.num  <- length (exons.start.end.order)
  exons.chr.sort       <- c(exons.chr, exons.chr)[exons.start.end.order]
  exons.start.end.sort <- c(exons.start, exons.end)[exons.start.end.order]
  exons.label.sort     <- c(rep("s", length(exons.start)), rep("e", length(exons.end)))[exons.start.end.order]
  exons.strand.sort    <- c(exons.strand, exons.strand)[exons.start.end.order]
  exons.ids.sort       <- c(exons, exons)[exons.start.end.order]

  exons.chr.start.end.non_ov <- array("", dim=exons.start.end.num)
  exons.start.end.non_ov     <- array(0, dim=exons.start.end.num)

  exon.count <- 0
  start      <- 1
  for (i in 1:exons.start.end.num) {

    if (i %% 100000 == 0 && printCounter) {
      print (paste(as.integer(i / 2), "exons processed."))
    }

    if (i > 1 && exons.chr.sort[i-1] != exons.chr.sort[i]) {
      stopifnot(exon.count == 0)
    }
    
    if (exons.label.sort[i] == "s") {
      if (exon.count == 0) {
        exons.chr.start.end.non_ov[i] <- exons.chr.sort[i]
        exons.start.end.non_ov[i] <- exons.start.end.sort[i]
      }
      exon.count <- exon.count + 1
    } else {
      exon.count <- exon.count - 1
      if (exon.count == 0) {
        exons.chr.start.end.non_ov[i] <- exons.chr.sort[i]
        exons.start.end.non_ov[i] <- exons.start.end.sort[i]
      }
    }

  }

  exons.chr.non_ov.empty.bool <- exons.chr.start.end.non_ov != ""

  exons.chr.non_ov.mat <- matrix (exons.chr.start.end.non_ov[exons.chr.non_ov.empty.bool], ncol=2, byrow=TRUE)
  stopifnot(sum(exons.chr.non_ov.mat[,1] != exons.chr.non_ov.mat[,2]) == 0)
  
  exons.start.end.non_ov.mat <- matrix (exons.start.end.non_ov  [exons.chr.non_ov.empty.bool], ncol=2, byrow=TRUE)
  stopifnot(sum(exons.start.end.non_ov.mat[,1] > exons.start.end.non_ov.mat[,2]) == 0)

  exons.chr.non_ov   <- exons.chr.non_ov.mat[,1]
  exons.start.non_ov <- exons.start.end.non_ov.mat[,1]
  exons.end.non_ov   <- exons.start.end.non_ov.mat[,2]
  exon.non_ov.num    <- length(exons.chr.non_ov)
  
  exon.overlap.check.bool <- (exons.chr.non_ov[1:(exon.non_ov.num-1)] == exons.chr.non_ov[2:exon.non_ov.num] &
                              exons.end.non_ov[1:(exon.non_ov.num-1)] >= exons.start.non_ov[2:exon.non_ov.num])

  stopifnot(sum(exon.overlap.check.bool) == 0)

  exons.non_ov     <- paste(exons.chr.non_ov, exons.start.non_ov, exons.end.non_ov, rep("+", exon.non_ov.num ), sep=separator)
  interExonRegions <- unlist(tapply(exons.non_ov, exons.chr.non_ov, computeInterExonRegionsChromosom, separator))

  return (interExonRegions)
  
}

