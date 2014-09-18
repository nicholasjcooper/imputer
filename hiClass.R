if(getwd()!= "/home/ncooper"){
  require(genoset)
}

require(NCmisc)

## depends on iFunctions ##

#' Class 'HiC'
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' SLOTS
#'  seqnames, seqnames.end ranges, ranges.end, [strand], seqinfo, seqinfo.end,
#'  build, and elementMetadata, containing columns:
#'  score, reads
#'  LISTS: promoters.bait, promoters.end
#' METHODS
#'  "[[", "[", show, print, length, initialize
#'  build, promo, score, reads, bait, otherEnd, chr2, start2, end2, genes, 
#'  chr, start, end, head, tail, show, print, readconvTo36, convTo37
#' COERCION
#'  can use 'as' to convert to : GRanges, RangedData, data.frame
setClass("HiC",
         contains="GRanges",
         representation(
           seqnames="Rle",
           seqnames.end="Rle",
           ranges="IRanges",
           ranges.end="IRanges",
           strand="Rle",
           elementMetadata="DataFrame",
           seqinfo="Seqinfo",
           seqinfo.end="Seqinfo",
           build="character",
           promoters.bait="list",
           promoters.end="list"
         ),
         prototype(
           seqnames=Rle(factor()),
           seqnames.end=Rle(factor()),
           strand=Rle(strand()),
           seqinfo=Seqinfo(),
           seqinfo.end=Seqinfo(),
           ranges=IRanges(),
           ranges.end=IRanges(),
           build=character(), 
           elementMetadata=DataFrame(reads=NULL, score=NULL),
           promoters.bait=list(),
           promoters.end=list()
         )
)

#' Length for HiC
#' 
#' Simply returns the number of SNPs (internally using 'nrow')
#' @param x a HiC object
#' @return integer
setMethod("length", "HiC", function(x) length(x@ranges))

setGeneric("bait", function(x) standardGeneric("bait"))
setGeneric("otherEnd", function(x) standardGeneric("otherEnd") )
setGeneric("ranges.bait", function(x) standardGeneric("ranges.bait"))
setGeneric("ranges.end", function(x) standardGeneric("ranges.end") )


#' Retrieve the Chip name for HiC
#' 
#' Simply returns the name of the chip, e.g, 'ImmunoChip'
#' @param x a HiC object
#' @return character
setMethod("bait", "HiC", function(x) emd.rmv(as(x,"GRanges")) )
setMethod("otherEnd", "HiC", function(x) { x@ranges <- x@ranges.end; return(emd.rmv(as(x,"GRanges"))) } )
if(!exists("build",mode="function")) {
 setGeneric("build", function(x) standardGeneric("build") )
}
# alternative names for these two GRanges access functions
setMethod("ranges.bait", "HiC", function(x) { bait(x) } )
setMethod("ranges.end", "HiC", function(x) { otherEnd(x) } )


#' Retrieve the build for HiC
#' 
#' Returns the UCSC build of the chip object, e.g, "hg18" or "hg19"
#' @param x a HiC object
#' @return character, "hg18" or "hg19"
setMethod("build", "HiC", function(x) x@build)
setGeneric("reads", function(x,b=TRUE) standardGeneric("reads") )

#' Access rs-ids for HiC
#' 
#' Returns the reads for the HiC object, e.g, 27, etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a HiC object
#' @return integer vector of read counts (or NAs)
setMethod("reads", "HiC", function(x,b=TRUE) { 
  u <- mcols(x) ;  
  if("reads" %in% colnames(u)) { 
    U <- u[,"reads"] 
    if(!b) { U <- gsub("b","",U) }
    return(U)
  } else { return(NULL) } 
})

if(!exists("score",mode="function")) {
  setGeneric("score", function(x) standardGeneric("score") )
}

#' Access score for HiC
#' 
#' Returns score column for the object
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a HiC object
#' @return numeric vector of scores (or NAs)
setMethod("score", "HiC", function(x) { u <- mcols(x) ;  if("score" %in% colnames(u)) { u[,"score"] } else { NULL } })

setGeneric("promo", function(x,bait=TRUE) standardGeneric("promo") )
setGeneric("pro.bait", function(x) standardGeneric("pro.bait") )
setGeneric("pro.end", function(x) standardGeneric("pro.end") )


#' Access promoter lists for HiC
#' 
#' Returns the letter for the second alleles for the chip object, 
#' e.g, "A","C","G","T", etc
#' Only if these are annotated internally, or else a vector of NAs
#' @param x a HiC object
#' @return character vector of allele codes (or NAs)
setMethod("promo", "HiC", function(x,bait=TRUE) { if(bait) { x@promoters.bait } else { x@promoters.end } })
# alternatives to ^
setMethod("pro.bait", "HiC", function(x) { x@promoters.bait })
setMethod("pro.end", "HiC", function(x) { x@promoters.end })



#' Subset HiC by chromosome
#' 
#' Returns the subset of the HiC object for which SNPs are on
#' the chromosome specified, by either number or character.
#' @param x a HiC object
#' @param i a chromosome number or letter, i.e, one of seqlevels(x)
#' @return HiC object for the subset of SNPs on chromosome i
setMethod("[[", "HiC", function(x,i,j,...) { 
  dotArgs <- list(...)
  if (length(dotArgs) > 0)
    dotArgs <- dotArgs[names(dotArgs) != "exact"]
  if (!missing(j) || length(dotArgs) > 0)
    stop("invalid subsetting")
  if (missing(i))
    stop("subscript is missing")
  if (!is.character(i) && !is.numeric(i)) 
    stop("invalid subscript type")
  if (length(i) < 1L)
    stop("attempt to select less than one element")
  if (length(i) > 1L)
    stop("attempt to select more than one element")
  cn <- chrNames(x)
  if (is.numeric(i) && !is.na(i) && (i < 1L || i > length(cn)))
    stop("subscript out of bounds")
  # do the selection #
  if(i %in% paste(chr(x))) {
    out <- x[chr(x)==i,]
  } else {
    if(is.numeric(i)) {
      out <- x[match(chr(x),chrNames(x))==i,]
    } else {
      stop("unknown index")
    }
  }
  out@build <- x@build
  return(out)
} )


setMethod("head", "HiC", function(x,...) { 
    .local <- function (x, n = 6L, ...) 
    {
      if (!isSingleNumber(n)) 
        stop("'n' must be a single integer")
      if (!is.integer(n)) 
        n <- as.integer(n)
      x_NROW <- NROW(x)
      if (n >= 0L) {
        n <- min(n, x_NROW)
      }
      else {
        n <- max(x_NROW + n, 0L)
      }
      x[1:n,]
    }
    .local(x, ...)
} )

setMethod("tail", "HiC", function(x,...) { 
  .local <- function (x, n = 6L, ...) 
  {
    if (!isSingleNumber(n)) 
      stop("'n' must be a single integer")
    if (!is.integer(n)) 
      n <- as.integer(n)
    x_NROW <- NROW(x)
    if (n >= 0L) {
      n <- min(n, x_NROW)
    }
    else {
      n <- max(x_NROW + n, 0L)
    }
    x[(x_NROW-n+1):x_NROW,]
  }
  .local(x, ...)
} )


setGeneric("start2", function(x) standardGeneric("start2"))
setGeneric("end2", function(x) standardGeneric("end2"))
setGeneric("otherChr", function(x) standardGeneric("otherChr"))
setGeneric("chr.end", function(x) standardGeneric("chr.end"))


setMethod("start2", "HiC", function(x) { 
  start(otherEnd(x))
} )

setMethod("end2", "HiC", function(x) { 
  end(otherEnd(x))
} )

setMethod("otherChr", "HiC", function(x) { chr(otherEnd(x)) } )
setMethod("chr.end", "HiC", function(x) { otherChr(x) } ) # duplicate function of otherChr


if(!exists("genes",mode="function")) {
 setGeneric("genes", function(x,bait=TRUE,end=TRUE,list=TRUE,missing=".") standardGeneric("genes"))
}

setGeneric("genes.bait", function(x,list=TRUE,missing=".") standardGeneric("genes.bait"))
setGeneric("genes.end", function(x,list=TRUE,missing=".") standardGeneric("genes.end"))

setMethod("genes", "HiC", function(x,bait=TRUE,end=TRUE,list=TRUE,missing=".") { 
  if(bait){ g1 <- pro.bait(x) }
  if(end) { g2 <- pro.end(x) }
  if(bait & end) {
    gg <- mapply(c,g1,g2)
    out <- extract.genes(gg,list=list,missing=missing)
  } else {
    if(bait) {
      #prv(g1)
      out <- extract.genes(g1,list=list,missing=missing)
    } else {
      #prv(g2)
      out <- extract.genes(g2,list=list,missing=missing)
    }
  }
  return(out)  
} )

setMethod("genes.bait", "HiC", function(x,list=TRUE,missing=".") { 
  genes(x,bait=TRUE,end=FALSE,list=list,missing=missing) } )
setMethod("genes.end", "HiC", function(x,list=TRUE,missing=".") { 
  genes(x,bait=FALSE,end=TRUE,list=list,missing=missing) } )

# internal function
stracta <- function(x) {
  chop <- strsplit(x,split = "-",fixed=TRUE)
  .local2 <- function(x) { lx <- length(x); if(lx>1) { paste(x[1:(lx-1)],collapse="") } else { x }}
  minus.suf <- unlist(lapply(chop,.local2))
  return(unique(minus.suf))
}

# internal function - requires 'comma'
extract.genes <- function(X,list=TRUE,missing=".") {
  if(!is.list(X)) { stop("this function is meant for a list object") }
  l.out <- lapply(X,stracta)
  l.out <- lapply(l.out,function(x) { y <- x[x!="."] ; if(length(y)>0) { return(y) } else { return(x[1]) } })
  l.out <- lapply(l.out,function(x) { x[x %in% c("",".",NA,"NA")] <- missing; return(x) })
  if(list) { 
    return(l.out) 
  } else {
    return(sapply(l.out,comma))
  }
}


#' Subset HiC by row (ranges)
#' 
#' Returns a row subset of the HiC object.
#' @param x a HiC object
#' @param i a chromosome number or letter, i.e, one of seqlevels(x)
#' @return HiC object for the subset of SNPs on chromosome i
setMethod("[", "HiC", function(x,i,j,...) { 
  ## all valid 'colname' values allowed for 'j' parameter
  cmd.list <- c("chr","chr.bait","seqnames","chr1","start","start.bait","start1","end","end.bait","end1",
                "chromosome","chromosome.bait","chromosome1","chromosome2","chromosome.end",
                "chr2","chr.end","seqnames.end","seqnames2","start.end","start2","end.end","end2",
                "promo","promo1","promo.bait","promoter","promoters","promoter.bait","promoters.bait",
                "promoter1","promoters1","pro","pro1","pro.bait",
                "promo2","promo.end","promoter.end","promoters.end","promoter2","promoters2","pro2","pro.end",
                "reads","n_reads","nreads","n.reads","score","scores","score.bait","score.end","scores.bait","scores.end")
  dotArgs <- list(...)
  if (length(dotArgs) > 0)
    dotArgs <- dotArgs[names(dotArgs) != "exact"] ##???
  if (length(dotArgs) > 0)
    stop("invalid subsetting")
  if(missing(i) & missing(j)) {
    return(x)
  }
  if (!missing(i)) {
    if(is.logical(i)) {
      i <- which(i)
    }
    if (!is.character(i) && !is.numeric(i)) 
      stop("invalid subscript type")
    if (length(i) < 1L)
      stop("attempt to select less than one element")
    if (is.numeric(i) && !is.na(i) && (i < 1L || i > length(x)))
      stop("subscript out of bounds")
    # do the selection #
    X <- IRanges:::extractROWS(x,i)
    n1 <- names(X@ranges)
    n2 <- names(x@ranges.end) # this seems a bit hacky but seems to work!
    ii <- match(n1,n2)
    if(any(is.na(ii))) { stop("mismatch between names in @ranges versus @ranges.end") }
    #cat("reduced length from ",length(n2),"to",length(ii),"\n")
    X@ranges.end <- x@ranges.end[ii,]
    X@promoters.bait <- x@promoters.bait[ii]
    X@promoters.end <- x@promoters.end[ii]
  } else { 
    X <- x # if no i index is used
  }
  if (!missing(j)) {
    lj <- length(j)
    out.list <- vector("list",lj)
    if(all(j %in% 1:10)) {
      for(cc in 1:lj) {
        if(j[cc]==1)  { out.list[[cc]] <- chr(X) }
        if(j[cc]==2)  { out.list[[cc]] <- start(X) }
        if(j[cc]==3)  { out.list[[cc]] <- end(X) }
        if(j[cc]==4)  { out.list[[cc]] <- otherChr(X) }
        if(j[cc]==5)  { out.list[[cc]] <- start2(X) }
        if(j[cc]==6)  { out.list[[cc]] <- end2(X) }
        if(j[cc]==7)  { tmp <- promo(X, TRUE) ; out.list[[cc]] <- if(lj>1) { sanitize.promo(tmp,TRUE) } else { tmp } }
        if(j[cc]==8)  { tmp <- promo(X, FALSE) ; out.list[[cc]] <- if(lj>1) { sanitize.promo(tmp,FALSE) } else { tmp } }
        if(j[cc]==9)  { out.list[[cc]] <- reads(X) }
        if(j[cc]==10) { out.list[[cc]] <- score(X) }
      }
    } else {
      if(all(tolower(j) %in% cmd.list)) {
        for(cc in 1:lj) {
          if(tolower(j[cc]) %in% c("chr","chr.bait","seqnames","chr1","chromosome","chromosome.bait","chromosome1"))  { 
            out.list[[cc]] <- chr(X) }
          if(tolower(j[cc]) %in% c("start","start.bait","start1"))  { out.list[[cc]] <- start(X) }
          if(tolower(j[cc]) %in% c("end","end.bait","end1"))  { out.list[[cc]] <- end(X) }
          if(tolower(j[cc]) %in% c("chr2","chr.end","seqnames.end","seqnames2","chromosome2","chromosome.end"))  {
            out.list[[cc]] <- otherChr(X) }
          if(tolower(j[cc]) %in% c("start.end","start2"))  { out.list[[cc]] <- start2(X) }
          if(tolower(j[cc]) %in% c("end.end","end2"))  { out.list[[cc]] <- end2(X) }
          if(tolower(j[cc]) %in% c("promo","promo1","promo.bait","promoter","promoters","pro","pro1","pro.bait",
                                   "promoter.bait","promoters.bait","promoter1","promoters1"))  { 
            out.list[[cc]] <- if(lj>1) { sanitize.promo(X,TRUE) } else { promo(X, TRUE) } }
          if(tolower(j[cc]) %in% c("promo2","promo.end","promoter.end","promoters.end","promoter2",
                                   "promoters2","pro2","pro.end"))  { 
            out.list[[cc]] <- if(lj>1) { sanitize.promo(X,FALSE) } else { promo(X, FALSE) } }
          if(tolower(j[cc]) %in% c("reads","n_reads","nreads","n.reads"))  { out.list[[cc]] <- reads(X) }
          if(tolower(j[cc]) %in% c("score","scores","score.bait","score.end","scores.bait","scores.end"))  { 
            out.list[[cc]] <- score(X) }
        }
      } else {
        if(any(j %in% 1:10) & any(tolower(j) %in% cmd.list)) { stop("cannot mix integer indexing with colname indexing") }
        stop("invalid column subscript, use integers 1:10, or standard column names: chr, chr.end, ranges, ranges.end,..., etc")
      }
    }
    if(lj==1) {
      X <- out.list[[1]] #leave as list
    } else {
      X <- do.call("cbind",args=out.list)
    }
  }
  return(X)
} )

setGeneric("convTo37", function(x) standardGeneric("convTo37"))
          
#' Convert HiC to build 37/hg19 coordinates
#' 
#' Returns the a HiC object with positions updated to build
#' 37 coordinates, assuming that the existing object was in build 36,
#' or already in build 37 coordinates, and that the build() slot was
#' entered correctly. Ensure that the value of build(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 37/hg19, but build(x)=="hg18" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
#' @param x a HiC object
#' @return HiC object with the build updated to hg19 coordinates
#' @seealso convTo36
setMethod("convTo37", "HiC", function(x) {
  if(!exists("conv.36.37",mode="function")) { warning("iFunctions.R not loaded, no conversion performed") ; return(x) }
  .local <- function(x) {
    if(ucsc.sanitizer(build(x))=="hg18") {
      u <- conv.36.37(ranges=as(x,"GRanges"))
      if(length(u)==length(x)) { 
        x@ranges <- u@ranges
        all.eq <- TRUE
        if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
        if(all.eq) { all.eq <- all(sort(seqlevels(x))==sort(seqlevels(u))) }
        if(!all.eq) {
          #print(seqlevels(x)); print(seqlevels(u))
          #warning("conversion altered the chromosomes"); 
          seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
        }
        xx <- as(x@seqnames,"character")
        uu <- as(u@seqnames,"character")
        if(any(xx!=uu)) { x@seqnames <- u@seqnames }
        x@build <- "hg19"
      } else { 
        stop("conversion to build37/hg19 failed, input had length ",length(x),"; output had length ",length(u)) 
      } 
    } else {
      if(ucsc.sanitizer(build(x))!="hg19") { 
        warning("input object was not tagged as hg18/build36 [@build], left unchanged") 
      } else {
        warning("object is already using hg19/build37, no change")
      }
    }
    return(x)
  }
  bb <- .local(bait(x))
  ee <- .local(otherEnd(x))
  x@build <- bb@build
  x@ranges <- bb@ranges
  x@ranges.end <- ee@ranges
  return(x)
})

setGeneric("convTo36", function(x) standardGeneric("convTo36"))


#' Convert HiC to build 36/hg18 coordinates
#' 
#' Returns the a HiC object with positions updated to build
#' 36 coordinates, assuming that the existing object was in build 37,
#' or already in build 36 coordinates, and that the build() slot was
#' entered correctly. Ensure that the value of build(x) is correct before
#' running this function for conversion; for instance, if the coordinates 
#' are already build 36/hg18, but build(x)=="hg19" (incorrect value), then
#' these coordinates will be transformed in a relative manner rendering the
#' result meaningless.
#' @param x a HiC object
#' @return HiC object with the build updated to hg18 coordinates
#' @seealso convTo37
setMethod("convTo36", "HiC", function(x) {
  if(!exists("conv.37.36",mode="function")) { warning("iFunctions.R not loaded, no conversion performed") ; return(x) }
  .local <- function(x) {
    if(ucsc.sanitizer(build(x))=="hg19") {
      u <- conv.37.36(ranges=as(x,"GRanges"))
      if(length(u)==length(x)) { 
        x@ranges <- u@ranges
        all.eq <- TRUE
        if(all.eq) { all.eq <- length(seqlevels(x))==length(seqlevels(u)) }
        if(all.eq) { all.eq <- any(sort(seqlevels(x))==sort(seqlevels(u))) }
        if(!all.eq) {
          #warning("conversion altered the chromosomes"); 
          seqlevels(x) <- c(seqlevels(x),seqlevels(u)[!seqlevels(u) %in% seqlevels(x)])  #x@seqinfo <- u@seqinfo 
        }
        xx <- as(x@seqnames,"character")
        uu <- as(u@seqnames,"character")
        if(any(xx!=uu)) { x@seqnames <- u@seqnames }
        x@build <- "hg18"
      } else { 
        stop("conversion to build36/hg18 failed") 
      } 
    } else {
      if(ucsc.sanitizer(build(x))!="hg18") { 
        warning("input object was not tagged as hg19/build37 [@build], left unchanged") 
      } else {
        warning("object is already using hg18/build36, no change")
      }
    }
    return(x)
  }
  bb <- .local(bait(x))
  ee <- .local(otherEnd(x))
  x@build <- bb@build
  x@ranges <- bb@ranges
  x@ranges.end <- ee@ranges
  return(x)
})



#' Display a HiC object
#' 
#' Returns a preview of a HiC object to the console. This
#' is similar to a GRanges preview, but the seqlevels are hidden, the UCSC
#' build and chip name are displayed, start and end are merged to the virtual
#' label 'position' (as it's assume we are dealing with SNPs, not ranges), the strand
#' by default is hidden, and the integer codes for pass/fail in QCcodes() are 
#' displayed as 'pass' or 'fail', even though this is not how they are represented internally. 
#' @param object a HiC object
#' @param up.to only SNPs at the start and end are generally displayed, however this
#' parameter specifies that when there are <= 'up.to' SNPs, then all SNPs will be displayed.
#' @param head.tail number of SNPs to display at start/end (only the head and tail are
#' shown as these objects are generally very large with >100K SNPs)
#' @param show.promo logical, by default the promoters are hidden, as these
#' take up a lot of space. Setting to TRUE will display them
#' @return displays a preview of the object on the console
setMethod("show", "HiC", 
     function(object) { showHiC(object,up.to=10,head.tail=5,show.promo=FALSE) } )

#' Print a HiC object to the console
#' 
#' See 'show' as the behaviour is very similar and ... are just arguments of 'show'.
#' The key difference with 'print' instead of 'show' is that by default the parameter
#' 'up.to' is set to 50, so that any HiC object (or subset) of less than or equal
#' to 50 rows will be displayed in its entirety, rather than just the top/bottom 5 rows. 
#' @param object a HiC object
#' @param ... further arguments to show()
#' @return displays a preview of the object on the console
setMethod("print", "HiC", 
          function(x,...) { showHiC(x,...) } )


#' Constructor for HiC annotation object
#' 
#' This class annotates a microarray SNP chip with data for each SNP including chromosome,
#' id, position, strand, 'rs' id, allele 1, allele 2 for each SNP of a microarray chip,
#' in either hg18 or hg19 (build 36/37) coordinates.
#' This package makes extension use of this class of annotation object for the working
#' microarray chip, e.g, default is ImmunoChip, but Metabochip is also built-in,
#' and you can also load your own annotation if using a different chip. The class
#' is basically a GRanges object, modified to always have columns for A1, A2 (alleles), 
#' rs.id, and a quality control flag. The default display is tidier than GRanges, it has
#' nice coersion to and frame data.frame and subsetting by chromosome using [[n]] has been
#' added, in addition to normal [i,j] indexing native to GRanges.
#' @param GRanges a GRanges object containing chromosome, start/end = position, and strand
#' information for the chip object to be created, also rownames should be used to code
#' the chip-ids for each SNP.
#' @param chr optional, alternative to using 'GRanges' to input SNP locations, enter here 
#' a vector of chromosome numbers/letters for each SNP. The recommended coding is: 
#' 1:22, X, Y, XY, MT
#' @param pos optional, vector of positions (integers), use in conjunction with 'chr' and
#'  'ids' as an alternative way to input SNP position information instead of GRanges.
#' @param ids optional, vector of SNP chip-ids, use in conjunction with 'chr' and
#'  'pos' as an alternative way to input SNP position information instead of GRanges.
#' @param chip character, name of the chip you are making this annotation for (only used
#' for labelling purposes)
#' @param build character, either "hg18" or "hg19". Will also accept build number, 36 or 37.
#' This indicates what coordinates the object is using, and will be taken into account by
#' conversion functions, and annotation lookup functions throughout this package.
#' @param rs.id 'rs' ids are standardized ids for SNPs, these usually differ from each chips'
#' own IDs for each snp. If you don't know these, or can't find them, they can be left blank,
#' but will render the functions 'rs.to.id()' and 'id.to.rs()' useless for this HiC object.
#' @param A1 the first allele letter code for each SNP, e.g, usually "A","C","G", or "T", but
#' you can use any scheme you like. Can be left blank.
#' @param A2, as for A1, but for allele 2.
#' @param QCcode optional column to keep track of SNPs passing and failing QC. You can completely
#' ignore this column. It works based on integer codes, 0,1,2, you may wish to use simple 0 and 1,
#' for pass and fail respectively, or else 0 can be pass, and 1,2,... can indicate failure for 
#' different criteria. 0 will always be treated as a pass and anything else as a fail, so you
#' can code fails however you wish.
HiC <- function(GRanges=NULL, GRanges.end=NULL, reads=NULL, score=NULL, build="",
                     promo.bait=NULL, promo.end=NULL,sep="[/||,]",fixed=FALSE) {
  if(build!="") { build <- ucsc.sanitizer(build) }
  LL <- length(GRanges)
  if(length(GRanges.end)!=LL) { stop("GRanges and GRanges.end must be the same size") }
  if(any(chr(GRanges)!=chr(GRanges.end))) { warning("GRanges and GRanges.end did not have identical chromosomes in each row") }
  if(length(reads)!=LL | length(score)!=LL) { reads <- score <- rep(NA,times=LL) }
  if(length(promo.bait)!=LL) { promo.bait <- rep(NA,times=LL) } 
  if(length(promo.end)!=LL) { promo.end <- rep(NA,times=LL) } 
  if(is.null(GRanges)) { stop("Granges was empty") }
  if(is(GRanges)[1]!="GRanges") { GRanges <- as(GRanges,"GRanges") }
  if(is.null(GRanges.end)) { stop("Granges.end was empty") }
  if(is(GRanges.end)[1]!="GRanges") { GRanges.end <- as(GRanges.end,"GRanges") }
  if(!is.list(promo.bait)) { promo1.list <- strsplit(paste(promo.bait), split = sep, fixed = fixed) } else { promo1.list <- promo.bait }
  if(!is.list(promo.end)) { promo2.list <- strsplit(paste(promo.end), split = sep, fixed = fixed) } else { promo2.list <- promo.end }
  df <- DataFrame(reads=reads,score=score)
  return(new("HiC", seqnames=GRanges@seqnames, ranges=GRanges@ranges, seqnames.end=GRanges.end@seqnames, ranges.end=GRanges.end@ranges, strand=GRanges@strand,
            elementMetadata=df, seqinfo=GRanges@seqinfo, seqinfo.end=GRanges.end@seqinfo,
            build=build, promoters.bait=promo1.list, promoters.end=promo2.list))
}


#' Initialize method for HiC
#' 
#' Please use the 'HiC' constructor
setMethod("initialize", "HiC",
              function(.Object, ...){
          		  callNextMethod(.Object, ...)
          	  })


setAs("HiC", "GRanges",
      function(from) { 
        #print(is(from)); print(from@seqnames)
        out <- GRanges(from@seqnames,ranges=from@ranges,strand=from@strand,
                       seqinfo=from@seqinfo,elementMetadata=from@elementMetadata) #,genome=build(from))
        return(out)
      }
)

setAs("HiC", "RangedData",
      function(from) { 
        out <- as(as(from,"GRanges"),"RangedData")
        if("strand" %in% colnames(out)) { out <- out[,-which(colnames(out) %in% "strand")] }
        return(out)
      }
)

setAs("HiC", "data.frame", function(from) { ranged.to.data.frame(as(from,"GRanges"),include.cols=TRUE) })


setValidity("HiC",
            function(object) {
              if (!is.character(build(object)) || length(build(object)) != 1 || is.na(build(object))) {
                return("'build' slot must be a single string") 
              } else {
                if(!build(object) %in% c("",ucsc.sanitizer(show.valid=T)[,1])) {
                  return("'build' must be a string, 36/37 or hg18/hg19") 
                }
              }
              if(length(object@promoters.bait)!=length(object@promoters.end)) { stop("mismatching lengths for promoter lists") }
              #print(length(object@promoters.bait)); print(length(object@ranges))
              #if(length(object@promoters.bait)!=length(object@ranges)) { stop("mismatching lengths for bait promoter list and ranges") }
              #if(length(object@promoters.end)!=length(object@ranges.end)) { stop("mismatching lengths for other-end promoter list and ranges") }
              if(nrow(mcols(object))!=length(object@ranges)) { stop("mismatching lengths between ranges and meta-data (ie, reads, score)") }
            }
)

#internal
showHiC <- function (x, margin = "", with.classinfo = FALSE, print.seqlengths = FALSE,...) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  score <- score(x)
  reads <- reads(x)
  bb <- build(x)
  if(bb=="") { bb <- "unknown" }
  #if(length(qc)==lx) { QC <- rep("pass",lx); QC[qc>0] <- paste0("fail",QC[qc>0]) ; x$QCcode <- QC }
  cat("HiC with ", lx, " ", ifelse(lx == 1L, "Entry", 
                   "entries")," using ",bb," coordinates",":\n", sep = "")
  out <- makePrettyMatrixForCompactPrinting2(x, .makeNakedMatFromHiC,...)
  if (nrow(out) != 0L) 
    rownames(out) <- paste0(margin, rownames(out))
  print(out, quote = FALSE, right = TRUE)
  if(ncol(out)==6) { cat("<promoter columns suppressed>\n") }
}


# convert the promoter lists into a 1 column summary
sanitize.promo <- function(x,bait=TRUE) {
  l1 <- promo(x,bait=bait)
  san.pro <- function(X) {
    ll <- length(X)
    if(ll>0) {
      if(ll<3) {
        out <- comma(X)
      } else {
        out <- paste(X[1],"+",ll-1,"more")
      }
    }
    return(out)
  }
  L1 <- paste(sapply(l1,san.pro))
  return(L1)
}

#internal
.makeNakedMatFromHiC <- function (x,show.promo=FALSE) 
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  x1 <- bait(x)
  x2 <- otherEnd(x)
  if(!show.promo) {
    ans <- cbind(seqnames = as.character(seqnames(x)), 
                 ranges = showAsCell(ranges(x1)), ranges.end = showAsCell(ranges(x2)))
  } else { 
    ans <- cbind(seqnames = as.character(seqnames(x)), 
                 ranges = showAsCell(ranges(x)), 
                 seqnames.end = as.character(seqnames(otherEnd(x))), 
                 ranges.end = showAsCell(ranges(otherEnd(x))),
                 pro.bait = sanitize.promo(x,T), 
                 pro.end = sanitize.promo(x,F) )
  }
  extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), lapply(GenomicRanges:::extraColumnSlots(x), 
                                              showAsCell)))
  }
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), showAsCell), 
                                 list(check.names = FALSE)))
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
    if(all(colnames(ans)[1:2]==c("seqnames","ranges"))) { colnames(ans)[1:2] <- c("chr","ranges.bait") }
    if(colnames(ans)[3]=="seqnames.end") { colnames(ans)[3] <- "chr.end" }
    ii <- which(colnames(ans) %in% c("pro.bait","pro.end"))
    if(length(ii)>0) { colnames(ans)[ii] <- paste(colnames(ans)[ii]," <list>") }
  }
  ans
}


#internal
makePrettyMatrixForCompactPrinting2 <- function (x, makeNakedMat.FUN,head.tail=6,up.to=50, show.promo=TRUE) 
{
  lx <- NROW(x)
  if(lx <= up.to) { head.tail <- up.to }
  nhead <- head.tail
  ntail <- head.tail
  if (lx < (nhead + ntail + 1L)) {
    ans <- makeNakedMat.FUN(x,show.promo=show.promo)
    ans_rownames <- .rownames3(names(x), lx)
  }
  else {
    top_idx <- 1:nhead
    if (nhead == 0) 
      top_idx <- 0
    bottom_idx = (lx - ntail + 1L):lx
    if (ntail == 0) 
      bottom_idx <- 0
    ans_top <- makeNakedMat.FUN(x[top_idx, , drop = FALSE],show.promo=show.promo)
    ans_bottom <- makeNakedMat.FUN(x[bottom_idx, , drop = FALSE],show.promo=show.promo)
    ans <- rbind(ans_top, matrix(rep.int("...", ncol(ans_top)), 
                                 nrow = 1L), ans_bottom)
    ans_rownames <- .rownames3(names(x), lx, top_idx, bottom_idx)
  }
  rownames(ans) <- format(ans_rownames, justify = "right")
  ans
}

#internal
.rownames3 <- function (names = NULL, len = NULL, tindex = NULL, bindex = NULL) 
{
  if (is.null(tindex) && is.null(bindex)) {
    if (len == 0L) 
      character(0)
    else if (is.null(names)) 
      paste0("[", seq_len(len), "]")
    else names
  }
  else {
    if (!is.null(names)) {
      c(names[tindex], "...", names[bindex])
    }
    else {
      s1 <- paste0("[", tindex, "]")
      s2 <- paste0("[", bindex, "]")
      if (all(tindex == 0)) 
        s1 <- character(0)
      if (all(bindex == 0)) 
        s2 <- character(0)
      c(s1, "...", s2)
    }
  }
}


# internal function to allow flexible input for the build parameter
ucsc.sanitizer <- function(build,allow.multiple=FALSE,show.valid=FALSE) {
  build.alt <- c("hg15","hg20","hg17","hg18","hg19","hg38",17,18,19,20,35,36,37,38,
                 "build35","build36","build37","build38","b35","b36","b37","b38")
  build.new <- c("hg15","hg20",rep(c("hg17","hg18","hg19","hg38"),times=5))
  if(show.valid) { return(cbind(valid=build.alt,mapsTo=build.new)) }
  build <- build.new[match(tolower(build),build.alt)]
  if(any(is.na(build))) { 
    warning("Illegal build parameter '",build[1],"', defaulting to hg18") 
    build[is.na(build)] <- "hg18" 
  }
  if(allow.multiple) {
    return(build)
  } else {
    return(build[1])
  }
}


# internal (from iFunctions)
# MAYBE8
# note that this will preserve only chr,start,end, nothing else including rownames
ranged.to.data.frame <- function(ranged,use.names=TRUE) {
  u <- as(ranged,"data.frame")
  cn <- tolower(colnames(u))
  if(is(ranged)[1]=="RangedData") {
    if("names" %in% cn) { 
      rownames(u) <- u[["names"]]
      u <- u[,-which(cn=="names")]
    } 
    if("space" %in% cn) { colnames(u)[which(cn=="space")] <- "chr" }
  } else {
    if(is(ranged)[1]=="GRanges") {
      if("seqnames" %in% cn) { colnames(u)[which(cn=="seqnames")] <- "chr" }
    } else {
      warning("'ranged' should be RangedData or GRanges, coercion could fail")
    }
  }
  return(u)
}



# from iFunctions internal?
# convenience function to use GRanges
# MAYBE7
data.frame.to.granges <- function(dat,...) {
  return(data.frame.to.ranged(dat=dat,...,GRanges=TRUE))
}

# internal from iFunctions iFunctions
# MAYBE6
## convert any data frame with chr,start,end, or pos data into a RangedData object
# not case sensitive
data.frame.to.ranged <- function(dat,ids=NULL,start="start",end="end",width=NULL,
                                 chr="chr",exclude=NULL,build=NULL,GRanges=FALSE) 
{
  ## abandon longer names as they clash with function names
  st <- paste(start); en <- paste(end); ch <- paste(chr); wd <- paste(width)
  if((!chr %in% colnames(dat)) & ("seqnames" %in% colnames(dat)) & GRanges) { ch <- "seqnames" }
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  must.use.package(c("genoset","IRanges"),T)
  g.illegal <- tolower(c("seqnames", "ranges", "strand", "seqlevels", "seqlengths", "isCircular", "start", "end", "width", "element"))
  if(is.matrix(dat)) { dat <- as.data.frame(dat,stringsAsFactors=FALSE) }
  if(!is.data.frame(dat)) { stop("Error: not a dataframe")}
  key.nms <- c(ids,st,en,ch,wd)
  tries <- 0
  #print(key.nms); print(colnames(dat))
  while(!all(key.nms %in% colnames(dat))) { 
    colnames(dat) <- tolower(colnames(dat)); key.nms <- tolower(key.nms)
    st <- tolower(st); en <- tolower(en); ch <- tolower(ch); wd <- tolower(wd)
    if(tries>2) {
      if((tolower(st)=="pos" | tolower(en)=="pos") & !(tolower(st)=="pos" & tolower(en)=="pos")) {
        st <- en <- "pos"
      } else {
        if(tolower(st)=="start" & tolower(en)=="end") { st <- en <- "pos" }
      }
    }
    key.nms <- c(ids,st,en,ch,wd)
    tries <- tries+1
    if(tries > 3) { if(!all(c(st,en,ch) %in% colnames(dat))) {
      warning("chromosome and position columns not found") } ; break }
  }
  if(!is.null(ids)) { 
    if(anyDuplicated(dat[[ids]])==0) { 
      id <- dat[[ids]] 
    } else { 
      key.nms <- key.nms[-match(ids,key.nms)] # allow non-unique ids as regular
      ids <- NULL
      warning("id must be unique to form rownames, will insert as a separate column") 
    }
  }
  if(is.null(ids)) { 
    if(!is.null(rownames(dat)) & all(rownames(dat)!=paste(1:nrow(dat)))) { 
      id <- rownames(dat)
    } else { 
      id <- paste(1:nrow(dat)) 
    }
  }
  ## not sure why here are adding 'chr' to X and Y?
  #this was here before? :  if(length(ch)>0) { ch1 <- gsub("Y","chrY",gsub("X","chrX",gsub("chr","",dat[[ch]],ignore.case=T))) } else { ch1 <- NULL }
  if(length(ch)>0) { ch1 <- gsub("chr","",dat[[ch]],ignore.case=T) } else { ch1 <- NULL }
  if(length(st)>0) { st1 <- as.numeric(dat[[st]]) } else { st1 <- NULL }
  if(length(en)>0) { en1 <- as.numeric(dat[[en]]) } else { en1 <- NULL }
  if(length(wd)>0) { en1 <- st1+as.numeric(dat[[wd]]) } # { en1 <- st1+dat[[wd]] }
  #print(length(st1)); print(length(en1)); print(length(id)); print(length(ch1))
  outData <- GRanges(ranges=IRanges(start=st1,end=en1,names=id),seqnames=ch1); genome(outData) <- build[1]
  #outData <- RangedData(ranges=IRanges(start=st1,end=en1,names=id),space=ch1,universe=build[1])
  ###  ###  ###  outData <- toGenomeOrder2(outData,strict=T)
  # note when adding data subsequently that 'RangedData' sorts by genome order, so need
  # to re-sort any new data before adding.
  if(is.null(rownames(outData))) { rownames(outData) <- paste(1:nrow(outData)) }
  reorder <- match(rownames(outData),id)
  more.cols <- colnames(dat)[!colnames(dat) %in% key.nms]
  more.cols <- more.cols[!more.cols %in% exclude]
  if(is(outData)[1]=="GRanges") { more.cols <- more.cols[!more.cols %in% g.illegal] }
  if(length(more.cols)>0) {
    for (cc in 1:length(more.cols)) {
      u <- dat[[more.cols[cc]]][reorder]; #prv(u)
      if(is(outData)[1]=="GRanges") {
        mcols(outData)[[more.cols[cc]]] <- u
      } else {
        outData[[more.cols[cc]]] <- u
      }
    }
  }
  if(GRanges) {
    return(as(outData,"GRanges"))
  } else {
    cncn <- colnames(mcols(outData))
    outData <- as(outData,"RangedData")
    if(any(cncn %in% "strand")) {
      outData <- outData[,-which(cncn=="strand")]
    }
    #outData <- toGenomeOrder2(outData,strict=T)
    return(outData)
  }
}
