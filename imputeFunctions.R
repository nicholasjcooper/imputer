options(ucsc="hg19")

# internal function
nam <- function(x) { y <- narm(x); return(y[y!="NA"]) }

##### contains functions #####
#parse.sexcheck - read a plink sexcheck file and derive the list of samples to exclude
# alphabetize.alleles - Set all alleles to alphabetical order by flipping the strand if required 
# join.alleles - Convert two separate allele columns into a single allele vector
# sep.alleles - Separate allele codes from a single field into two columns
# flip.strand - Flip the strand for allele codes
# dup.pairs - Obtain an index of all members of values with duplicates (ordered)
# sync.id.by.pos - workout which snp ids are equivalent because they have the same chr, pos
# get.strands - compare a SnpMatrix and support file to a reference, e.g, hapmap, and determine strands
# prob.major - probability that an allele with the given frequency is the major allele
# duplicate.report - for a RangedData object, report on any multiple listings for the same gene
# quick.mat.format - to quickly extract whether a delimited matrix
# get.type - get data type of a vector, numeric or character, 
# coerce.chrname - coerce a chr list like (chr1, chr2, chr6_extra_bit,chrx) to (1,2,6,X)
##############################


# read a plink sexcheck file and derive the list of samples to exclude
parse.sexcheck <- function(fn, excl=0.5, write=TRUE) {
  sc <- read.table(fn,header=T)
  sc <- sc[sc$STATUS!="OK",]
  sc <- sc[!is.na(sc$F),]
  somefails <- ((sc$PEDSEX==1 & sc$SNPSEX==2) | (sc$PEDSEX==2 & sc$SNPSEX==1))
  morefails <- sc$F>excl
  failers <- paste(sc[somefails | morefails,"IID"])
  if(write & length(failers)>0) {
    out.fn <- cat.path(dirname(fn),basename(fn),suf="fail",ext="txt")
    writeLines(failers,con=out.fn)
    cat("wrote samples failing sexcheck to",out.fn,"\n")
  }
  return(failers)
}



#internal
comma <- function(...) {
  paste(...,collapse=",")
}




#' Set all alleles to alphabetical order by flipping the strand if required 
#' 
#' The default behaviour of the popular GWAS package snpStats is
#' to store import genotypes in alphabetical order. So if the data
#' was imported and constructed using snpStats functions, you would
#' be able to see in your annotation that the reference allele for the
#' 0,1,2 coding was always the allele later in the alphabet. This
#' function helps to flip the strand for any allele pair not already in
#' alphabetical order, which can help synchronise a support or annotation
#' file with your SnpMatrix object. For instance you may have a support file
#' produced by another program and want to synchronise your SnpMatrix annotation.
#' @param a1 the vector of allele1, or else a 3 column matrix of allele1, allele2 and 
#' the reference allele
#' @param a2 the vector of allele2 values, or NULL if using 3 column input for a1
#' @param ref the reference allele code for each genotype, or NULL if using 3 
#' column input in a1
#' @param missing the value to allocate to missing alleles
#' @param check logical, if TRUE, then will run through the whole process again
#' and if nothing changed we'll know that the alphabetization was successful
#' for all genotypes
#' @return a new set of A1,A2 and reference allele values, now alphabetized, as
#' a 3 column data.frame, and A2 should consistently be the reference allele.
#' @export
#' @examples
#' alphabetize.alleles(c("A","T","G"),c("G","C","C"),ref=c("A","C","C"))
alphabetize.alleles <- function(a1, a2=NULL, ref=NULL, missing=NA, check=TRUE) {
  # ensure snpStats compatible alphabetic snp order,
  # by strand flipping any where the reference is the lower letter
  alphabet <- unlist(strsplit("abcdefghijklmnopqrstuvwxyz","",fixed=TRUE))
  if(!is.null(dim(a1))) {
    if(ncol(a1)==3) {
      ref <- a1[,3]
      a2 <- a1[,2]
      a1 <- a1[,1]
    } else {
      stop("a1 can only be a data frame containing 3 columns for a1, a2 and ref, or else a vector alongside a2 and ref vectors")
    }
  }
  if(length(a1)!=length(a2) | length(ref)!=length(a2)) { stop("a1 must have same length as a2 and ref") }
  a1 <- paste(a1); a2 <- paste(a2); ref <- paste(ref)
  if(!is.character(a1) | !is.character(a2) | !is.character(ref)) { stop("all arguments must be character format") }
  all.codes <- c(a1,a2,ref)
  if(any(nchar(nam(all.codes))>1)) { ii <- which(nchar(nam(all.codes))>1); 
                                     stop("a1,a2, ref allele codes must be no more than 1 character long, e.g, ",comma(nam(all.codes)[ii]))}
  pair.alph <- function(ab) { 
    a <- ab[1]; b <- ab[2]
    aa <- match(tolower(a),alphabet)
    bb <- match(tolower(b),alphabet)
    if(is.na(aa))  { aa <- 0 }
    if(is.na(bb))  { bb <- 0 }
    ret <- TRUE # NA?
    if(bb>aa) { ret <- TRUE } 
    if(aa>bb) { ret <- FALSE } 
    return(ret)
  }
  mini <- cbind(a1,a2)
  if(any(ref!=a2)) {
    cat(length(which((ref!=a2))),"did not have the reference allele as a2\n")
    mini[ref!=a2,] <- mini[ref!=a2,c(2,1)]
   # ref[ref!=a2] <- a2[ref!=a2]
  }
  tfs <- apply(mini,1,pair.alph)
  if(any(!tfs)) {
    cat(length(which((!tfs))),"did not have the reference allele as a2\n")
    ref[!tfs] <- flip.strand(ref[!tfs])
    fs <- flip.strand(join.alleles(mini[!tfs,,drop=FALSE]))
    mini[!tfs,] <- as.matrix(sep.alleles(fs,missing=missing))
  }
  out <- cbind(mini,ref);  colnames(out) <- c("A1","A2","REF")
  if(check) {
    ## see whether we acheived the objective for all alleles
    check <- alphabetize.alleles(a1=out, missing=missing, check=FALSE)
    if(all(out[,3]==check[,3])) {
      cat("check passed, all alleles alphabetized and ref alleles are all a2\n")
    } else {
      failers <- length(which(out[,3]!=check[,3]))
      warning("check failed for ",failers," alleles; allele strand switching did not fix structure consistency")
    }
  }
  return(out)
}


#' Convert two separate allele columns into a single allele vector
#' 
#' In genetic datasets allele pairs are often stored as two separate 
#' columns. This function provides a simple and robust function to 
#' join them together, e.g, join "A", "C" to give "A/C", the separator
#' can be customized. Obviously this is a trivial operation but
#' this function is just for convenience, and provides built-in assurance
#' that the codes are all valid.
#' @param a1 can be a vector of allele1 codes (characters), or a matrix/data.frame with two columns
#' @param a2 a vector of allele2 codes (characters), leave as NULL if a1 has 2 columns
#' @param missing character, the code to use for missing alleles
#' @param sep character, the character to put between the allele codes when joining (can be "")
#' @param allow.IUPAC logical, whether to allow ambiguous allele codes, following the IUPAC
#' standard, ie., R means G or A, H means A or C or T, etc.
#' @export
#' @return returns a vector of 2-part genotypes, will replace any invalid characters with 
#' the missing code.
#' @seealso sep.alleles
#' @examples
#' join.alleles(c("A","G","T","C","-"), c("C","A","A","-","-"))
#' join.alleles(c("A","G","T","C","K"), c("C","A","A","W","N"))
#' join.alleles(c("A","G","T","C","K"), c("C","A","A","W","N"),sep="|", allow.IUPAC=FALSE)
join.alleles <- function(a1,a2=NULL,missing="-", sep="/", allow.IUPAC=TRUE) {
  std <- c("A","T","G","C"); std <- c(toupper(std),tolower(std))
  iupac <- c("R","Y","M","K","S","W","H","B","V","D","N")
  iupac <- c(toupper(iupac),tolower(iupac))
  if(allow.IUPAC) { valid <- c(std,iupac) } else { valid <- std }
  if(!is.null(dim(a1)) & is.null(a2)) {
    if(ncol(a1)==2) { a2 <- a1[,2]; a1 <- a1[,1] } else { stop("'a1' must have 2 columns for matrix input") }
  }
  a1 <- paste(a1) ; a2 <- paste(a2)
  if(!is.character(a1)) { stop("a1 must be coercible to character type")}
  if(!is.character(a2)) { stop("a2 must be coercible to character type")}
  all.codes <- c(a1,a2)
  if(any(nchar(nam(all.codes))>1)) { ii <- which(nchar(nam(all.codes))>1); 
                  stop("a1,a2 codes must be no more than 1 character long, e.g, ",comma(nam(all.codes)[ii]))}
  if(!any(all.codes %in% valid)) { stop("could not find any valid allele codes, is 'allow.IUPAC' set correctly?") }
  if(!allow.IUPAC) { if(any(all.codes %in% iupac)) {
      warning("'allow.IUPAC' is FALSE, so R,Y,M,K, etc, will be set as missing (",missing,")" ) } }
  a1[!a1 %in% valid] <- missing
  a2[!a2 %in% valid] <- missing
  if(length(a1)!=length(a2)) { stop("a1 must be the same length as a2") }
  joined <- paste(a1,a2,sep=sep)
  return(joined)
}

#' Separate allele codes from a single field into two columns
#' 
#' Allele codes are often coded into a single vector containing a separator,
#' e.g, A/C, G/T, AC, A>C. This function will separate these into two
#' separate columns, automatically detecting what the separator is, as long
#' as it is one of: /|-.,:;><
#' @param genotype character, the genotypes to separate into 2 separate columns
#' of allele codes, e.g, A/C, G/T, AC, A>C, etc
#' @param missing character, the code to assign to missing or illegal values
#' in the output data.frame
#' @return As long as one of the supported separators is used, it will be
#' automatically detected and the alleles will be split into two columns,
#' and returned as a data.frame with column names, "A1" and "A2".
#' @seealso join.alleles
#' @examples
#' sep.alleles(c("A/C","G/T","C/A","T/-","G/-"))
#' sep.alleles(c("A>C","G>T","C<A","T","G"))
#' sep.alleles(c("AC","gT","cA","Ta","G"))
#' sep.alleles(c("AC","ZZ","XY","zxxzx","shame","xzxx","xzzz")) # fails!
sep.alleles <- function(genotype,missing=NA) {
  std <- c("A","T","G","C"); std <- c(toupper(std),tolower(std))
  iupac <- c("R","Y","M","K","S","W","H","B","V","D","N")
  delims <- c("","/","\\","|","-",".",",",":",";",">","<")
  max.check <- 100
  iupac <- c(toupper(iupac),tolower(iupac))
  valid <- c(std,iupac)
  genotype <- paste(genotype)
  genotype <- gsub(" ","",genotype)
  if(!is.character(genotype)) { stop("genotype must be coercible to character type") }
  num.del <- list()
  genotype <- gsub(">","<",genotype)
  for (cc in 1:length(delims)) {
    fff <- nchar(delims[[cc]]) == 1
    num.del[[cc]] <- sapply(strsplit(head(genotype,max.check), delims[[cc]], 
                                     fixed = TRUE), length)
  }
  if (all(unlist(num.del) == 1)) {
    stop("no separator was found in text 'genotype' ; allowed =",delims)
  }
  n2s <- sapply(num.del,function(X) { length(which(X==2)) })
  bads <- sapply(num.del,function(X) { length(which(X>2)) })
  maxn <- max(n2s,na.rm=TRUE)
  theone <- which(n2s==maxn)[1]
  #prv(bads,maxn,theone,n2s,num.del)
  if(bads[theone]>n2s[theone]) { stop("no separator was found in text 'genotype' ('sep') that could split the majority of entries into just 2 pairs of allele codes") }
  sep <- delims[theone]
  sepz <- grep(sep,genotype,fixed=TRUE)
  nn <- length(genotype)
  if(maxn < (0.5*length(num.del[[1]]))) { warning("less than half of genotypes contained a separator (guessed it was: '",sep,"')") }
  if(length(sepz)==0) { stop("no separator was found in text 'genotype' ('sep')") }
  all <- strsplit(genotype,split=sep,fixed=TRUE)
  A1 <- sapply(all, "[", 1)
  A2 <- sapply(all, "[", 2)
  A1[!A1 %in% valid] <- missing
  A2[!A2 %in% valid] <- missing
  return(data.frame(A1=A1,A2=A2))
}
  
  


# put into NCmisc
#dup.pairs <- function(x,pairs.only=FALSE) {


#internal
#'
#' @examples
#' ichip.ids <- sync.id.by.pos(rr,sup.for.impute,FALSE)
#' ref.ids <- sync.id.by.pos(rr,sup.for.impute,TRUE)
sync.id.by.pos <- function(ref.snp,snp.inf,return.ref.ids=FALSE) {
  if(!"position" %in% colnames(ref.snp)) { stop("ref.snp must have column 'position'") }
  if(!is(snp.inf)[1] %in% "ChipInfo") { stop("snp.inf must be a ChipInfo object") }
  posz <- startSnp.inf)
  rsz <- rownames(snp.inf)
  rsv <- rmv.trail(id.to.rs(rsz))
  icv <- rs.to.id(rsz)
  fic.a <- find.id.col(ref.snp,ids=rsv)
  fic.b <- find.id.col(ref.snp,ids=icv)
  if(fic.b$maxpc>=fic.a$maxpc) { 
     col <- fic.b$col;  rsz <- icv
  } else {
     col <- fic.a$col;  rsz <- rsv
  }
#   a.pos.match <- match(posz,ref.snp$position)
#   b.pos.match <- match(id.torsz,ref.snp$position)
  if(!return.ref.ids) {
    # return the source ids (e.g, ids from snp.inf)
    pos.match <- match(ref.snp$position,posz) 
    badz <- posz[is.na(pos.match)]
    ichip.ids <- rsz[pos.match]
    one <- "snp.inf"; two <- "ref.snp"
  } else {
    # return the reference ids (e.g, ids from ref.snp)
    pos.match <- match(posz,ref.snp$position)  
    badz <- ref.snp$position[is.na(pos.match)]
    ichip.ids <- ref.snp[pos.match,col]
    one <- "ref.snp"; two <- "snp.inf"
  }
  if(any(is.na(pos.match))) { warning(one," should contain every position found in ",two," :",
                              "could not find ",length(badz),"; positions:",comma(head(badz,50))," ...") }

  return(ichip.ids)
}

# compare a SnpMatrix and support file to a reference, e.g, hapmap, and determine strands
# snp.inf can be a dataframe, but only if it has columns chr, pos/start+end, A1, A2, rs.id
get.strands <- function(ref.snp, snp.mat, snp.inf, aaf.col="eur.aaf", maf.col="eur.maf", id.col="id") {
  ####### extensive checking for valid input! ########
  if(!is(snp.mat)[1]=="SnpMatrix") { stop("snp.mat must be a SnpMatrix") }
  if(!is(ref.snp)[1] %in% c("data.frame")) { stop("ref.snp must be a data.frame") }
  ccols <- c("position","a0","a1",aaf.col,maf.col)
  #print(colnames(ref.snp)); print(ccols)
  if(!all(ccols %in% colnames(ref.snp))) { stop("ref.snp must contain columns: ", paste(ccols,collapse=","),
                                           ", suggest taking it from the *.legend files downloaded to use IMPUTE2, e.g, \n",
                                           "https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference") }
  mm <- match(id.col,colnames(ref.snp))
  if(!is.na(mm)) { ref.ids <- ref.snp[[mm]] } else { ref.ids <- rownames(ref.snp) }
  if(!is(snp.inf)[1] %in% c("data.frame","ChipInfo")) { stop("snp.inf must be a data.frame or ChipInfo object") }
  if(is(snp.inf)[1] == "data.frame") { 
    ## allow data.frame but only if it will smoothly convert to a ChipInfo object:#'
    if(!(all(c("CHR", "A1","A2") %in% toupper(colnames(snp.inf))))) { stop("snp.inf must have columns chr, A1, A2") }
    if(!("POS" %in% toupper(colnames(snp.inf))) & !(all(c("START","END") %in% toupper(colnames(snp.inf))))) { stop("snp.inf must have columns pos, or start + end") }  
    if(is.null(rownames(snp.inf))) { stop("snp.inf rownames should be SNP ids") }
    if(any(duplicated(rownames(snp.inf)))) { stop("snp.inf rownames should be unique (should be SNP ids)") }
    snp.inf <- as(snp.inf,"ChipInfo") 
  }
  if(!all(colnames(snp.mat) %in% rownames(snp.inf))) {
    stop("snp.inf must have a rowname corresponding to every colname of snp.mat [e.g, must cover all SNP-ids]")}
  ref.snp[["chip.id"]] <- sync.id.by.pos(ref.snp,snp.inf)
  if(!all(colnames(snp.mat) %in% ref.ids)) {
    stop("snp.ref must have a rowname/id.col corresponding to every colname of snp.mat [e.g, must cover all SNP-ids]")}
  #######
  cs <- col.summary(snp.mat)
  nsamp <- nrow(snp.mat)
  p.maj.sm <- prob.major(p=1-cs$MAF,n=cs$Calls)
  p.maj.rf <- prob.major(p=1-ref.snp[,maf.col],n=1092)
  conf <- (1-p.maj.sm)*(1-p.maj.rf) # per snp probability of matching on maj.vs.min
  return(conf)
}


prob.major <- function(p,n) {
  z=(p-.5)/(sqrt((p*(1-p))/n)); 
  return(Z.to.p(z))
}


duplicate.report <- function(ga,full.listing=F,colname="ids",pos="position") {
  # for a RangedData object, report on any multiple listings for the same gene
  if(colname %in% colnames(ga)) { 
    gene.col <- which(colnames(ga)==colname) 
  } else { 
    stop("'colname' not found in ga") 
  } 
  if(length(gene.col)>0) { gene.col <- gene.col[1] } else { warning("no 'gene' column"); return(NULL) }
  colnames(ga)[gene.col] <- "gene" #force this colname
  duplicate.report <- T  ### when would this be FALSE???
  if(duplicate.report) {
    culprits <- unique(ga$gene[which(duplicated(ga$gene))])
    n.gene.multi.row <- length(culprits)
    culprit.ranges <- ga[ga$gene %in% culprits,]
    total.culprit.rows <- nrow(culprit.ranges)
    start.same.ct <- 0; which.ss <- NULL
    for (cc in 1:length(culprits)) { 
      mini <- (ga[ga$gene %in% culprits[cc],]) 
      if(full.listing) {
        cat(colname,":",culprits[cc],"# same pos:",anyDuplicated((mini[,pos])),"\n") }
      start.same.ct <- start.same.ct+anyDuplicated((mini[,pos]))
      if(anyDuplicated((mini[,pos])) ) { which.ss <- c(which.ss,cc) }
    }
    cat(" ",colname,"s with split ranges:\n"); print(culprits,quote=F); cat("\n")
    cat(" ",colname,"s with same pos:\n"); print(culprits[which.ss],quote=F); cat("\n")
    cat(" total ",colname,"-segments with same pos",start.same.ct,"\n")
  }
  return(culprits)
}


# internal? function to quickly extract whether a delimited matrix
# has row or column names, and what the column names are
quick.mat.format <- function(fn) {
  temp.fn <- "twmpwerw123.txt"
  txtx <- readLines(fn,n=11) # change this to like 5 when reader is updated to save memory
  writeLines(txtx,con=temp.fn)
  first2 <- reader(temp.fn)
  unlink(temp.fn) # remove this temporary file straight away
  typz <- apply(first2,2,get.type)
  if(is.null(dim(first2))) {
    ## assuming no header row #
    return(list(rownames=F,colnames=F,ncol=1,cnames=NULL,ctypes=typz))
  }
  rn <- rownames(first2); if(rn[1]=="1") { rn <- NULL }
  if(!is.null(rn)) { ern <- T } else { ern <- F }
  cn <- colnames(first2); if(cn[1]=="V1") { cn <- NULL }
  if(!is.null(cn)) { ecn <- T } else { ecn <- F }
  ncl <- ncol(first2)
  return(list(rownames=ern,colnames=ecn,ncol=ncl,cnames=cn,ctypes=typz))
}


#internal 
# get data type of a vector, numeric or character, 
# or if detect.integers=TRUE, then search also for integer type
# when column has less than def.pc% valid numbers, but more than 0%, fall back 
# to default type 'def'.
# when column has > def.pc% valid numbers, then will call numeric
get.type <- function(X,def.pc=0.5,def="character",detect.integers=FALSE) {
  nm <- suppressWarnings(as.numeric(paste(X))) 
  def.pc <- force.percentage(def.pc)
  if(all(is.na(nm))) { 
    typ <- "character"
  } else {
    if(length(which(is.na(nm)))>(def.pc*length(X))) { typ <- def } else { typ <- "numeric" }
  }
  if(typ=="numeric" & detect.integers) {
    if(all(nm==round(nm))) { typ <- "integer" }
  }
  return(typ)
}



#fn <- "/chiswick/data/ncooper/imputation/HiC/monocyte_upto1M_out.txt"
#hic <- read.hic("/chiswick/data/ncooper/imputation/HiC/monocyte_upto1M_out.txt",FALSE)
#myHiC2 <- read.hic("/chiswick/data/ncooper/imputation/HiC/monocyte_upto1M_out.txt")
read.hic <- function(fn, return.S4=TRUE, build=37) {
  hic <- read.table(fn,header=F)
  nr <- nrow(hic)
  if(nr<2) { warning("empty hiC file"); return(NULL) }
  evs <- 2*(1:round(nr/2))
  odds <- (1:nr)[!(1:nr) %in% evs]
  if(length(odds)!=length(evs)) { stop("file should contain an equal number of odd and even lines") }
  oddlines <- hic[odds,]
  evlines <- hic[evs,]
  all <- cbind(oddlines[,1:6],evlines[,1:4])
  colnames(all) <- c("chr","start.bait","end.bait","promoters.bait","N.reads","score","chr.end",
                     "start.end","end.end","promoters.end")
  rownames(all) <- paste(1:nrow(all))
  hic <- all[,c(1,2,3,7,8,9,5,6,4,10)]
  #print(head(hic))
  if(return.S4) {
    gr <- makeGRanges(hic$chr,start=hic[[2]],end=hic[[3]])
    gr2 <- makeGRanges(hic$chr.end,start=hic[[5]],end=hic[[6]])
    #colnames(hic) <- make.names(colnames(hic))
    source("~/github/imputer/hiClass.R")
    hic <- HiC(gr,gr2,hic$N.reads,hic$score,build,hic[[9]],hic[[10]])
  }
  return(hic)    
}



#cd4,eryth,macro,mega,monoc
# jj <- which(genes.bait(cd4)=="CTSH")
# jj0 <- jj
# jj1 <- which(genes.bait(eryth)=="CTSH")
# jj2 <- which(genes.bait(macro)=="CTSH")
# jj3 <- which(genes.bait(mega)=="CTSH")
# jj4 <- which(genes.bait(monoc)=="CTSH")
# print(cd4[jj0,])
# print(eryth[jj1,])
# print(macro[jj2,])
# print(monoc[jj4,])

# snp.list overlap with tissue HiC set (ends)
#for big list

snp.end.overlap <- function(snps, hic, index.only=FALSE, combine=TRUE, chr.only=NULL, range.only=NULL, bait=FALSE) {
  if(!is.null(rownames(snps))) { rn <- rownames(snps); nonames <- FALSE  } else { nonames <- TRUE }
  rownames(snps) <- init <- paste(1:nrow(snps))
  if(bait) { HIC <- bait(hic) } else { HIC <- otherEnd(hic) }
  if(!is.null(range.only)) { range.only <- as.numeric(range.only) ; if(length(narm(range.only))!=2) { stop("invalid range, should be integer length 2") } }
  if(!is.null(chr.only)) {
    chr.only <- chr.only[1] # only allowed 1
    if(!all(paste(chr.only) %in% chrNames(snps))) { stop("invalid chromosome entered, ",chr.only," not found in 'snps': ",comma(chrNames(snps))) }
    if(is.null(range.only)) {
      if(nchar(paste(chr.only))<4) { ch <- paste0("chr",chr.only) } else { ch <- chr.only }
      range.only <- c(1,(get.chr.lens(mito=TRUE,names=TRUE)[ch]))
    } 
    range.limit <- makeGRanges(chr=chr.only,start=range.only[1],end=range.only[2])
    #prv(rownames(snps))
    snps <- subsetByOverlaps(snps,range.limit)
    #print(as.numeric(rownames(snps)))
    rn <- (rn[as.numeric(rownames(snps))])
  #  rn <- rn[rownames(snps) %in% init]
  }  
  ff <- findOverlaps(HIC,snps)
  ss <- subjectHits(ff)
  qq <- queryHits(ff)
  subs <- unique(ss)
  res.list <- lapply(subs, function(x) { qq[which(ss %in% x)] } )
  names(res.list) <- paste(subs); res.list <- res.list[order(as.numeric(names(res.list)))]
  if(combine) {
    index <- unique(unlist(res.list))
    if(index.only) { return(index) } else { return(hic[index,]) }
  } else {
    if(!nonames) { nms <- rn[as.numeric(names(res.list))] ; names(res.list) <- nms }
    if(index.only) { 
      return(res.list) 
    } else { 
      out <- lapply(1:length(res.list),function(x) { hic[res.list[[x]]] }) 
      names(out) <- names(res.list)
      return(out)
    }
  }
}




### HERE ##
if(F) {

fmsnps <- read.table("/chiswick/data/ncooper/imputation/T1D/FMSNPS.csv",stringsAsFactors=FALSE)
colnames(fmsnps) <- c("chr","start","end","width","strand","band","p","rs.id")
nz <- which(fmsnps$rs.id=="NONE")
fmsnps <- fmsnps[-nz,]
fmSnps <- data.frame.to.GRanges(fmsnps)
pv <- mcols(fmSnps)[,"p"]
  
fmSnps2 <- fmSnps[!duplicated(start(fmSnps)),]
nms <- paste(fmsnps$chr,fmsnps$start,sep="_")
rs.list <- tapply(paste(fmsnps$rs.id),factor(nms),c)
nms2 <- paste(chr(fmSnps2),start(fmSnps2),sep="_")
rownames(fmSnps2) <- nms2
rs.flat <- lapply(rs.list,paste,collapse=",")
rs.all <- unlist(rs.list)

print(load("/chiswick/data/ncooper/imputation/COMMON/SNPREFWorking/lookupTableAllRSids.RData"))
index <- match(rs.all,bigRS$rs.id)
result <- with(bigRS,cbind(coerce.chrname(chr[index]),pos[index]))
rownames(result) <- rs.all
colnames(result) <- c("chr","pos")
result[is.na(result[,1]),2] <- as.numeric(Pos(paste(rs.all[is.na(result[,1])])))
result[is.na(result[,1]),1] <- paste(Chr(paste(rs.all[is.na(result[,1])])))
extras <- reader("~/Downloads/markers.csv") # from https://mart.immunobase.org/martanalysis/#!/Markers/maker_id_upload?
result[is.na(result[,1]),1] <- paste(extras[[1]][match(rs.all[is.na(result[,1])],rownames(extras))])
result[is.na(result[,1]),2] <- paste(extras[[2]][match(rs.all[is.na(result[,1])],rownames(extras))])
# last three manually extracted from immunobase
result[c("rs61281301","rs201801601","rs150652635","rs57151617","rs35072264","rs60969753"),1] <- c(12,17,3,17,7,7)
result[c("rs61281301","rs201801601","rs150652635","rs57151617","rs35072264","rs60969753"),2] <- c(54755110,38774852,46345342,35286803,50324168,50295306)
extras <- reader("~/Downloads/markers2.csv") # from https://mart.immunobase.org/martanalysis/#!/Markers/maker_id_upload?
result[is.na(result[,2]),1] <- paste(extras[[1]][match(rs.all[is.na(result[,2])],rownames(extras))])
result[is.na(result[,2]),2] <- paste(extras[[2]][match(rs.all[is.na(result[,2])],rownames(extras))])
gresult <- data.frame.to.GRanges(result,build=37) # last one is 'NONE'
rm(bigRS) # free up the memory
####### READ IN THE HiC DATA #######
cd4 <- read.hic("/chiswick/data/ncooper/imputation/HiC/cd4_upto1M_out.txt")
eryth <- read.hic("/chiswick/data/ncooper/imputation/HiC/erythroblast_upto1M_out.txt")
macro <- read.hic("/chiswick/data/ncooper/imputation/HiC/macrophage_upto1M_out.txt")
mega <- read.hic("/chiswick/data/ncooper/imputation/HiC/megakaryocyte_upto1M_out.txt")
monoc <- read.hic("/chiswick/data/ncooper/imputation/HiC/monocyte_upto1M_out.txt")
ii1 <- (snp.end.overlap(fmSnps2,cd4,combine=FALSE))
ii2 <- (snp.end.overlap(fmSnps2,eryth,combine=FALSE))
ii3 <- (snp.end.overlap(fmSnps2,macro,combine=FALSE))
ii4 <- (snp.end.overlap(fmSnps2,mega,combine=FALSE))
ii5 <- (snp.end.overlap(fmSnps2,monoc,combine=FALSE))
All.Tissues <- list(CD4=ii1,Erythrocytes=ii2,Macrophages=ii3,Megakaryocytes=ii4,Monocytes=ii5)
hic.list <- list(CD4=cd4,Erythrocytes=eryth,Macrophages=macro,Megakaryocytes=mega,Monocytes=monoc)

#plot.across.tissues(fmSnps2, gresult, hic.list, fn="tissuePlot", bait=FALSE, score.as.y=TRUE,snp.pch=17, ylim=c(0,20))
#plot.across.tissues(fmSnps2, gresult, list(CD4=cd4), fn="FulltissuePlotHiCChevronsWideCD4", bait=FALSE, score.as.y=TRUE,
#                    snp.pch=17, ylim=c(0,20), end.overlap=FALSE, plot.hic=TRUE, exp.pc=1.5, pdf.width=25)
#cdT <- read.hic("/chiswick/data/ncooper/imputation/HiC/TEST_cd4_otherCHR.txt")
#plot.hic(cdT[[6]],fn="cdTChr1Test",ylim=c(0,22),lty.join=c("dashed","dotted"))
#cd4_un <- cd4[!duplicated(start(cd4)),]
#cd4_un_ov <- subsetByOverlaps(cd4_un,fmSnps2)
#plot.across.tissues(cd4_un_ov, gresult, list(CD4=cd4), fn="BaitWisetissuePlotHiCChevronsWideCD4", bait=FALSE, score.as.y=TRUE,
#                    snp.pch=17, ylim=c(0,20), end.overlap=FALSE, plot.hic=TRUE, exp.pc=1.5, pdf.width=25)
#plot.across.tissues(bait(cd4_un_ov), gresult, list(CD4=cd4), fn="BaitWisetissuePlotHiCChevronsWideCD4", bait=FALSE, score.as.y=TRUE,
#                    snp.pch=17, ylim=c(0,20), end.overlap=FALSE, plot.hic=TRUE, exp.pc=.9, pdf.width=25, fit.both=TRUE)
#plot.across.tissues(bait(cd4_un_ov), gresult, list(CD4=cd4), fn="BaitWisetissue45PlotHiCChevronsWideCD4", bait=FALSE,
#                    score.as.y=TRUE, deg45=TRUE,snp.pch=17, ylim=c(0,20), end.overlap=FALSE, 
#                    plot.hic=TRUE, exp.pc=.9, pdf.width=25, fit.both=TRUE)

}

if(F) {

## check just difs
cur.loc <- "/chiswick/data/ncooper/imputation/COMBINED/sst-aligned"
ichip.t1d <- reader("/chiswick/data/ncooper/iChipData/compiledTableAllResultsPassingQC5.RData")
ibase <- reader("/home/ncooper/Downloads/IC_RESULTS_15_02_2013.tab")
ichip.ra <- ibase[,"ra_eyre",drop=F]
ichip.atd <- ibase[,"cooper_aitd",drop=F]
ichip.celiac <- ibase[,"celiac_trynka",drop=F]
#ichip.rsids <- rownames(ichip.ra)
#mmrt <- ibase$marker_mart
#rs.maybe <- paste0("rs",abs(mmrt))
#rss <-(rs.to.id(rs.maybe))
#ic.rsids <- character(); for (cc in 1:10) { ic.rsids[cc] <- substr(ichip.rsids[cc],1,nchar(ichip.rsids[cc])-17) } #slow
#notgot <- (which(!ic.rsids %in% rownames(chip.support())))
#ic.rsids[notgot] <- rss[notgot]
(load("ic.rsids.backup.RData")) # instead of commented out lines above

options(ucsc="hg19")
ic.chr.pos <- paste(Chr(ic.rsids),Pos(ic.rsids)) # e.g, 1_21345434 - unique and universal snp id!
ic.t1d.chr.pos <- paste(ichip.t1d$Chr,ichip.t1d$Pos)

dis.codes <- c("t1d","ra","atd","celiac")

all.files <- cat.path(cur.loc,list.files(cur.loc))
corz.e <- numeric(length(all.files))
corz <- list(corz.e,corz.e,corz.e,corz.e,corz.e,corz.e)[1:length(dis.codes)]
difzs.e <- vector("list",length(all.files))
names(difzs.e) <- basename(all.files)

difzs <- list(difzs.e,difzs.e,difzs.e,difzs.e,difzs.e,difzs.e,difzs.e)[1:length(dis.codes)]
rps <- ips <- vector("list",length(dis.codes)); names(rps) <- dis.codes; names(ips) <- dis.codes
names(difzs) <- names(corz) <- dis.codes
RESULTS <- NULL
for (dd in 1:length(all.files)) {
  (load(all.files[dd]))
  print(basename(all.files[dd]))
  RESULTS <- rbind(RESULTS,results)
  chris.chr.pos <- paste(results$chromosome,results$position)
  ips[["t1d"]] <- as.numeric(ichip.t1d[match(chris.chr.pos,ic.t1d.chr.pos),"p.value"])
  mm <- match(chris.chr.pos,ic.chr.pos)
  ips[["ra"]] <- as.numeric(ichip.ra[mm,1])
  ips[["atd"]] <- as.numeric(ichip.atd[mm,1])
  ips[["celiac"]] <- as.numeric(ichip.celiac[mm,1])
  
  rps[["t1d"]] <- as.numeric(results$p.t1d)
  rps[["ra"]] <- as.numeric(results$p.ra)
  rps[["atd"]] <- as.numeric(results$p.atd)
  rps[["celiac"]] <- as.numeric(results$p.celiac)
 # rps[["ms"]] <- results$p.ms
 # rps[["jia"]] <- results$p.jia
  for (ee in 1:length(dis.codes)) {
    corz[[ee]][dd] <- (cor((-log10(ips[[ee]])),(-log10(rps[[ee]])),use="pairwise.complete"))
    difz <- -(-log10(ips[[ee]])) + (-log10(rps[[ee]]))
    difzs[[ee]][[dd]] <- (summary(difz))
  }
  
}

}

followup.worst <- function(difzs,dis=1) {
  mins <- sapply(difzs[[dis]],"[",1)
  bad.mins <- head(sort(mins[(mins< -10)]),3)
  if(length(bad.mins)>0) {
   #print(bad.mins); print(mins)
    min.nos <- narm(match(bad.mins,mins))
    for(cc in 1:length(min.nos)) {
      difz <- one.run(dd=min.nos[cc],ee=dis)
      #print(min.nos[cc]); print(bad.mins[cc])
      cat("lower than",round(-(exp(log(abs(bad.mins[cc]))-.3))),"\n")
      #print(summary(difz))
      #print(which(difz < -(exp(log(abs(bad.mins[cc]))-.3))))
      cat(comma(ic.rsids[which(difz < -(exp(log(abs(bad.mins[cc]))-.3)))]),"\n")
      cat("\n")
    }
  } else {
    cat("none lower\n\n")
  }
  
  maxs <- sapply(difzs[[dis]],"[",6)
  bad.maxs <- head(rev(maxs[(maxs>10)]),3)
  if(length(bad.maxs)>0) {
    max.nos <- narm(match(bad.maxs,maxs))
    for(cc in 1:length(max.nos)) {
      difz <- one.run(dd=max.nos[cc],ee=dis)
      #print(dd)
      cat("higher than",round(exp(log(bad.maxs[cc])-.3)),"\n")
      cat(comma(ic.rsids[which(difz > (exp(log(bad.maxs[cc])-.3)))]),"\n")
      cat("\n")
    }
  } else {
    cat("none higher\n\n")
  }
}


one.run <- function(dd,ee) {
  (load(all.files[dd]))
  print(basename(all.files[dd]))
  chris.chr.pos <- paste(results$chromosome,results$position)
  ips[["t1d"]] <- as.numeric(ichip.t1d[match(chris.chr.pos,ic.t1d.chr.pos),"p.value"])
  mm <- match(chris.chr.pos,ic.chr.pos)
  ips[["ra"]] <- as.numeric(ichip.ra[mm,1])
  ips[["atd"]] <- as.numeric(ichip.atd[mm,1])
  ips[["celiac"]] <- as.numeric(ichip.celiac[mm,1])
  
  rps[["t1d"]] <- as.numeric(results$p.t1d)
  rps[["ra"]] <- as.numeric(results$p.ra)
  rps[["atd"]] <- as.numeric(results$p.atd)
  rps[["celiac"]] <- as.numeric(results$p.celiac)
  # rps[["ms"]] <- results$p.ms
  # rps[["jia"]] <- results$p.jia
  corz[[ee]][dd] <- (cor((-log10(ips[[ee]])),(-log10(rps[[ee]])),use="pairwise.complete"))
  difz <- -(-log10(ips[[ee]])) + (-log10(rps[[ee]]))
  difzs[[ee]][[dd]] <- (summary(difz))
  return(difz)
}





#' @param hic 'HiC' object, must only contain bait from 1 chromosome. Can use [[n]] to select only
#' chromosome 'n', which can be the element number or name of the chromosome. Must also have
#' length less than 1000, or if xlim is used, length within the specified window should be less than
#' 1000.
#' @param alt.y what values to use for the y-axis, can be "score" or "reads", or another 
#' variable in the HiC object. Alternatively a vector of the same length as HiC can be passed into
#' this parameter, although be careful with lengths if using in conjunction with 'xlim'.
#' @param exp.pc numeric, expansion percentage, length to expand the graph either side of the range
#' of the bait(s); although will be overridden if 'xlim' is in use
#' @param fn filename for the output pdf, or leave blank to plot to the terminal
#' @param build ucsc build, default is getOption("ucsc"), e.g, 36/hg18, 37/hg19, etc.
#' @param verbose logical, whether to display a message when xlim is in use, and whether to use
#' a progress bar for plots with more than 20 entries
#' @param xlim limit for the x-axis, see plot(), but will also filter the 'hic' 
#' object for only baits in this range
#' @param ylim limit for the y-axis, see plot()
#' @param y.FUN function to scale the y-axis values (e.g, 'score'). Default is 
#' .default.score.fun(), see plot.across.tissues() for more details. Use 'c' (combination function)
#'  to effectively turn this off (or NULL).
#' @param join.scl percentage, vertical scale factor for the chevrons that join the bait and other
#' end ranges, i.e, to control the height
#' @param deg45 logical, if TRUE, make all chevrons join at 45 degree angles, else (FALSE), all
#' chevrons will have the same height, and angle will be determined via the distance spanned.
#' @param plot.new whether to draw the ranges specified in a new plot or add to existing
#' @param col.bait colour for the bait ranges, can be length 1 (all the same), or the length of 
#' hic (filtered), to specify specific colours for each range
#' @param col.end colour for the other-end ranges, same restrictions as col.bait
#' @param col.join colour for the joining chevrons, same restrictions as col.bait
#' @param lwd.ranges line width for bait and other-end ranges
#' @param lty.join line type for the joining chevrons, can use a character vector length 2 to
#' specify a different line type from bait to centre, versus centre to 'other-end'.
plot.hic <- function(hic, alt.y="score", exp.pc=.5, fn=NULL, build=NULL, verbose=TRUE,
                     xlim=NULL, ylim=NULL, y.FUN=.default.score.fun, join.scl=1, deg45=FALSE, plot.new=TRUE,
                     col.bait="black", col.end="red", col.join="black", lwd.ranges=2, lty.join="dotted") {
  if(is(hic)[1]!="HiC") { stop("hic must be an object of type 'HiC', e.g, obtained reading a SeqMonk Hi-C file using read.hic()") }
  main.chr <- unique(chr(hic))
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(is.null(y.FUN)) { y.FUN <- c }
  if(length(main.chr)>1) { stop("cannot plot bait data from multiple chromosomes") }
  if(length(hic)>1000) { stop("hic has maximum length of 1000, set had ",length(hic),". Suggest plotting subsets for large datasets") }
  if(!is.null(xlim)) {
    xlim <- narm(as.numeric(xlim))
    if(length(xlim)!=2) { stop("xlim should be an x-axis range, a numeric vector of length 2")}
    range.limit <- makeGRanges(chr=main.chr,start=xlim[1],end=xlim[2])
    if(is.null(rownames(hic))) { rownames(hic) <- paste(1:length(hic)) }
    rn <- rownames(hic)
    hic.set <- subsetByOverlaps(bait(hic),range.limit)
    rn.sub <- rownames(hic.set)
    to.keep <- match(rn.sub,rn)
    ltk <- length(to.keep)
    if(ltk<1) { 
      warning("no hic entries found on ",ranges.to.txt(range.limit))
      if(plot.new){
        plot(1,1,ylim=ylim, xlim=xlim, col="white", xlab=paste("Chromosome",main.chr),ylab=alt.y[1]) }
      return() 
    }
    hic <- hic[to.keep,]
    if(ltk>0 & verbose) { message(ltk," hic entries found on ",ranges.to.txt(range.limit)) }
  }
  if(length(col.bait)!=1 & length(col.bait)!=length(hic)) { stop("col.bait must have length 1, or the same length as hic (including any filters)") }
  if(length(col.end)!=1 & length(col.end)!=length(hic)) { stop("col.end must have length 1, or the same length as hic (including any filters)") }
  all.chr <- c(chr(hic),chr2(hic))
  if(is.character(fn) & length(fn)==1) { pdf(cat.path("",fn,ext="pdf")) } 
  chrz <- unique(all.chr)
  n.chr <- length(chrz)
#  if(n.chr>2) { stop("cannot plot bait data that transitions across more than 1 other chromosome") }
  all.locs <- c(start(hic[[main.chr]]), end(hic[[main.chr]]), start2(hic[[main.chr]]), end2(hic[[main.chr]]))
  main.rng <- extend.50pc(range(all.locs),Chr=main.chr,snp.info=chip.support(build=build),pc=exp.pc) # extend 50pc is from FunctionsCNVAnalysis.R
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(hic)) {
        yy <- alt.y
        y.column <- substitute(alt.y)
      } else {
        warning("alt.y ignored, must be same length as hic, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        cn <- colnames(mcols(hic)); 
        if(!alt.y %in% cn) { stop("y.axis column name ",alt.y," not found in 'ranged'") }
        vec.y.FUN <- function(X) { out <- unlist(sapply(X,y.FUN)) }
        y.column <- alt.y
        mcols(hic)[,y.column] <- vec.y.FUN(mcols(hic)[,y.column]) # transform 'score' or 'reads'
        df <- mcols(hic)
        yy <- df[,alt.y]; rm(df); 
      } else { 
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(is.null(ylim)) { if(exists("yy")) { ylim <- range(c(0,yy)) } else { yy <- 1:length(hic); y.column <- "" } }
  #prv(ylim,yy,alt.y)
  # create a blank plotting area #
  if(is.null(xlim)) { xlim <- main.rng }
  if(plot.new) {
    plot(1,1,ylim=ylim, xlim=xlim, col="white", xlab=paste("Chromosome",main.chr),ylab=y.column)
  }
  n.sets <- length(hic)
  for (cc in 1:n.sets) {
    p1 <- bait(hic[cc,]); p2 <- otherEnd(hic[cc,])
    plotRanges(p1,alt.y=alt.y,lwd=lwd.ranges, do.labs=F, col=col.bait,skip.plot.new=T)
    c1 <- abs(end(p1)-start(p1))/2 + min(start(p1),end(p1))
    c2 <- abs(end(p2)-start(p2))/2 + min(start(p2),end(p2))
    yoffs <- (diff(ylim)/20)*join.scl
    cp2 <- chrNums(p2); cp1 <- chrNums(p1)
    if(chr(p2)!=chr(p1)) { 
      sw <- if(cp2>cp1) { 1 } else { -1 }
      xoffs <- (diff(xlim)/10)*join.scl*sw
      coords <- rbind(c(c1,yy[cc]),
                      c(c1+xoffs,yy[cc]+yoffs*2),
                      c(c1+(xoffs*1.5),yy[cc]+yoffs*2))
      lines(coords[1:2,],col=col.join,lty=lty.join[1])
      lines(coords[2:3,],col=col.join,lty=lty.join[1])
      arrows(coords[3,1],coords[3,2],coords[3,1]+xoffs/100,coords[3,2],length=.1)
      text(coords[3,,drop=F],labels=paste0(ranges.to.txt(p2)),cex=.5,pos=3+sw)
#      message("for region ",cc,", 'other-end' had a different chromosome to the bait, plotting this is not yet implemented!")
    } else {
      midpt <- abs(c2-c1)/2 + min(c1,c2)
      if(deg45) {
        middst <- abs(c2-midpt)
        yoffs <- middst*(diff(ylim)/diff(xlim)) # distance to make chevrons at fixed 45 degree angle
      }
      plotRanges(p2,alt.y=alt.y,lwd=lwd.ranges, do.labs=F, col=col.end,skip.plot.new=T)
      if(length(lty.join)==1) { lty.join <- rep(lty.join,2) }
      lines(c(c1,midpt),yy[cc]+c(0,yoffs),col=col.join,lty=lty.join[1])
      lines(c(midpt,c2),yy[cc]+c(yoffs,0),col=col.join,lty=lty.join[2])
    }
    if(n.sets>20 & verbose) { loop.tracker(cc,n.sets) }
  }
  ## make dual plot as mini window? or a mfrow thing? or on same plot?
  if(is.character(fn) & length(fn)==1) { dev.off() } 
}


# default function to scale the HiC 'score' for plotting
.default.score.fun <- function(x, scale.all=c, above=10, scale.above=sqrt, max=20) {
  x <- scale.all(x)
  u <- x[x>above]
  lu <- length(u)
  if(lu>0) { 
    u <- u-above; u <- scale.above(u)
    if(length(u)!=lu) { stop("invalid scale.above function, output changed length") } 
    x[x>above] <- above+u
  }
  if(any(x>max)) { x[x>max] <- max }
  return(x)
}

# internal
midpoint <- function(x) {
  u <- range(x,na.rm=T)
  return(min(u) + (max(u)-min(u))/2)
}

#' @param regions GRanges object containing phenotype regions of interest to plot over, essentially
#' an object specifying chr, start, end, and optionally 'band', e.g, 8q23.1, etc. You can create
#' easily from standard vectors using makeGRanges()
#' @param snps GRanges object containing SNPs tested for the study, essentially an object 
#' specifying chr and pos, You can create easily from standard vectors using makeGRanges()
#' @param hic.list list of HiC objects, see 'HiC()', assuming a separate HiC type element for each
#' tissue. These objects can be imported from  Seqmonk browser 
#' files (http://www.bioinformatics.babraham.ac.uk/projects/seqmonk/), or created using HiC().
#' These contain start and end position information for the 'bait' and 'other end' of HiC data.
#' They also contain read counts and scores, and lists of intersected promoters. The names
#' of the elements of hic.list will be used as the labels for each tissue type in the plot
#' @param fn character, name of the output pdf file, if null, will plot to standard graphics
#' @param bait logical, whether to plot the bait regions for each tissue, or the 'other end's.
#' @param end.overlap logical, whether to choose bait-end segments where the end overlaps the
#' 'region' range (TRUE), or where the bait does (set to FALSE)
#' @param fit.both logical, whether to ensure the xlim range covers both ends of bait-end pairs,
#' (set TRUE) or whether just to fit the x axis range to the current 'region' from 'regions,
#' plus the 'exp.pc' expansion factor.
#' @param build character, ucsc build, e.g, 37/hg19 or 36/hg18, to use
#' @param score.as.y character, which HiC column to take the y-axis values from, can be
#' either 'reads' or 'score'. TRUE will use 'score', FALSE will use 'reads'.
#' @param plot.hic logical, whether to plot bait and end joined by a chevron using the
#' 'plot.hic()' function, or to just plot the bait/end as specified.
#' @param deg45 logical, if plot.hic=TRUE, whether to use a fixed angle (TRUE) or a fixed
#' height (FALSE) for the chevrons
#' @param by.bait logical, if TRUE, then a separate plot will be made for each 'bait' rather
#' than for each region ['regions']
#' @param only.if.snp logical, if TRUE, only baits/regions overlapping at least one SNP
#' shall be plotted, if FALSE, all baits intersection a 'region' will be plotted ['regions'].
#' @param tissue.cols character, vector of colour names to use for each tissue type, should
#'  be the same  length as hic.list, or longer (in which case only the first 'length(hic.list)'
#'  colours will be used)
#' @param genes.col character, single string to indicate the colour for the gene annotation,
#' default is grey. Set to NULL to avoid plotting genes.
#' @param reg.col character, single string to indicate the colour for the 'region' line,
#' @param snp.pch see 'pch' in points(); symbol to use to plot SNPs (taken from 'snps' & 
#' overlapping 'regions')
#' @param ylim numeric, range, length 2, see ?plot parameter of the same name. Default is 0 to 20,
#' assuming this is a good range for values of 'score', when scaled using the score.FUN().
#' @param exp.pc percentage (eg, .5 = 50%) to expand plot either side of each 'region'. 
#' @param pdf.width create a super-wide pdf to allow detail to be more easily discernible on the 
#' screen (7 is the default)
#' @param score.FUN function, transformation function for the variable specified in 'score.as.y',
#'  e.g, score.
#' For instance it may have infinite values, or you may wish to log transform, etc. Should take a 
#' parameter x (length 1) and return a single value. The default function will scale anything 
#' above 10 using sqrt() and caps the max at 20. Use NULL or 'c' (combination function) to 
#' effectively turn this off.
#' @examples
#' plot.across.tissues(fmSnps2, gresult, hic.list, fn = "BaitWisetissue45PlotHiCChevronsWide5Tis3", 
#'                       bait = FALSE, score.as.y = TRUE, deg45 = TRUE, snp.pch = 17, 
#'                       ylim = c(0, 20), end.overlap = FALSE, plot.hic = TRUE, exp.pc = 0.9, 
#'                       pdf.width = 25, fit.both = TRUE, by.bait = TRUE)
plot.across.tissues <- function(regions, snps, hic.list, fn=NULL, bait=FALSE, end.overlap=!bait, fit.both=FALSE,
                                build=NULL, score.as.y=TRUE, plot.hic=FALSE, deg45=FALSE, by.bait=FALSE, only.if.snp=TRUE,
                                tissue.cols=c("red","blue","orange","green","brown"), genes.col="grey", 
                                reg.col="black", snp.pch=17, ylim=c(0,20), exp.pc=.9, pdf.width=7,
                                score.FUN=.default.score.fun) {
  if(is(hic.list)[1]=="HiC") { stop("HiC object, rather than a list of HiC objects overlaps was entered") }
  if(length(tissue.cols)<length(hic.list)) { stop("tissue.cols must have the same length as hic.list (or longer)") }
  if(is.null(names(hic.list))) { names(hic.list) <- paste0("tissue.",1:length(hic.list)) ; 
  warning("the list elements of hic.list did not have 'names()', so tissue types will be displayed in the legend simply
          as tissue.1, tissue.2, ..., tissue.n")}
  n.tissue <- length(hic.list)
  tissue.cols <- tissue.cols[1:n.tissue] # take the first 'n' colours
  if(!score.as.y) { y.column <- "reads" } else { y.column <- "score" }
  if(is.character(fn) & length(fn)==1) { pdf(cat.path("",fn,ext="pdf"),width=pdf.width) } 
  if(is.null(build)) { build <- getOption("ucsc") }
  build <- ucsc.sanitizer(build)
  if(bait) { end.extract <- ranges.bait } else { end.extract <- otherEnd }
  gs <- get.gene.annot(build=build) # save loading each time in the loop
  if(is.null(score.FUN)) { score.FUN <- c }
  vec.score.FUN <- function(X) { out <- unlist(sapply(X,score.FUN)) }

  regions.orig <- regions # if not by.bait=TRUE, this is just a copy
  if(by.bait) {
    #### in the middle of implementing THIS #####
    fit.both <- TRUE
    regions <- Baits(as(hic.list,"HiCList"))
#    regions <- hic.list[[1]][!duplicated(start(hic.list[[1]])),]
    ### maybe this? ### regions <- subsetByOverlaps(regions.orig,regions) # original regions overlapped by the baits
    regions <- subsetByOverlaps(regions,regions.orig) # baits overlapped by original regions
    ## not sure if this should go here, or be more general and go after or before ##
    if(only.if.snp) {
      regions <- subsetByOverlaps(regions,snps)
    }
  }
  # just moved these from above the looop...
  All.Tissues <- as(hic.list,"list")
  All.Tissues <- lapply(hic.list, function(x) { snp.end.overlap(regions,x,combine=FALSE,bait=!end.overlap) })
  up.to <- length(regions)
  #return(regions)
  # LOOP ONCE FOR EACH REGION (or bait if that option selected) #
  for(cc in 1:up.to){
    fi <- 1
    nxt.chr <- chr(regions[cc,])
    rng0 <- c(start(regions[cc,]), end(regions[cc,]))
    if(by.bait) { 
      this.region <- subsetByOverlaps(regions.orig,regions[cc])[1,]
      rng2 <- c(start(this.region), end(this.region)) 
      if(all(!is.na(rng2))) {
        rng0 <- range(c(rng0,rng2),na.rm=T)
      }
    }
    rng <- extend.50pc(rng0,Chr=nxt.chr,snp.info=chip.support(build=build),pc=exp.pc) # extend 50pc is from FunctionsCNVAnalysis.R
    hi.ind <- match(rownames(regions)[cc],names(All.Tissues[[fi]]))
    #return(All.Tissues)
    #print(rownames(regions)[cc]); prv(names(All.Tissues[[fi]])); prv(fi)
    while(is.na(hi.ind) & fi<n.tissue) { fi <- fi+1; hi.ind <- match(rownames(regions)[cc],names(All.Tissues[[fi]])) }
    if(!is.na(hi.ind)) {
      ## only run if this index is found ##
      if(fit.both) {
        # make sure xlim is expanded to include maximum bait-end span in the region #
        big.range <- NULL
        for (ee in (fi):n.tissue) {
          #### think i FIXed THIS #####
          hi.ind <- match(rownames(regions)[cc],names(All.Tissues[[ee]]))
          if(!is.na(hi.ind)) {
            AA <- All.Tissues[[ee]] # which is broken? how?
            #print(length(AA)); print(head(AA))
            aa <- AA[[hi.ind]] # which is broken?
            ba <- bait(aa)
            oe <- otherEnd(aa); oe[(chr(aa)!=chr2(aa)),] <- ba[(chr(aa)!=chr2(aa)),]  # i.e, don't allow another chr range into mix
            big.range <- range(c(big.range,start(ba),end(ba),start(oe),end(oe)),na.rm=TRUE)
          }
        }
        if(all(!is.na(big.range))) { rng <- big.range } # set rng to maximum width of any bait/end in the set
      }
      # check for possible errors due to missing/incomplete sets:
      if(is.na(fi) | is.na(hi.ind)) {  next }
      if(length(All.Tissues)<fi) {  next }
      if(length(All.Tissues[[fi]])<hi.ind) {  next }
      #
      hic.rngz <- All.Tissues[[fi]][[hi.ind]]
      if(is.null(hic.rngz)) { next } # if 
      rngz <- end.extract(hic.rngz)
      mcols(rngz)[,y.column] <- vec.score.FUN(mcols(rngz)[,y.column]) # transform 'score' or 'reads'
      ### PLOT A BLANK CANVAS ###
      plot(NA,NA,col="white",ylim=ylim, xlim=rng,xlab=paste("Chromosome",nxt.chr),ylab=y.column)
      ## PLOT BAIT AND END FOR THE FIRST TISSUE ##
#       if(plot.hic) {
#         plot.hic(hic.rngz, alt.y=y.column, exp.pc=exp.pc, fn=NULL, build=build, verbose=FALSE,
#                              xlim=rng, ylim=ylim, y.FUN=c, join.scl=1, plot.new=TRUE, deg45=deg45,
#                              col.bait=tissue.cols[fi], col.end=tissue.cols[fi], col.join="black",
#                              lwd.ranges=2, lty.join=c("dashed","dotted")[2])
#       } else {
#         plotRanges(rngz,alt.y=y.column,lwd=2,
#                   ylim=ylim, xlim=rng,do.labs=F, col=tissue.cols[fi], xlab=paste("Chromosome",nxt.chr),ylab=y.column)
#       }
      if(n.tissue>0 & fi<n.tissue) {
        ## PLOT BAIT AND END FOR EACH TISSUE ##
        for (dd in (fi):n.tissue) {
          hi.ind <- match(rownames(regions)[cc],names(All.Tissues[[dd]]))
          if(!is.na(hi.ind)) {
            hic.rngz <- All.Tissues[[dd]][[hi.ind]]
            rngz <- end.extract(hic.rngz)
            #print(dd)
            mcols(rngz)[,y.column] <- vec.score.FUN(mcols(rngz)[,y.column]) # transform 'score' or 'reads'
            if(plot.hic) {
              plot.hic(hic.rngz, alt.y=y.column, build=build, verbose=FALSE,
                       xlim=rng, y.FUN=c, join.scl=1, plot.new=FALSE, deg45=deg45,
                       col.bait=tissue.cols[dd], col.end=tissue.cols[dd], col.join="black", 
                       lwd.ranges=2, lty.join=c("dashed","dotted")[2])
            } else {
              plotRanges(rngz,alt.y=y.column,lwd=2,
                        do.labs=F, col=tissue.cols[dd],skip.plot.new=T)
            }
          }
        }
      }
    } else {
      print(regions[cc,])
      warning("region ",cc," [",ranges.to.txt(regions[cc,]),"] had no rownames in any tissue of 'hic.list'")
      plot(1,1,col="white",main=paste("region",cc,"was empty"))
      #text(midpoint(rng),midpoint(ylim), labels=paste("region",cc,"was empty"))
    }
    lblz <- if("band" %in% colnames(mcols(regions))) { "band" } else { NULL }
    ### PLOT the REGION ###
   ## print(ranges.to.txt(this.region))
    if(by.bait) {
      plotRanges(this.region,skip.plot.new=T,col=reg.col,lwd=3, labels=lblz, pos=2) # alt.y=rep(2,length(cc))
    } else {
      plotRanges(regions[cc,],skip.plot.new=T,col=reg.col,lwd=3, labels=lblz, pos=2) # alt.y=rep(2,length(cc))
    }
    ### PLOT the overlapping GENES ###
    if(!is.null(genes.col)) {
      plotGeneAnnot(chr=nxt.chr,pos=rng,scl="b",y.ofs=.5,build=build,gs=gs, box.col=genes.col)
    }
    #### PLOT ANY SNPS in RANGE ###
    range.limit <- makeGRanges(chr=nxt.chr,start=rng0[1],end=rng0[2],build=build)
    sub.snps <- subsetByOverlaps(snps,range.limit)
    myseq <- c(start(regions[cc,]),start(toGenomeOrder(sub.snps)),end(regions[cc,]))
    #print(myseq)
    dd <- diff(myseq)
    pcs <- round(dd/diff(range(myseq)),3)
    #print(cc); print(pcs)
    if(any(pcs<.05)) { leg.snps <- TRUE } else { leg.snps <- FALSE } # are any really close together?
    if(nrow(sub.snps)>0) {
      lss <- length(sub.snps); nCl <- 1; snp.PCH <- snp.pch
      if(lss>22) { 
        cx <- .5
        cz <- "darkgrey"; snp.PCH <- 2
        if(lss>40) { nCl <- (lss %/% 40)+1 }
      } else { 
        cx <- .75
        cz <- get.distinct.cols(lss)
      }
      if(!leg.snps) { rownames(sub.snps) <- paste0(". ",rownames(sub.snps)) } # hack to offset the names a little
      plotRanges(sub.snps,skip.plot.new=T,col=cz,
                  pch=snp.pch,srt=90,alt.y=rep(1.5,lss),do.labs=!leg.snps, pos=4) 
      if(leg.snps) { legend("topleft",legend=rownames(sub.snps), col=cz, pch=snp.PCH, cex=cx, bty="n", ncol=nCl) }
    } else { warning("no overlapping snps in ",ranges.to.txt(range.limit),"\n") }
    legend("topright",legend=names(All.Tissues),lwd=2,col=tissue.cols,bty="n")
    #### END PLOT SNPS ####
    loop.tracker(cc,up.to)
  }
  if(is.character(fn) & length(fn)==1) { dev.off() } 
}




## REDUNDANT!!! ?? ##
# removes all but the first of any overlapping ranges within the same sample
# should be a list of cnvs with a column 'id', e.g, sample id
remove.duplicated.ranges <- function(X) {
  if(!is(X)[1]=="GRanges") { stop("X wasn't GRanges") }
  prn <- nrow(X)
  XX <- vector("list",length(chrNames2(X)))
  do.one.chr <- function(X) {
    idz <- start(X); 
    idz <- idz[duplicated(idz)]
    udz <- unique(idz)
    n.id <- length(udz)
    #prv(udz)
    for (cc in 1:n.id) {
      select <- which(start(X) %in% udz[cc])
      U <- X[select,]
      ov <- findOverlaps(U)
      qq <- queryHits(ov)
      ss <- subjectHits(ov)  
      tt <- table(qq,ss)
      diag(tt) <- 0
      tt[lower.tri(tt)] <- 0
      to.kill <- colnames(tt)[colSums(tt)>0]
      if(length(to.kill)>0) {
        X <- X[-c(select[as.numeric(to.kill)]),]
      }
      # loop.tracker(cc,n.id); #cat("x")
    }
    #print(X)
    return(X)
  }
  ct <- 1
  all.chr <- chr2(X)
  chrn <- chrNames2(X)
  names(XX) <- chrn
  for(chrz in chrn) {
    XX[[chrz]] <- do.one.chr(chrSel(X,chrz))
    ct <- ct + 1
  }
  X <- do.call("rbind",args=XX)
  ND <- prn - nrow(X)
  if(ND>0) { cat("NOTE: removed",ND,"ids that were duplicates\n") }
  return(X)
}

# somehow is slower than reader?

read.snp.info <- function(input.fn=NULL, dir=getwd(), 
                            row.names=NULL, col.names=NULL,
                            ram.gb=2, verbose=TRUE, enforce.types=TRUE)
{
  if( all(dir=="") | is.null(dir) ) { dir <- getwd() }
  ## for compatibility with plumbCNV directory object
  ## Define data types for big.matrix
  # label
  input.fn <- cat.path(dir,input.fn)
  HIMEM <- 0.01 # GB
  # get the specs of the input file
  del <- get.delim(input.fn)
  inf <- quick.mat.format(input.fn)
  cls <- inf$ncol
  cN <- inf$cnames
  rws <- file.nrow(input.fn)
  mem <- estimate.memory(rws*cls)
  if(mem>HIMEM & verbose) { tracker <- TRUE } else { tracker <- FALSE }
  if(mem>ram.gb) { stop("file was larger than 'ram.gb' size specified (",ram.gb,"), aborting read") }
  row.mode <- "ALL"
  # make sure 'col.names' are valid
  if(is.null(col.names)) { col.names <- 1:inf$ncol }
  headerRow <- TRUE
  if(!inf$colnames) {
    if(is.character(col.names)) { 
      stop("input.fn must have column names as the first row ",
       "of the text file if 'col.names' are to be selected as character type (by names)") 
    } else {
      cN <- paste0("COL",1:inf$ncol) # insert dummy column names
      headerRow <- FALSE
    }
  } 
  if(inf$colnames & is.character(col.names)) { 
     col.names <- match(col.names,inf$cnames)
     if(all(is.na(col.names))) { stop("no specified column names were found") }
     if(any(is.na(col.names))) { warning("some specified column names were not found",
                                         ", continuing with the remainder") }
     col.names <- narm(col.names)
  }
  if(is.numeric(col.names)) {
    if(length(col.names)>cls) { 
      stop("more numeric column number ['col.name'] indexes were entered than the number of columns in 'input.fn'") }
    if(max(col.names)>cls) { 
      stop("some numeric column number ['col.name'] indexes were larger than the number of columns in 'input.fn'") }
    if(any(col.names!=round(col.names))) { stop("col.names (as column numbers) must be integers") }
    out.cols <- length(col.names) 
  }
  if(length(col.names)<1) { stop("no columns remain") }
  cnames.out <- cN[col.names] # output names are the chosen subset of column names
  ctypes.out <- inf$ctypes[col.names]
  
  # make sure 'row.names' are valid
  if(!inf$rownames) {
    if(!is.null(row.names)) { 
      stop("input.fn must have row names in the first column",
         " of the text file if 'row.names' are to be selected by names or row-numbers\n",
         "set 'row.names' to NULL, or correct the input file, to proceed") 
    } else {
      row.mode <- "NONE"
    }
  }
  dat.rows <- rws - (if(headerRow) {1} else {0})
  if(!is.null(row.names)) { 
    if(is.numeric(row.names)) { 
      if(length(row.names)>rws) { 
        stop("more numeric row number ['row.name'] indexes were entered than the number of rows in 'input.fn'") }
      last.row <- max(row.names)
      if(last.row>rws) { 
        stop("some numeric row number ['row.name'] indexes were larger than the number of rows in 'input.fn'") }
      if(any(row.names!=round(row.names))) { stop("row.names (as row numbers) must be integers") }
      out.rows <- length(row.names) 
      if(row.mode=="NONE") { row.mode <- "NINDEX" } else { row.mode <- "INDEX" }
    } else {
      if(!is.character(row.names)) {
        stop("row.names must be NULL, character, or numeric")
      }
      row.mode <- "NAMED"
    }
    out.rows <- length(row.names) # if above criteria have passed
  } else {
    out.rows <- dat.rows
  }

  ### MAIN IMPORT LOOP ###
  if(verbose) { cat("Creating data.frame object to store imported data\n") }
  impDat <- as.data.frame(matrix(character(),nrow=out.rows,ncol=out.cols), stringsAsFactors=FALSE)
  colnames(impDat) <- cnames.out
  if(row.mode=="NAMED") {  rownames(impDat) <- row.names  }
  rn.list <- character(length(dat.rows))
  if(row.mode %in% c("NINDEX","INDEX")) { irn.list <- character(length(out.rows)) }
  ext <- get.ext(input.fn)
  
  print(row.mode)
  
  if(ext %in% c("gz","zip","tar")) {
    dat.file <- gzfile(input.fn)
  } else {
    dat.file <- file(input.fn)
  }

  open(con=dat.file,open="r")

  ## read from matrix format txt file
  if(row.mode %in% c("NINDEX","INDEX","NAMED")) { cnt <- 0 } # set rows imported to 0
  if(headerRow) { ignore.me <- readLines(dat.file,n=1) } # ignore header row
  for (cc in 1:dat.rows) {
    #cat(cc,".")
    if(tracker) { loop.tracker(cc,dat.rows) }
    next.line <- readLines(dat.file,n=1)
    next.row <- strsplit(next.line,del,fixed=T)[[1]]
    if(row.mode=="NINDEX" | row.mode=="NONE") {
      next.row <- next.row[col.names] # select the column subset
    } else {
      rn.list[cc] <- next.row[1] # the 'rowname' (first column)
      next.row <- next.row[col.names+1]
    }
    #print(next.row); print(rn.list)
    # ALL all.rows - just import all, keep track of the row vec, add names
    # NAMED specific named rows - ignore if not in the named list, else place in its place
    # INDEX numbered row indexes - ignore if not in the numbered list
    # NINDEX numbered row indexes but no names - ignore if not in the numbered list, and ignore rownames
    # NONE none.rows - import all, but there are no rownames
    
    if(row.mode=="NONE" | row.mode=="ALL") {
      # just import everything
      impDat[cc,] <- next.row
    } else {
      if(row.mode=="NAMED") {
        # import a set of named rows
        ii <- match(rn.list[cc],row.names)
        if(!is.na(ii)) {
          impDat[ii,] <- next.row
          cnt <- cnt + 1 # successfully imported 1 more row
          if(cnt>=out.rows) { break } # all named rows imported, no need to continue
        } else {
          # next # row name was not in the import list
        }
      } else {
        #indexed, i.e, row.mode=="NINDEX" or "INDEX"
        if(cc %in% row.names) {
          cnt <- cnt + 1 # get the next row number
          impDat[cnt,] <- next.row
          irn.list[cnt] <- rn.list[cc]
        } else {
          # next # row number was not in the import list
          if(cc>=last.row) { break } # highest numbered row has passed, no need to continue
        }
      }
    }
  }
  close(dat.file)
  if(row.mode %in% c("NINDEX","INDEX","NAMED")) {
    if(verbose) { cat(out.of(cnt,length(row.names)),"'row.names' specified were imported\n") }
    if(!row.mode %in% "NAMED") {
      if(cnt <= dat.rows) { 
        rownames(impDat) <- irn.list 
      }  else {
        warning("more rows imported than were in file??")
      }
    }
  } else {
    # ie, all rows were imported NONE/ALL
    if(row.mode=="NONE") {
      rownames(impDat) <- paste(1:nrow(impData))
    } else {
      rownames(impDat) <- rn.list 
    }
  }
  # set appropriate column data types # 
  if(enforce.types) {
    if(verbose) { cat("converting columns to correct types ..")}
    CHAR_FRAC <- .2 # if number of different values is less than CHAR_FRAC% of nrow, make factor
    CHAR_N <- 20 # if number of different values is less than CHAR_N of nrow, make factor
    for(cc in 1:ncol(impDat)) {
      if(ctypes.out[cc]=="character") {
        if(length(unique(impDat[[cc]]))<20) {
          if((length(unique(impDat[[cc]]))/nrow(impDat))<.2) {
            impDat[[cc]] <- as.factor(impDat[[cc]])
          } else { next }
        } else {
          next # should already be character
        }
      } else {
        impDat[[cc]] <- as(impDat[[cc]],ctypes.out[cc])
      }
    }
    if(verbose) { cat(". done\n") }    
  }
  return(impDat)
}



if(F) {
  results <- reader("/chiswick/data/ncooper/imputation/T1D/sanity-1.RData")
  
  results <- results[,c(1:3,13)]
  results[["chr"]] <- 21
  
  results2 <- data.frame.to.GRanges(results,chr="chr",start="position",end="position")
  
  ichip.data <- reader("/chiswick/data/ncooper/imputation/COMMON/iChipPaperAnalysisResultsFeb2014.RData")
  
  ichip.data2 <- data.frame.to.GRanges(ichip.data[,1:5])
  
  regions <- reader("/chiswick/data/ncooper/imputation/COMMON/iChipFineMappingRegionsB37.RData")
  regions <- as(regions,"GRanges")
  res <- vector("list",188); st <- proc.time()[3]
  names(res) <- mcols(regions)[,"Band"]
  for (cc in 1:188) {
    res[[cc]] <- check.region(regions[180,],results2,ichip.data2)  
    loop.tracker(cc,188)
  }
  
  all.c <- sort(c(paste0("chr",3:22),"chr2-part-1","chr2-part-2"))
  locs <- cat.path("/chiswick/data/ncooper/imputation/COMBINED/single-snp-tests",all.c,ext="RData")
  results <- reader("/chiswick/data/ncooper/imputation/COMBINED/single-snp-tests/chr1.RData")
  results <- results[,c(1:4,15)]
  for (cc in 1:22) {
    X <- reader(locs[cc])
    X <- X[,c(1:4,15)]
    results <- rbind(results,X)
  }
  colnames(results)[5] <- "p"
  rownames(results) <- results$rs_id
  results$chr <- paste(results$chr)
  results$chr[results$chr=="2-part-1"] <- 2
  results$chr[results$chr=="2-part-2"] <- 2
  results2 <- data.frame.to.GRanges(results,chr="chr",start="position",end="position")
  nR <- nrow(regions)
  res <- vector("list",nR); 
  names(res) <- mcols(regions)[,"Band"]
  ichip.data3 <- select.autosomes(toGenomeOrder(ichip.data2))
  results2 <- toGenomeOrder2(results2)
  ######## Making plots for each Dense Region ########
  print(Dim(results2[rownames(all.snps)[!filt2],])); print(Dim(results2))
  results22 <- toGenomeOrder(results2[rownames(all.snps)[!filt3],])
  st <- proc.time()
  pdf("/chiswick/data/ncooper/imputation/COMBINED/T1Dtests/big188ComparisonB.pdf")
  for (cc in 1:nR) {
    #next.reg <- cat.path("/chiswick/data/ncooper/imputation/COMBINED/T1Dtests",
    #   "T1Dregion",suf=cc,ext="pdf")
    #pdf(next.reg)
#    par(mfrow=c(1,2))
    res[[cc]] <- check.region2(regions[cc,],results2,ichip.data3)
#    check.region2(regions[cc,],results22,ichip.data3)
    #dev.off()
    loop.tracker(cc,nR,st)
#    save(res,file="resultsOfRawRunT1D.RData")
  }
  dev.off()
  
  # check these results against impute stats #
  ref.dir <- "/chiswick/data/ncooper/imputation/COMBINED/impute_output/"
  each.chr <- list.files(ref.dir)
  X1 <- reader(cat.path(ref.dir,"imputed-2-part-1",ext="RData"))
  X2 <- reader(cat.path(ref.dir,"imputed-2-part-2",ext="RData"))
  X <- cbind2(X1$q,X2$q)
  ch2 <- snps(X)
  for (dd in c(1:22)[-2]) {
    XX <- reader(cat.path(ref.dir,"imputed-",suf=dd,ext="RData"))
    if(!is.null(XX$q)) {
      if(!is.null(XX$p)) {
        X <- rbind2(snps(XX$p),snps(XX$q))
      } else {
        X <- snps(XX$q)
      }
    } else {
      X <- snps(XX$p)
    }
    assign(paste0("ch",dd),X)
    loop.tracker(dd,22)
  }
  all.snps <- rbind(ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10,ch11,
            ch12,ch13,ch14,ch15,ch16,ch17,ch18,ch19,ch20,ch21,ch22)
  save(all.snps,file="/chiswick/data/ncooper/imputation/COMBINED/allSnpsInfoFromChris2.RData")
  all.dif <- unlist(sapply(res,"[",1))
  all.outly <- unlist(sapply(res,"[",2))
  ## density plot overlaying overall distribution excluding these two lists, then dists for the lists
  ## do for each stat
  main.set <- which(!rownames(all.snps) %in% c(all.outly,all.dif))
  dif.set <- which(rownames(all.snps) %in% all.dif)
  out.set <- which(rownames(all.snps) %in% all.outly)
  statz <- c("exp_freq_a1","info","certainty","type",
              "info_type0","concord_type0","r2_type0") #,"A1","A2")
  pdf("mytestplot.pdf")
  filt <- rep(TRUE,nrow(all.snps))
 # filt2 <- with(all.snps,type>=0 & info>.825 & exp_freq_a1>.001 & exp_freq_a1<.006)
 # filt3 <- with(all.snps,type>=0 & info>.825  & exp_freq_a1<.006)
  for (dd in 1:length(statz)) {
    X <- all.snps[[statz[dd]]]
    plot(density(X[main.set][filt[main.set]]),col="grey",main=statz[dd])
    lines(density(X[dif.set][filt[dif.set]]),col="red")
    lines(density(X[out.set][filt[out.set]]),col="orange")
    legend("topright",
        legend=c("large diff to ichip","high p value","all others"),
        col=c("red","orange","grey"),lwd=1)
  }
  dev.off()
  categ2 <- rep(0,nrow(all.snps))
  categ2[out.set] <- 1
  all.snps[["category"]] <- categ2
  quick.lda(Y="category",preds=c("exp_freq_a1","info","certainty"),
       labs=c("main","difs"),data=all.snps,prior=c(0.65,0.35))
 
  library(party)
  ctree("category ~ .",data=all.snps)
  pdf("ctree.test.pdf")
  ctree(as.formula(paste("categ2 ~ ",paste(colnames(all.snps)[-c(1:3,11:13)],collapse=" + "))),data=all.snps)
  dev.off()
 
  categ <- rep(0,nrow(all.snps))
  categ[out.set] <- 1
  all.snps[["category"]] <- categ
  quick.lda(Y="category",preds=c("exp_freq_a1","info","certainty"),
       labs=c("main","outs"),data=all.snps,prior=c(0.65,0.35))
  
 # Node 54,55, 56 type>0, info>.825, info_type0<.538, exp_freq_a1>.001 & <.004, especially concord_type0<.994
 # Node 60 type>0, info>.825, info_type0>.212 & <.538, exp_freq_a1>.001 & <.006, r2_type0<.2
 # Node 70 type>0, info>.825, info_type0>.538, exp_freq_a1>.001 & <.003, r2_type0>.914
 # Node 72 type>0, info>.825, info_type0>.538, exp_freq_a1>.001 & <.004, info_type0 < .678
  
  #filt1 <- with(all.snps,type>0 & info>.825 & exp_freq_a1>.001 & exp_freq_a1<.006)
  
}

#16700000, 16840000
#40450000, 45700000

#quick.lda(Y="categ",preds=c("exp_freq_a1","info","certainty"),labs=c("main","difs"))

# runs LDA and automatically calculates sensitivity, specificity, etc
quick.lda <- function(Y,preds,data,labs=c("grp1","grp2"),prior=c(0.2,0.8)) {
  require(MASS)
  Y_var <- Y
  next.form <- as.formula(paste(Y_var,"~",paste(preds,collapse="+")))
  next.mod <- lda(next.form,data=data ,prior=prior)
  uu <- predict(next.mod)$class; vv <- labs[(1+data[[Y_var]])]
  resulto <- table(uu,vv)
  X <- predict(next.mod)$x
  Y <- predict(next.mod)$class
  senso <- resulto[2,2]/sum(resulto[,2])
  speco <- resulto[1,1]/sum(resulto[,1])
  ppp <- resulto[2,2]/sum(resulto[2,])
  npp <- resulto[1,1]/sum(resulto[1,])
  acco <- sum(diag(resulto))/(sum(resulto))
  rownames(resulto) <- c(paste("Pred",labs,sep="-")) ; print(resulto); 
  cat("\nSensitivity:",round(senso,3)," Specificity",round(speco,3),"\n")
  cat("PPP:",round(ppp,3),"         NPP",round(npp,3),"\n")
  cat("Overall Accuracy:",round(acco,3),"\n")
  return(next.mod)
}


# check 1 T1D region between ichip results and new results
check.region <- function(region,result,chip) {
  region <- region[1,] # only check 1 region
  chip  <- subsetByOverlaps(chip,region)
  resultA <- result
  ## overlapping SNP correlations ##
  equal <- subsetByOverlaps(result,chip)
  equal.pos <- start(equal);   equal.p <- mcols(equal)[,"p"]
  chip.pos <- start(chip);   chip.p <- mcols(chip)[,"p.value"]
  chip.ps <- chip.p[match(equal.pos,chip.pos)]
  main.cor <- cor(-log10(chip.ps),-log10(equal.p),use="pairwise.complete")
  ## NON-overlapping SNP correlations ##
  result <- resultA[countOverlaps(resultA, chip) <= 1L] # get only non overlapping snps
  en <- nrow(result)
  test.regs <- makeGRanges(chr=chr(result)[-1],start=start(result)[-en],end=start(result)[-1])
  p.set <- vector("list",nrow(test.regs))
  for (cc in 1:nrow(test.regs)) {
    p.set[[cc]] <- mcols(subsetByOverlaps(chip,test.regs[cc,]))[["p.value"]]
  }
  s.set <- vector("list",nrow(result))
  s.set[[1]] <- p.set[[1]]; s.set[[length(s.set)]] <- tail(p.set,1)[[1]]
  for(cc in 2:length(p.set)) {
    s.set[[cc]] <- -log10(c(p.set[[cc-1]],p.set[[cc]]))
  }
  #return(s.set)
  set.mins <- sapply(s.set,minna)
  set.means <- sapply(s.set,meanna)
  set.medians <- sapply(s.set,medianna)
  rez <- -log10(mcols(result)[,"p"])
  badz <- is.na(rez) | is.infinite(rez) | sapply(set.mins,length)==0
  min.cor <- cor(rez[!badz],set.mins[!badz],use="pairwise.complete")
  mean.cor <- cor(rez[!badz],set.means[!badz],use="pairwise.complete")
  median.cor <- cor(rez[!badz],set.medians[!badz],use="pairwise.complete")
  return(list(main=main.cor,min=min.cor,mean=mean.cor,median=median.cor))
}

# check.region(regions[180,],results2,ichip.data2)


# check 1 T1D region between ichip results and new results
# region: dense regions
# result: snp results
# chip: all ichip snps
#'  pdf("testmyplot.pdf"); check.region2(regions[180,],results2,ichip.data2); dev.off()
check.region2 <- function(region,result,chip,sd.thr=5) {
  region <- region[1,] # only check 1 region
  chip  <- subsetByOverlaps(chip,region) # all snps in region
  result  <- subsetByOverlaps(result,region)
  resultA <- result
  ## overlapping SNP correlations ##
  equal <- subsetByOverlaps(result,chip) # snp results for ichip snps
  equal.pos <- start(equal);   equal.p <- mcols(equal)[,"p"]
  chip.pos <- start(chip);   chip.p <- mcols(chip)[,"p.value"]; chipL10 <- -log10(chip.p)
  chip.ids <- rownames(chip)
  chip.ps <- chip.p[match(equal.pos,chip.pos)]
  main.dif <- (-log10(equal.p)) - (-log10(chip.ps))
  #print(summary(main.dif))
  ## NON-overlapping SNP correlations ##
  result <- resultA[countOverlaps(resultA, chip) <= 1L] # get only non overlapping snps
  # find any with large difs
  ind <- which.outlier(main.dif,low=FALSE)
  outly.main.snps <- rownames(equal)[ind]
  # select set with small difs as a reliable set
  ind2 <- which.outlier(main.dif,low=FALSE,thr=1)
  tight.main.snps <- rownames(equal)[-ind2]
  # within the region, are any particularly large? 3SD or some thresh?
  tight.ps <- -log10(mcols(equal[tight.main.snps,])[,"p"])
  mn <- meanna(tight.ps); sD <- sdna(tight.ps)
  thr <- mn+(sd.thr*sD)
  # plot region 3 colz, ichip and ichip-impute, and new-impute,
  # color anything close to being an outlier, use identify()
  mcols(resultA)[["log10"]] <- -log10(mcols(resultA)[["p"]])
  L10s <- -log10(mcols(resultA)[["p"]])
  if(length(L10s)>1) {
    yl <- extendrange(r=range(L10s,na.rm=T),f=.2); if(any(is.na(yl))) { print(L10s) ; yl <- NULL }
    plotRanges(resultA,alt.y="log10",do.labs=F,col="grey",pch=".",ylim=yl, 
                main=paste("Chr",mcols(region)[,"Band"]),ylab="-log10 p-values",xlab="position")
  } else { warning(paste("no data for region")); return(list(dif=NULL,outly=NULL)) }
  if(length(chip.pos)>0 & length(chipL10)==length(chip.pos)) {
    points(chip.pos,chipL10,col="blue",pch=".")
  }
  big.ones <- which(L10s>thr)
  if(length(big.ones)>0) {
    plotRanges(resultA[big.ones,],alt.y="log10",do.labs=F,col="orange",skip.plot.new=TRUE) 
  }
  big.ones2 <- which(L10s>(thr*2))
  if(length(big.ones2)>0) {
    plotRanges(resultA[big.ones2,],alt.y="log10",do.labs=T,col="orange1",skip.plot.new=TRUE) 
  }
  more.big.ones <- which(chipL10>thr)
  if(length(more.big.ones)>0) {
    points(chip.pos[more.big.ones],chipL10[more.big.ones],col="green") 
    text(chip.pos[more.big.ones],chipL10[more.big.ones],labels=chip.ids[more.big.ones],cex=.7,col="darkgreen",pos=4) 
  }
  if(length(outly.main.snps)>0) {
    plotRanges(resultA[outly.main.snps,],alt.y="log10",do.labs=T,col="red",skip.plot.new=TRUE) 
    points(chip.pos[chip.ids %in% outly.main.snps],chipL10[chip.ids %in% outly.main.snps],col="purple") 
    for(dd in 1:length(outly.main.snps)) {
      one <- resultA[outly.main.snps[dd],]
      if(any(chip.pos==start(one))) {
        YYy <- c(-log10(mcols(one)[,"p"]),chipL10[chip.pos==start(one)][1])
        if(length(YYy)==2 & all(!is.na(YYy))) {
          lines(rep(start(one),2),YYy,col="red")
        }
      }
    }
  }
  legend("bottom",legend=c(paste0("imputed > ",sd.thr,"SD of reliable"),"> 3SD vs ichip",
         paste0("ichip > ",sd.thr,"SD of reliable"),"ichip pvalues","imputation pvalues"),
         col=c("orange","red","green","blue","grey"),pch=1,bty="n",ncol=3, cex=.75)
  return(list(dif=outly.main.snps,outly=rownames(resultA)[big.ones],ichip=chip.ids[more.big.ones])) 
}


# ichip.sets <- c("T1D","COELIAC","GRAVES","JIA","MS","RA1","RA2","RA3","RA4","RA5","RA6","IC.CASE","IC.CTRL")
# gwas.sets <- c("Cases_UKC","Cases_UKN","Cases_UKP","Cases_UKW","Controls_Illu58C","Controls_IlluNBS")

# personal function to return a SnpMatrixList for any of the project disease groups/sub-groups #
get.sml <- function(disease="T1D",sub="",ms.sub="",dir.form=FALSE) {
  diseasez <- c("T1D","COELIAC","GRAVES","JIA","MS","RA")
  subz <- c("IC.CASE","IC.CTRL","WTCCC2",paste0("RA",1:6))
  sub.d <- c(rep("MS",3),rep("RA",6))
  ms.subz <- c("Cases_UKC","Cases_UKN","Cases_UKP","Cases_UKW","Controls_Illu58C","Controls_IlluNBS")
  dis <- disease[1]
  disease <- toupper(dis); sub <- toupper(sub[1]); ms.sub <- ms.sub[1]
  if(is.na(sub)) { sub <- "" }
  if(is.na(ms.sub)) { ms.sub <- "" }
  if(is.null(sub)) { sub <- "" }
  if(is.null(ms.sub)) { ms.sub <- "" }
  if(!disease %in% diseasez) { 
    if(!disease %in% subz) {
      if(!disease %in% toupper(ms.subz)) {
        stop("disease must be one of: ",comma(diseasez)) 
      } else {
        ms.sub <- dis; sub <- "WTCCC2"; disease <- "MS"
      }
    } else {
      sub <- disease
      disease <- sub.d[match(sub,subz)]
    }
  }
  if(disease %in% sub.d) {
    ii <- which(sub.d %in% disease)
    if(!sub %in% subz[ii]) { stop("disease was ",disease," which needs one of the following values for sub:",comma(subz[ii]),"\n") }
  } else { sub <- "" }
  if(disease!="MS") { ms.sub <- "" }
  if(!sub %in% c(subz,"")) { stop("sub must be one of: ",comma(subz)) } else {
    if(sub=="WTCCC2") {
      if(!ms.sub %in% c(ms.subz)) { stop("ms.sub is required when sub='WTCCC2'' and must be one of: ",comma(ms.subz)) }
    }
  }
  dd <- paste0("/chiswick/data/ncooper/imputation/",disease,"/PQDATA/")
  if(sub!="") {  
    dd <- paste0(dd,sub,"/")
    if(ms.sub!="") {  dd <- paste0(dd,ms.sub,"/") }
  }
  sml <- list.files(dd)
  must.use.package("gtools")
  sml <- gtools::mixedsort(sml)
  if(!dir.form) {
    sml <- cat.path(dd,sml)
  }
  sml <- as.list(sml)
  #sml <- sml[c(18:19,26:40,1:17,20:24)]
  if(dir.form) {
    return(list(sml=sml,dir=dd))
  } else {
    return(sml)
  }
}


#
# PCs <- explore.PCs.13("T1D",cutoff=10^-7,PC.out=T,file="PCAT1D"); save(PCs,file="PCsT1D.RData")
# PCs <- explore.PCs.13("GRAVES",cutoff=10^-4.3,PC.out=T,file="PCAGRAVES"); save(PCs,file="PCsGRAVES.RData")
# PCs <- explore.PCs.13("COELIAC",cutoff=10^-6,PC.out=T,file="PCACOELIAC"); save(PCs,file="PCsCOELIAC.RData")
# PCs <- explore.PCs.13("JIA",cutoff=10^-6,PC.out=T,file="PCAJIA"); save(PCs,file="PCsJIA.RData")
# PCs <- explore.PCs.13("RA1",cutoff=10^-6,PC.out=T,file="PCARA1"); save(PCs,file="PCsRA1.RData")
# PCs <- explore.PCs.13("RA2",cutoff=10^-5.5,PC.out=T,file="PCARA2"); save(PCs,file="PCsRA2.RData")
# PCs <- explore.PCs.13("RA3",cutoff=10^-5,PC.out=T,file="PCARA3"); save(PCs,file="PCsRA3.RData")
# PCs <- explore.PCs.13("RA4",cutoff=10^-4.5,PC.out=T,file="PCARA4"); save(PCs,file="PCsRA4.RData")
# PCs <- explore.PCs.13("RA5",cutoff=10^-6.5,PC.out=T,file="PCARA5"); save(PCs,file="PCsRA5.RData")
# PCs <- explore.PCs.13("RA6",cutoff=10^-6,PC.out=T,file="PCARA6"); save(PCs,file="PCsRA6.RData")
# PCs <- explore.PCs.13("IC.CTRL",cutoff=10^-6,PC.out=T,file="PCAIC.CTRL"); save(PCs,file="PCsIC.CTRL.RData")
# PCs <- explore.PCs.13("IC.CASE",cutoff=10^-6,PC.out=T,file="PCAIC.CASE"); save(PCs,file="PCsIC.CASE.RData")

explore.PCs.13 <- function(disease,cutoff=10^-5,plot=FALSE, PC.out=FALSE,
                           file="",pch.out=19, pch.clean=1, out34in12=FALSE) {
  pca.pred.a.tg <- reader("/chiswick/data/ncooper/imputation/THOUSAND/PCs1000g.RData")
  sample.info <- reader("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData")
  anc <- factor(sample.info$ancestry[match(rownames(pca.pred.a.tg),rownames(sample.info))])
  (load(paste0("PCs",disease,".RData")))
  ii <- Moutlier(PCs[,1:4],quantile=1-cutoff,plot=plot)
  ii2 <- Moutlier(PCs[,1:2],quantile=1-(cutoff*100),plot=plot)
  ww <- (which(ii$md>ii$cutoff))
  ww2 <- (which(ii2$md>ii2$cutoff))
  if(nchar(paste(file))>1) { pdf(cat.path(getwd(),file,ext="pdf"),width=15) }
  par(mfrow=c(1,2))
  ### PC 1 vs 2 ###
  plot(pca.pred.a.tg[,1],pca.pred.a.tg[,2],col=get.distinct.cols(14)[as.numeric(anc)],
       main="PCs 1-2",xlab="PC-1",ylab="PC-2")
  points(PCs[,1],PCs[,2],pch=pch.clean)
  points(pca.pred.a.tg[,1],pca.pred.a.tg[,2],col=get.distinct.cols(14)[as.numeric(anc)])
  if(out34in12) { 
    #points(PCs[ww,1],PCs[ww,2],col="white",pch=pch.clean) 
    points(PCs[ww,1],PCs[ww,2],col="orange",pch=pch.out) 
  }
  #points(PCs[ww2,1],PCs[ww2,2],col="white",pch=pch.clean)
  points(PCs[ww2,1],PCs[ww2,2],col="red",pch=pch.out)
  legend("topleft",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=3,bty="n")
  
  ### PC 3 vs 4 ###
  plot(pca.pred.a.tg[,3],pca.pred.a.tg[,4],col=get.distinct.cols(14)[as.numeric(anc)],
       main="PCs 3-4",xlab="PC-3",ylab="PC-4")
  points(PCs[,3],PCs[,4],pch=pch.clean) 
  points(pca.pred.a.tg[,3],pca.pred.a.tg[,4],col=get.distinct.cols(14)[as.numeric(anc)])
  if(out34in12) {  
    #points(PCs[ww2,3],PCs[ww2,4],col="white",pch=pch.clean) 
    points(PCs[ww2,3],PCs[ww2,4],col="red",pch=pch.out) 
  }
  #points(PCs[ww,3],PCs[ww,4],col="white",pch=pch.clean)
  points(PCs[ww,3],PCs[ww,4],col="orange",pch=pch.out)
  legend("topleft",legend=c("outliers PC 1 vs 2","outliers PC 3 vs 4",
                            paste("Non-outlier",disease),"thousand genomes (all other colours, see adjacent plot)"),
         col=c("red","orange","black","grey"),pch=c(pch.out,pch.out,pch.clean,1),bty="n")
  mtext(paste(disease,"ancestry PCA"), side = 3, line = -2, outer = TRUE,cex=1.5)
  #####
  if(nchar(file)>0) { dev.off() }
  outliers <- rownames(PCs)[ww]
  outliers2 <- rownames(PCs)[ww2]
  outliers <- outliers[!outliers %in% outliers2]
  if(PC.out) {
    PCs[["outlier.code"]] <- 0
    PCs[["outlier.code"]][rownames(PCs) %in% outliers] <- 1
    PCs[["outlier.code"]][rownames(PCs) %in% outliers2] <- 2
    return(PCs)
  } else {
    return(list(definite=outliers2,borderline=outliers))
  }
}

