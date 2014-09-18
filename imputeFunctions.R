# internal function
nam <- function(x) { y <- narm(x); return(y[y!="NA"]) }


##### contains functions #####
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
# annot.sep.support - Create aSnpMatrix from SnpMatrix plus SNP/sample support files
# split.pq - split an aSnpMatrix into a long and short arm (arm.p and arm.q)
##############################



#internal
comma <- function(...) {
  paste(...,collapse=",")
}

# used for when chip has duplicate positions for which you only want to keep 1 of each, and you 
# have a continuous criteria to maximise, e.g, call rate that you can use to decide between them
choose.dups.to.drop <- function(snp.info, snp.excl=NULL, ind.fn="../iChipDupPairs.RData", col.name="call.rate") {
  dup.mat <- reader(ind.fn)
  if(!is.null(snp.excl)) {
    dup.mat <- apply(dup.mat,1,function(x) { x[which(x %in% snp.excl)] <- NA; return(x) })
  }
  dup.mat <- narm(t(dup.mat))

  cal.mat <- matrix(numeric(),nrow=nrow(dup.mat),ncol=ncol(dup.mat))
  if(!any(colnames(snp.info) %in% col.name)){ stop("column name ",col.name," was not found in snp.info") }
  ind1 <- match(dup.mat[,1],rownames(snp.info))
  ind2 <- match(dup.mat[,2],rownames(snp.info))
  notfound <- c((dup.mat[,1][which(is.na(ind1))]),(dup.mat[,2][which(is.na(ind2))]))
  if(any(is.na(ind1)) | any(is.na(ind2))) { stop("some snps: ",comma(notfound)," in the duplicate matrix were not in the snp.info, add 'snp.excl' entries for such") }
  cal.mat[,1] <- snp.info[ind1,][[col.name]]
  cal.mat[,2] <- snp.info[ind2,][[col.name]]
  j <- apply(cal.mat,1,function(x) { head(which(x==min(x,na.rm=T)),1) })
  i <- 1:nrow(cal.mat)
  DUPPOS.EXCL <- dup.mat[cbind(i,j)]
  return(DUPPOS.EXCL)
}


# you want to include certain SNPs missing from a datafile as columns explicitly set to all missing (00)
add.missing.snps <- function(X, snpSup) {
  all.snps <- rownames(snpSup)
  to.add <- all.snps[!all.snps %in% colnames(X)]
  add.m <- matrix(raw(),nrow=nrow(X),ncol=length(to.add))
  rownames(add.m) <- rownames(X)
  colnames(add.m) <- to.add
  add.sm <- as(add.m,"SnpMatrix")
  add.asm <- new("aSnpMatrix",.Data=add.sm, snps=snpSup[to.add,],
    samples=X@samples,phenotype=X@phenotype,alleles=X@alleles)
  m2 <- cbind2(X,add.asm)
  m2 <- m2[,all.snps] #reorder
  return(m2)
}


# a set of SnpMatrix objects with mostly the same SNPs, but different omitted, this helps
# synchronise them by inserting missing vectors of snps not missing in at least one dataset
# into the datasets where they are not present
sync.asnp.mats <- function(myMat) {
  nn <- length(myMat)
  sMat <- vector("list",nn)
  for (cc in 1:nn) {
    sMat[[cc]] <- myMat[[cc]]@snps
  }
  snpso <- do.call(rbind,args=sMat)
  snpso <- snpso[!duplicated(snpso$snp.name),]
  for (cc in 1:nn) {
    myMat[[cc]] <- add.missing.snps(myMat[[cc]],snpso)
  }
  return(myMat)
}


## for an aSnpMatrix, extract all snps from a given chromosome
get.chr <- function(x,chr) {
  select <- which(x@snps$chromosome==chr)
  if(length(select)>0) {
    return(x[,select])
  } else {
    return(NULL)
  }
}

## from a list of SnpMatrix objects with all SNPs but different samples,
# extract all samples for a given chromosome
Get.chr <- function(X,CHR) {
  print(CHR) # abysmal memory usage - big fail!
  each.chr <- lapply(X,get.chr,chr=CHR)
  comb.chr <- do.call("rbind3",args=each.chr)
  return(comb.chr)
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
  
  
#' Flip the strand for allele codes
#' 
#' When strands are being converted sometimes strand flipping is necessary,
#' whereby C goes to G, A goes to T, etc. This function takes care of this
#' in the most intuitive way possible, and can also handle IUPAC ambiguity 
#' codes, and separators within the text.
#' @param acgt character, allele codes, a series of letters, may include
#' IUPAC ambiguity codes, "R","Y","M","K","S","W","H","B","V","D","N", in
#' addition to the standard A, C, G, T codes. Can also involve a separator
#' in the text, e.g, C/T, A/T, etc, as long as the seperator being used
#' is passed to the 'sep' parameter.
#' @param sep optional separator (ignore if entering single letters), that
#' goes between pairs of alleles. The function will work if more that 2
#' alleles are used as long as multiple separators are used
#' @export
#' @return returns a result in the same character format as that entered,
#' including the same separator if any, but with the strands flipped
#' @examples
#' flip.strand(c("A","C","G","M"))
#' flip.strand(c("A/C","C/T","G/A","M/B"))
#' flip.strand(c("A|C|T|G","C|T","G|A|C","M|B"),sep="|")
flip.strand <- function(acgt,sep="/") {
  new <- paste(acgt)
  nn <- length(acgt)
  sepz <- grep(sep,acgt,fixed=TRUE)
  if(length(sepz)> 0.5*nn) {
    two.sets <- strsplit(acgt,split=sep,fixed=TRUE)
    lenz <- sapply(two.sets, length)
    #if(any(lenz>2)) {
    # warning("A maximum of 2 allele codes should be entered with a separator, e.g, A/C",
    #          ". If more types are possible, use IUPAC ambiguity codes") }
    result <- sapply(lapply(two.sets, flip.strand),paste,sep="",collapse=sep)
    return(result)
  } else {
    if(length(which(nchar(acgt)>1)) > 0.5*nn) { 
      warning("allele codes should only be 1 character long, the majority of 'acgt' entries were longer than this")
    } else {
      ## valid, and continue
    }
  }
  if(is(new)[1]!="character") { stop("invalid input, needs to be coercible to characters") }
  # define recoding characters, including uncertain codes #
  pre <- c("A","T","G","C","R","Y","M","K","S","W","H","B","V","D","N")
  post <- c("T","A","C","G","Y","R","K","M","S","W","D","V","B","H","N")
  pre <- c(toupper(pre),tolower(pre)) # create upper and lower case versions
  post <- c(toupper(post),tolower(post)) # create upper and lower case versions
  if(!any(acgt %in% pre)) { warning("no characters in 'acgt' were valid allele codes") }
  if(length(pre)!=length(post)) { stop("internal function error, lists invalid, please report bug") }
  for (cc in 1:length(pre)) {
    new[acgt==pre[cc]] <- post[cc]
  }
  return(new)
}



#' Obtain an index of all members of values with duplicates (ordered)
#' 
#' The standard 'duplicated' function, called with which(duplicated(x)) will 
#' only return the indexes of the extra values, not the first instances. For instance
#' in the sequence: A,B,A,C,D,B,E; it would return: 3,6. This function will also
#' return the first instances, so in this example would give: 1,3,2,6 [note it
#' will also be ordered]. This index can be helpful for diagnosis if duplicates 
#' are unexpected, for instance in a data.frame, and you wish to compare the differences
#' between the rows with the duplicate values occuring.
#' @param x a vector that you wish to extract duplicates from
#' @param pairs.only logical, whether to assume that there are exactly 2 copies
#' for each duplicate - this will be slightly faster, and result in a slightly
#' easier to parse output as you know every second value will be the index.
#' If this is used when there are more than 2 copies of some 'x', additional
#' copies will be left out of the result. If pairs.only=FALSE, then sets 
#' of any length can be returned.
#' @return vector of indices of which values in 'x' are duplicates (including
#' the first observed value in pairs, or sets of >2), ordered by set, then
#' by appearance in x. If pairs.only=FALSE, then sets can have length >=2,
#' or if pairs.only = TRUE, then sets will all have length =2, and any
#' in-between duplicates will be left out of the listing. Only use
#' @examples
#' dup.pairs(c(1,1,2,2,3,4,5,6,2,2,2,2,12,1,3,3,1))
dup.pairs <- function(x,pairs.only=FALSE) {
  if(pairs.only) {
    # fast if you know the max is pairs #
    df <- cbind(which(rev(duplicated(rev(x)))),which(duplicated(x)))
    return(as.vector(t(df)))
  } else {
    dx <- duplicated(x)
    other.dups <- which(dx)
    not.and.first <- which(!dx)
    ind.dups <- not.and.first[which(x[not.and.first] %in% x[other.dups])]
    xo <- x[other.dups]
    vc <- vector()
    for (cc in 1:length(ind.dups)) {
      vc <- c(vc,ind.dups[cc],other.dups[which(xo %in% x[ind.dups[cc]])])
    }
    return(vc)
  }
}

#internal
#'
#' @examples
#' ichip.ids <- sync.id.by.pos(rr,sup.for.impute,FALSE)
#' ref.ids <- sync.id.by.pos(rr,sup.for.impute,TRUE)
sync.id.by.pos <- function(ref.snp,snp.inf,return.ref.ids=FALSE) {
  if(!"position" %in% colnames(ref.snp)) { stop("ref.snp must have column 'position'") }
  if(!is(snp.inf)[1] %in% "ChipInfo") { stop("snp.inf must be a ChipInfo object") }
  posz <- start(snp.inf)
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


# split an aSnpMatrix into a long and short arm (arm.p and arm.q)
split.pq <- function(aSnpMat,build=37,pqvec=FALSE,verbose=TRUE) {
  si <- snp.from.annot(aSnpMat)
  if(!is(si)[1]=="RangedData") { si <- as(si,"RangedData") }
  if(!is(si)[1]=="RangedData") { stop("couldn't extract snp.info from aSnpMat") }
  cent <- get.centromere.locs(build=build)
  cnt.ch <- chr2(cent)
  chrz <- chrNames2(si)
  chrz2 <- gsub("chr","",paste(chrz),ignore.case = TRUE)
  arm <- rep("p",nrow(si))
  for(cc in 1:length(chrz)) {
    ii <- which(cnt.ch==chrz2[cc])
    if(length(ii)>0) { 
      if(length(ii)>1) {
        warning("more than one centromere entry for a single chromosome - import failure likely")
      }
      nxt.chr <- chr.sel(si,paste(chrz[cc]))
      cent.loc <- start(cent[ii[1],])
      if(verbose) { cat("centromere location:",cent.loc,"\n") }
      posz <- start(nxt.chr)
      indx <- match(rownames(nxt.chr),rownames(si))
      #prv(cent.loc,posz)
      arm[indx][posz > cent.loc] <- "q"
    } else {
      if(!chrz2[cc] %in% c("XY","MT","M")) {
        warning("centromere not found for chromosome labelled:",chrz2[cc])
      }
    }
  }
  #bb <- Band(si,build=37)
  #arm <- rep("p",nrow(si))
  #arm[grep("q",bb)] <- "q"
  if(pqvec) { return(arm) }
  if(length(unique(arm))==1) { warning("only arm values of ",unique(arm)," were found"); return(aSnpMat) }
  mat1 <- aSnpMat[,arm=="p"]
  mat2 <- aSnpMat[,arm!="p"]
  if(verbose) { cat("split aSnpMatrix into",ncol(mat1),"on arm 'p' and",ncol(mat2),"on arm 'q'\n")  }
  out <- list(arm.p=mat1,arm.q=mat2)
  return(out)
}

#' Create aSnpMatrix from SnpMatrix plus SNP/sample support files
#' 
#' @param snpMat a SnpMatrix or XSnpMatrix
#' @param snp.info annotation with rownames as snp labels, then chromosome, position and allele 
#' codes, can be RangedData, GRanges, ChipInfo or data.frame type
#' @param sample.info data.frame with rownames of sample ids, may also contain 'sex' and some
#' phenotype variable, such as 'pheno' or 'affected' or 'case' (autodetected), may also contain
#' pedigree information, such as a plink .fam/ped file
#' @param snp.excl text filename or character vector containing SNP ids to exclude
#' @param samp.excl text filename or character vector containing sample ids to exclude
#' @export
#' @return returns an annotSnpStats::aSnpMatrix object with @samples and @snps info taken
#' from the inputted support files. Will only include samples and SNPs occuring in both
#' the dataset and the support files.
#' @examples
#' ## my own example for now, only works @DIL ##
#' t1.dr <- "/chiswick/data/ncooper/imputation/T1D/"
#' setwd(t1.dr)
#' smp <- reader("/chiswick/data/ncooper/iChipData/sample.info.RData")
#' snp.excl <- "/chiswick/data/ncooper/imputation/T1D/snpsExcluded.txt"
#' smp.excl <- "/chiswick/data/ncooper/imputation/T1D/sampsExcluded.txt"
#' print(load("CHR20.RData"))
#' x <- gt.for.impute
#' si <- sup.for.impute
#' aSnpMat <- annot.sep.support(x,si,smp, snp.excl=snp.excl, samp.excl=smp.excl)
annot.sep.support <- function(snpMat, snp.info, sample.info, snp.excl=NULL, 
                              samp.excl=NULL, warn.int=TRUE, genome.order=TRUE, verbose=TRUE) {
  #### check for valid input types ####
  if(!is(snpMat)[1] %in% c("XSnpMatrix","SnpMatrix")) { 
    if(snp.mat.list.type(snpMat)!="error") {
      snpMat <- snpMatLst.collate(snpMat)
    } else {
      stop("snpMat must be of type 'SnpMatrix' or 'XSnpMatrix' or a 'snpMatrixList'") 
    }
  }
  if(is(snp.info)[1] %in% c("RangedData","GRanges","ChipInfo")) { 
    snp.info <- as(snp.info,"GRanges") 
    ii <- grep("elementMetadata.",colnames(mcols(snp.info)))
    if(length(ii)>0) {
      colnames(mcols(snp.info)) <- gsub("elementMetadata.","",colnames(mcols(snp.info)))
    }
  } else { 
    snp.info <- force.frame(snp.info)
    snp.info <- column.salvage(snp.info,"pos",c("Start","Pos","pos","POS","Pos37","pos37","Pos36","pos36"),ignore.case=FALSE)
    snp.info <- column.salvage(snp.info,"chr",c("chr","chromosome","chrom","chromo","space","seqnames"),ignore.case=TRUE)
    #snp.info <- column.salvage(snp.info,"end",c("end"),ignore.case=TRUE)
    print(head(snp.info))
    snp.info <- data.frame.to.ranged(snp.info,GRanges=TRUE)
  }
  if(!is(snp.info)[1] %in% c("GRanges")) { stop("snp.info must be coercible to a data.frame") }
  if(genome.order) { snp.info <- toGenomeOrder(snp.info) }
  sample.info <- force.frame(sample.info)
  if(!is(sample.info)[1] %in% c("data.frame")) { stop("sample.info must be coercible to a data.frame") }
  ###############################
  # find the right columns in the sample file
  sample.info <- column.salvage(sample.info,"affected",c("affected","phenotype","pheno","case","cases","phenot","phen","aff","ph","controls","ctrls","ctrl","control"),ignore.case=TRUE)
  # extract optional fields
  suppressWarnings({
    sample.info <- column.salvage(sample.info,"sex",c("sex","gender","mf"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"father",c("father","paternal","pat"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"mother",c("mother","maternal","mat"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"pedigree",c("pedigree","family","fam"),ignore.case=TRUE);
    sample.info <- column.salvage(sample.info,"id",c("id","sampleid","caseid","member"),ignore.case=TRUE)
    sample.info <- column.salvage(sample.info,"plate",c("plate","array.plate","plateid","plate.id"),ignore.case=TRUE)
  })
  if(!is.null(snp.excl)) { snp.excl <- force.vec(snp.excl) }
  if(!is.null(samp.excl)) { samp.excl <- force.vec(samp.excl) }
  ## now can assume valid input parameters have been checked ##
  samps.in.dat <- rownames(snpMat)
  snps.in.dat <- colnames(snpMat)
  samps.in.support <- rownames(sample.info)
  snps.in.support <- rownames(snp.info)
  if(warn.int) {
    # warn when sample/snp support doesn't perfectly match the SnpMatrix row/colnames #
    if(!all(samps.in.dat %in% samps.in.support)) { warning("sample list in support differs from data, will take intersection")}
    if(!all(snps.in.dat %in% snps.in.support)) { warning("SNP list in support differs from data, will take intersection")}
  }
  final.samps <- samps.in.dat[(samps.in.dat %in% samps.in.support) & (!samps.in.dat %in% samp.excl)]
  final.snps <- snps.in.dat[(snps.in.dat %in% snps.in.support) & (!snps.in.dat %in% snp.excl)]
  if(length(final.snps)<1) { warning("after filtering, no SNPs remain") ; return(NULL) }
  if(length(final.samps)<1) { warning("after filtering, no samples remain") ; return(NULL) }
  if(genome.order) { snp.ord <- rownames(snp.info); final.snps <- snp.ord[snp.ord %in% final.snps] }
  pre <- Dim(snpMat)
  snpMat <- snpMat[final.samps,final.snps]
  post <- Dim(snpMat)
  if(verbose) { cat("snpMat original size: ",comma(pre),"; new size:",comma(post),"\n") }
  pre <- Dim(sample.info) ; sample.info <- sample.info[final.samps,] ;  post <- Dim(sample.info)
  if(verbose) { cat("sample.info original size: ",comma(pre),"; new size:",comma(post),"\n") }
  #prv(snp.info)
  pre <- Dim(snp.info) ; snp.info <- snp.info[final.snps,] ;  post <- Dim(snp.info)
  if(verbose) { cat("snp.info original size: ",comma(pre),"; new size:",comma(post),"\n") }
  #################################
  ##### CREATE THE @SNPS SLOT #####
  # find alleles, and other support in the snp.info object
  df.info <- ranged.to.data.frame(snp.info,include.cols=TRUE)
  #
  df.info <- column.salvage(df.info,"allele.1",c("allele.1","allele_1","allele-1","a1","allele1",
                  "allele-A","alleleA","allele.a","allele_a","forward","fwd","allele.f","allele.fwd"),ignore.case=TRUE)
  df.info <- column.salvage(df.info,"allele.2",c("allele.2","allele_2","allele-2","a2","allele2",
                  "allele-B","alleleB","allele.b","allele_b","backward","bk","bck","allele.bck"),ignore.case=TRUE)
  # extract optional fields
  suppressWarnings({
    df.info <- column.salvage(df.info,"cM",c("CM","distance","map"),ignore.case=TRUE);
    df.info <- column.salvage(df.info,"snp.name",c("snp.name","snpid","snp.id","dbSNP","id","label","rs.id","rsid"),ignore.case=TRUE)
  })
  snp.rn <- rownames(snp.info)
  snp.ch <- chr2(snp.info)
  if("snp.name" %in% colnames(df.info)) { snp.nm <- df.info[,"snp.name"] } else { snp.nm <- snp.rn }
  if("cM" %in% colnames(df.info)) { snp.cm <- df.info[,"cM"] } else { snp.cm <- rep(NA,length(snp.rn)) }
  snp.ps <- start(snp.info)
  if("allele.1" %in% colnames(df.info)) { snp.a1 <- df.info[,"allele.1"] } else { snp.a1 <- rep(NA,length(snp.rn)) }
  if("allele.2" %in% colnames(df.info)) { snp.a2 <- df.info[,"allele.2"] } else { snp.a2 <- rep(NA,length(snp.rn)) }
  if(!all(length(snp.rn)==c(length(snp.ch),length(snp.nm),length(snp.cm),length(snp.ps),length(snp.a1),length(snp.a2)))) {
    stop("Import of snp.info failed, columns had different lengths") # this should be almost impossible to happen HERE
  } else {
    if(verbose) { cat("@snps slot created successfully\n") }
  }
  snps.slot <- data.frame(chromosome=snp.ch,snp.name=snp.nm,cM=snp.cm,position=snp.ps,allele.1=snp.a1,allele.2=snp.a2)
  rownames(snps.slot) <- snp.rn
  allele.slot <- c("allele.1","allele.2")
  #################################
  #### CREATE THE @SAMPLES SLOT ###
  smp.rn <- rownames(sample.info)
  if("pedigree" %in% colnames(sample.info)) { smp.pd <- sample.info[,"pedigree"] } else { smp.pd <- smp.rn }
  if("member" %in% colnames(sample.info)) { smp.mb <- sample.info[,"member"] } else { smp.mb <- smp.rn }
  if("father" %in% colnames(sample.info)) { smp.ft <- sample.info[,"father"] } else { smp.ft <- rep(NA,length(smp.rn)) }
  if("mother" %in% colnames(sample.info)) { smp.mt <- sample.info[,"mother"] } else { smp.mt <- rep(NA,length(smp.rn)) }
  if("sex" %in% colnames(sample.info)) { smp.sx <- sample.info[,"sex"] } else { smp.sx <- rep(NA,length(smp.rn)) }
  if("plate" %in% colnames(sample.info)) { smp.pl <- sample.info[,"plate"] } else { smp.pl <- rep(NA,length(smp.rn)) }
  smp.ph <- sample.info[,"affected"]
  if(!all(length(smp.rn)==c(length(smp.pd),length(smp.mb),length(smp.ft),length(smp.mt),length(smp.sx),length(smp.ph),length(smp.pl)))) {
    stop("Import of sample.info failed, columns had different lengths") # this should be almost impossible to happen HERE
  } else {
    if(verbose) { cat("@samples slot created successfully\n") }
  }
  samps.slot <- data.frame(pedigree=smp.pd,member=smp.mb,father=smp.ft,mother=smp.mt,sex=smp.sx,affected=smp.ph) #,plate=snp.pl)
  rownames(samps.slot) <- smp.rn
  pheno.slot <- c("affected")
  #################################
  
  ## aSnpMatrix format SPECIFICATION: ##
  #@.Data a SnpMatrix
  #@snps  
  #              chromosome     snp.name cM position allele.1 allele.2
  # imm_1_898835          1 imm_1_898835 NA   898835     <NA>        A
  #@samples
  #          pedigree             member father mother sex affected
  # 2447     2447 120434_A01_BLOOD320736     NA     NA   1        1
  #@phenotype character() [phenotype colname]
  #@alleles = c("allele.1", "allele.2") [allele colnames]
  #######################################
  aSnpMat <- new("aSnpMatrix", .Data = snpMat, snps = snps.slot, 
      samples = samps.slot, alleles = allele.slot, 
      phenotype = pheno.slot)
  return(aSnpMat)
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
    gr <- make.granges(hic$chr,start=hic[[2]],end=hic[[3]])
    gr2 <- make.granges(hic$chr.end,start=hic[[5]],end=hic[[6]])
    #colnames(hic) <- make.names(colnames(hic))
    source("~/github/imputer/hiClass.R")
    hic <- HiC(gr,gr2,hic$N.reads,hic$score,build,hic[[9]],hic[[10]])
  }
  return(hic)    
}

# fmsnps <- read.table("/chiswick/data/ncooper/imputation/T1D/FMSNPS.csv")
#colnames(fmsnps) <- c("chr","start","end","width","strand","band","p","rs.id")
#fmSnps <- data.frame.to.granges(fmsnps)
# pv <- fmSnps[["p"]]

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
#  subsetByOverlaps(otherEnd(cd4),fmSnps)
#within set region for big list
#  # promoters ends in region #
#  reg.limit <- make.granges(chr=15,start=79*10^6,end=80*10^6)
#  subsetByOverlaps(subsetByOverlaps(otherEnd(cd4),reg.limit),fmSnps)
#  # one at a time
#  subsetByOverlaps(otherEnd(cd4),subsetByOverlaps(fmSnps,reg.limit))
#one at a time
#  hit.subset <- subsetByOverlaps(fmSnps,otherEnd(cd4))
#  hit.list <- lapply(as.list(1:length(hit.subset)),function(x) { subsetByOverlaps(otherEnd(cd4),hit.subset[x,]) })
#  ff <- findOverlaps(otherEnd(cd4),fmSnps)
#  ss <- subjectHits(ff)
#  qq <- queryHits(ff)
#  subs <- unique(ss)
#  res.list <- lapply(subs, function(x) { qq[which(ss %in% x)] } )
#  names(res.list) <- paste(subs); res.list <- res.list[order(as.numeric(names(res.list)))]
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
    range.limit <- make.granges(chr=chr.only,start=range.only[1],end=range.only[2])
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
    if(!nonames) { nms <- rn[as.numeric(names(res.list))] ; print(nms) ; names(res.list) <- nms }
    if(index.only) { 
      return(res.list) 
    } else { 
      out <- lapply(1:length(res.list),function(x) { hic[res.list[[x]]] }) 
      names(out) <- names(res.list)
      return(out)
    }
  }
}


# test that it worked! (it did)
#for(cc in 1:length(iii)) {
#  nxt <- names(iii)[cc]
#  fmSnps[nxt,]
#  iii[[nxt]]
#  if(unique(chr(iii[[nxt]]))!=unique(chr(fmSnps[nxt,]))) { stop(cc,"failed")}
#  loop.tracker(cc,length(iii))
#}

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
