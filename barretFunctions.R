library(gtools)

#get all the output RData files from a given chromsome, return as snp.matrix list
get.chr.sml <- function(chr,location="~/barrett/OUTPUT/",common.txt="SnpMatrix") {
  SMLa <- list.files(location,pattern=common.txt,full.names=T)
  SMLa <- SMLa[grep(paste0("-",chr,"-"),SMLa,fixed=T)]
  SML <- as.list(SMLa)
  return(SML)
}


# for aSnpMatrix, a better 'write.impute' function
asm.write.impute <- function(X,fn="testmat22") {
  if(!is(X)[1] %in% c("aSnpMatrix","aXSnpMatrix")) { stop("X must an aSnpMatrix object") }
  acoln <- alleles(X)
  dir <- dirname(fn)
  fn <- basename(fn)
  curdir <- getwd()
  if(!dir %in% c(".","./","",getwd())) { setwd(dir) } # this function doesn't seem to be able handle paths
  poz <- snps(X)$position
  if(any(poz!=sort(poz))) { message("warning: IMPUTE2 will not read files correctly if not sorted by SNP position") }
  if(length(unique(snps(X)$chromosome))>1) { stop("You should only have 1 chromosome per file for this function (or for IMPUTE2)") }
  write.impute(X,a1=snps(X)[[acoln[1]]],a2=snps(X)[[acoln[2]]],bp=snps(X)$position, pedfile=fn, snp.id=colnames(X))
  setwd(curdir) # restore original working directory
}

cohort.alignment.check <- function(snpmat,cohort=NULL,sample.info=NULL) {
  ## enter snpMatrix, and either 'sample.info' or 'cohort' vector
  # checks RAFs of each cohort for consistency, prints report to terminal output #
  typ <- is(snpmat)[1]
  if(!typ %in% c("SnpMatrix","aSnpMatrix","XSnpMatrix","aXSnpMatrix")) { stop("snpmat must be a SnpMatrix or aSnpMatrix") }
  if(!is.null(sample.info) & is.null(cohort)) {
    if(length(Dim(sample.info))!=2) { stop("invalid sample.info object") }
    ii <- match(rownames(samples(snpmat)),rownames(sample.info))
    cohort <- sample.info$cohort[ii]
  } else {
    if(!length(Dim(cohort))==1) { stop("cohort must be a vector") }
    if(!is.null(sample.info)) { stop("cohort overrides argument 'sample.info'")}
    if(nrow(snpmat)!=length(cohort)) { stop("cohort must be the same length as the number of rows in snpmat")}
  }
  ph.tab <- table(cohort,exclude=NULL)
  if(any(is.na(names(ph.tab)))) { message("NAs found in phenotype, corresponding rows will be ignored") }
  phenoz <- narm(names(ph.tab))
  raf.list <- vector("list",length(phenoz))
  for (cc in 1:length(phenoz)) {
    ind <- which(cohort==phenoz[cc])
    if(length(narm(ind))<2) { warning("each phenotype should have at least 2 samples, skipping"); next }
    sm <- snpmat[ind,]
    if(nrow(sm)==0 | ncol(sm)==0) { warning("phenotype had no non-missing data, skipping"); next } 
    rsm <- row.summary(sm)
    if(all(rsm$Call.rate==0)) {  warning("phenotype had no non-missing data, skipping"); next } 
    nxt.sum <- col.summary(sm)
   # prv(nxt.sum)
    raf.list[[cc]] <- nxt.sum$RAF
    if(cc>1) {
      raf.dif <- raf.list[[cc]] - raf.list[[cc-1]]
      cat("----------------\ncohort",phenoz[cc],"versus",phenoz[cc-1],":\n----------------\n")
      cat(length(which(is.na(raf.dif))),"RAFs comparisons were NA (100% missing in at least 1 dataset\n")
      cat(length(which(abs(raf.dif)>=.00 & abs(raf.dif)<.05)),"RAFs differed by less than 0.05\n")
      cat(length(which(abs(raf.dif)>=.05 & abs(raf.dif)<.1)),"RAFs differed by 0.05 - 0.10\n")
      cat(length(which(abs(raf.dif)>=.1 & abs(raf.dif)<.5)),"RAFs differed by 0.10 - 0.50\n")
      cat(length(which(abs(raf.dif)>=.5 & abs(raf.dif)<.75)),"RAFs differed by 0.50 - 0.75\n")
      cat(length(which(abs(raf.dif)>=.75)),"RAFs differed by more than 0.75\n----------------\n\n")
    }
  }
  return(NULL)
}


# create a phenotype vector to match an  aSnpMatrix
# if the result should be a dependent variable in a logitistic regression
# (typical of phenotype) then DV should equal true to convert 1,2 coded
# data to 0,1 coded data if necessary. If it's just an independent variable
# being imported, then set DV=FALSE.
get.pheno <- function(X,sample.info,verbose=TRUE,col.name="phenotype",DV=TRUE) {
  if(is(X)[1] %in% c("aSnpMatrix","aXSnpMatrix")) {
    ii <- match(rownames(samples(X)),rownames(sample.info))
  } else {
    ii <- match(rownames(X),rownames(sample.info))
  }
  if(!col.name %in% colnames(sample.info)) { stop("col.name ",col.name," not found in sample.info") }
  ph <- (sample.info[,"phenotype"][ii])
  ll <- length(narm(ph))
  if(length(which(narm(as.numeric(paste(ph))) %in% c(1,2)))/ll > .9 & DV) { ph <- ph-1 }
  if(verbose) {  print(table(ph,exclude=NULL))  }
  mmm <- which(is.na(ph))
  if(length(mmm)>0) { 
    if(verbose) { message("found ",length(mmm)," missing values in ",col.name) }
  }   
  return(ph)
}

# basic GWAS #
do.one <- function(X,sample.info=sample.info) {
  # conduct a basic GLM analysis, yield OR and p-value on an aSnpMatrix
  #X <- X[ ,1:100]
  #   ii <- match(rownames(samples(X)),rownames(sample.info))
  #   ph <- sample.info$phenotype[ii]-1
  #   print(table(ph,exclude=NULL))
  #   mmm <- which(is.na(ph))
  #   if(length(mmm)>0) { 
  #     message("found ",length(mmm)," missing values in phenotype")
  #   } 
  ph <- get.pheno(X,sample.info)
  dat <- aSnpMatrix.to.df(X)
  snups <- colnames(dat)
  result <- vector("list",length(snups))
  for (ss in 1:length(snups)) {
    fm <- paste("ph ~",snups[ss])
    nxt <- glm(as.formula(fm), family = binomial(logit), data=dat)
    result[[ss]] <- mysumfun(nxt,p.digits=250,ci=T)[[1]]
    #  loop.tracker(ss,ncol(dat))
  }  
  res <- do.call("rbind",args=result)
  colnames(res) <- c("OR","OR-low","OR-hi","p-value")
  res <- as.data.frame(res)
  return(res)
}

# convert aSnpMatrix to data.frame
aSnpMatrix.to.df <- function(X) {
  bb <- as(X@.Data,"matrix")
  dd <- dim(bb)
  ee <- as.numeric(bb)
  ee[ee==0] <- NA
  dim(ee) <- dd
  BB <- as.data.frame(ee)
  cn <- colnames(X); rn <- rownames(X)
  if(!is.null(cn)) { colnames(BB) <- cn }
  if(!is.null(rn)) { rownames(BB) <- rn }
  return(BB)
}

# nice summary from GLM output
mysumfun <- function(glmr,o.digits=3,p.digits=6,lab=TRUE,ci=FALSE)
{
  #return(summary(glmr))
  co <- summary(glmr)$coefficients
  predz <- rownames(co)[-1]
  label <- paste(summary(glmr)$call)[2]
  p <- co[2:nrow(co),4]; #print(p)
  o.r <- exp(co[2:nrow(co),1]); #print(o.r)
  p <- round(p,p.digits)
  
  # outlist <- list(round(o.r,o.digits),p)
  # names(outlist) <- c("OR","p-value")
  if(ci) {
    co1 <- co[2:nrow(co),1]-(1.96*co[2:nrow(co),2])
    co2 <- co[2:nrow(co),1]+(1.96*co[2:nrow(co),2])
    if(sum(co1)<=sum(co2)) { co.l <- co1; co.h <- co2 } else {  co.h <- co1; co.l <- co2 }
    co.l <- exp(co.l); co.h <- exp(co.h)
    out.fr <- cbind(round(o.r,o.digits),round(co.l,o.digits),round(co.h,o.digits),p)
    colnames(out.fr) <- c("OR","OR-low","OR-hi","p-value")
  } else {
    out.fr <- cbind(round(o.r,o.digits),p)
    colnames(out.fr) <- c("OR","p-value")
  }
  if(lab) { out.fr <- list(out.fr); names(out.fr) <- label }
  return(out.fr)
}


# convert the output object of snp.rhs.estimates to a nice results dataframe
estimates.to.results <- function(s1,suf="") {
  beta1 <- unlist(sapply(sapply(s1,"[",1),"[",1))
  Var.beta1 <- sqrt(unlist(sapply(sapply(s1,"[",2),"[",1)))
  N1 <- unlist(sapply(sapply(s1,"[",3)  ,"[",1))
  flat <- as(as(s1,"GlmTests"),"GlmTestsScore")  
  pp <- p.value(flat)
  nn1 <- gsub(".N","",names(N1),fixed=TRUE)
  rezult1 <- cbind(exp(beta1),Var.beta1,N1)
  rownames(rezult1) <- nn1
  #print(head(names(flat)))
  #prv(rezult1,nn1,flat,names(flat),pp)
  ppp <- pp[match(nn1,names(flat))]
  #prv(ppp,nn1,pp)
  #return(list(ppp=ppp,rezult1=rezult1))
  rezult1 <- as.data.frame(rezult1) 
  rezult1[["p.value"]] <- ppp
  colnames(rezult1) <- paste(c("OR","SE","N","p"),suf,sep="")
  return(rezult1)
}


# pool the results of two gwas analyses, using meta.me()
my.pool <- function(s1,s2,method=c("beta","z.score","sample.size")[2]) {
  beta1 <- unlist(sapply(sapply(s1,"[",1),"[",1))
  Var.beta1 <- sqrt(unlist(sapply(sapply(s1,"[",2),"[",1)))
  N1 <- unlist(sapply(sapply(s1,"[",3)  ,"[",1))
  nn1 <- gsub(".N","",names(N1),fixed=TRUE)
  rezult1 <- cbind(exp(beta1),Var.beta1,N1)
  #prv(nn1,beta1,N1)
  rownames(rezult1) <- nn1
  beta2 <- unlist(sapply(sapply(s2,"[",1),"[",1))
  Var.beta2 <- sqrt(unlist(sapply(sapply(s2,"[",2),"[",1)))  # why is this needed?
  N2 <- unlist(sapply(sapply(s2,"[",3)  ,"[",1))
  nn2 <- gsub(".N","",names(N2),fixed=TRUE)
  # if(length(beta1)!=length(beta2)) { stop("betas had different lengths")}
  # if(length(N1)!=length(N2)) { stop("betas had different lengths")}
  # if(length(Var.beta1)!=length(Var.beta2)) { stop("betas had different lengths")}
  rezult2 <- cbind(exp(beta2),Var.beta2,N2)
  rownames(rezult2) <- nn2
  all.pair.names <- unique(c(nn1,nn2)[duplicated(c(nn1,nn2))])
  rezult1 <- rezult1[all.pair.names,]
  rezult2 <- rezult2[all.pair.names,]
  #prv(rezult1,rezult2)
  comb <- data.frame(OR1=rezult1[,1],SE1=rezult1[,2],N1=rezult1[,3],OR2=rezult2[,1],SE2=rezult2[,2],N2=rezult2[,3])
  print(head(comb))
  return(meta.me(comb,OR1="OR1",OR2="OR2",SE1="SE1",SE2="SE2",N1="N1",N2="N2",method=method))
}

# create a vector to select any samples from a list in a given dataset
sampsIn <- function(X,samps,snps=FALSE) {
  if(is.null(colnames(X)) & is.null(rownames(X))) { rn <- names(X) } else {
    if(snps) { rn <- colnames(X) } else { rn <- rownames(X) }
  }
  return(narm(match(samps,rn)))
}

# create a vector to select any SNPs from a list in a given dataset
snpsIn <- function(X,snps) { sampsIn(X,snps,snps=TRUE) }


## impute2 1000gn outputted names in my dataset have a set format allowing a snp.info object to
# be extrapolated from the rownames. this function automates this extrapolation.
impute.rn.to.info <- function(nms,sep=":") {
  if(!is.null(dim(nms))) { nms <- colnames(nms) }
  #prv(nms)
  ll <- strsplit(nms,split=sep,fixed=TRUE)
  pos <- sapply(ll,"[",2)
  a1 <- sapply(ll,"[",3)
  a2 <- sapply(ll,"[",4)
  id <- sapply(ll,"[",1)
  tt <- table(id)
  tt <- tt[names(tt) %in% paste(1:30)]
  prob.chr.name <- names(tt)[which(tt==max(tt,na.rm=T))]
  chr <- rep(as.numeric(prob.chr.name),length(pos))
  id[id %in% paste(1:30)] <- NA
  #prv(chr,pos,a1,a2,id)
  new.frame <- data.frame(chr=chr,pos=pos,allele.1=a1,allele.2=a2,rs.id=id)
  rownames(new.frame) <- nms
  return(new.frame)
}


## give a directory and chromosome number and this will extract all of the relevant impute output file names
get.impute.chr.fn <- function(dir,chr=22,full.names=TRUE,info.only=FALSE,by.sample=FALSE) {
  ii <- list.files(dir,pattern=paste0("impute-",chr,"-"),full.names=full.names)
  if(info.only) {
    if(by.sample) {
      ii <- ii[c(grep("info_by_sample",ii))]
    } else {
      ii <- ii[-c(grep("sample",ii),grep("warnings",ii),grep("summary",ii))]
      ii <- ii[c(grep("out_info",ii))]
    }
  } else {
    mm <- c(grep("sample",ii),grep("info",ii),grep("warnings",ii),grep("summary",ii))
    if(length(mm)>0) { ii <- ii[-mm] }
    #ii <- ii[-c(grep("sample",ii),grep("info",ii),grep("warnings",ii),grep("summary",ii))]
  }
  if(length(grep("SnpMatrix",ii))>0) { ii <- ii[-c(grep("SnpMatrix",ii))] }
  ii <- gtools::mixedsort(ii)
  return(rev(ii))
}


# starting with a filename and a sample info file for the whole dataset, read in an impute file (assuming for 1-5MB) and 
# run simple GWAS OR/SE analysis for each SNP, optionally returning with annotation too.
# ... further arguments to snp.rhs.estimates, e.g, subset=<SAMPLE_SUBSET>
analyse.impute.file <- function(fn,sample.info=NULL,add.info=TRUE,add.interp=FALSE,samp.subset=NULL,smp="~/barrett/sampleIDs.txt",
                                save.snpmat=TRUE,restore=TRUE,nPCs=0,PC.source.fn="/chiswick/data/ncooper/barrett/result.quick.A.RData") { 
  if(is.null(sample.info)) { stop("must provide sample.info") }
  #smp <- readLines("~/barrett/sampleIDs.txt")
  if(!file.exists(smp)) { stop("could not read sample file",smp) } 
  smp <- readLines(smp)
  new.fn <- cat.path(dirname(fn),rmv.ext(basename(fn),F),suf="SnpMatrix",ext="RData")
  #oo <- file.exists(new.fn); prv(new.fn,oo,restore)
  if(file.exists(new.fn) & restore) { 
    s1 <- reader(new.fn); cat("read",new.fn,"\n") 
  } else {
    #smp <- readLines("~/barrett/sampleIDs.txt")
    success <- T
    success <- tryCatch(s1 <- read.impute(fn,rownames=smp),error=function(e) { F } )
    if(!is.logical(success)) { success <- T }
    if(!success) { warnings("read.impute failed for file: ",fn); stop("impute file not read") }
    #s1 <- read.impute(fn,rownames=smp)
    if(any(nchar(colnames(s1))>127)) { colnames(s1) <- substr(colnames(s1),1,127) }
    if(save.snpmat) {
      X <- s1
      save(X,file=new.fn)
      rm(X)
    }
  }
  if(!is(s1)[1]=="SnpMatrix") { message("read failed for file: ",new.fn) }
  ph <- get.pheno(s1,sample.info,verbose=FALSE)
  sdat <- as(s1,"SnpMatrix")
  if(any(nchar(colnames(sdat))>127)) { colnames(sdat) <- substr(colnames(sdat),1,127) }
  if(!is.null(samp.subset)) {
    iii <- which(rownames(sdat) %in% samp.subset)
    ph <- ph[iii] ; sdat <- sdat[iii,]
    #print(Dim(sdat))
  } else {
    # print("no subsetting")
  }
  if(nPCs>0 | T) {
    if(!file.exists(PC.source.fn)) { stop("nPCs was set >0 but PC file 'PC.source.fn' didn't exist") }
    rqa <- reader(PC.source.fn)
    PCs <- rqa$PCs
    mtch <- (rownames(sdat) %in% rownames(PCs))
    if(length(which(!mtch))>length(which(mtch))) { stop("most samples in matrix were not in the PCs element of the PCA specified by argument PC.source.fn") }
    sdat <- sdat[mtch,]; ph <- ph[mtch]
    myDat <- as.data.frame(cbind(ph,PCs[rownames(sdat),]))
    #prv(myDat,PCs)
    #print(table(ph))
    if(nPCs==0) {
      ## this here just to be sure we're comparing equal n's
      estz2 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=sdat, uncertain = TRUE )
    } else {
      form <- paste0("ph ~ ",paste(colnames(myDat)[-1][1:nPCs],collapse="+"))
      estz2 <- snp.rhs.estimates(as.formula(form), family = "binomial", snp.data=sdat, data=myDat, uncertain = TRUE)    
    }
  } else {
    estz2 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=sdat, uncertain = TRUE )
  }
  #prv(sdat)
  res <- estimates.to.results(estz2)
  if(add.info) {
    info <- impute.rn.to.info(sdat)
    #prv(res,info)
    res <- cbind(res,info[rownames(res),])
  }
  if(add.interp) {
    # do caseway and majmin, adjust directions for ORs
    res <- add.dirs.to.result(res,sdat,ph)
  }
  #prv(res)
  return(res)
}

# just like 'analyse.impute.file' but does a comparitive col.summary between phenotypes
summarize.impute.file <- function(fn,sample.info=NULL,samp.subset=NULL,
                                save.snpmat=TRUE,restore=TRUE) {
  if(is.null(sample.info)) { stop("must provide sample.info") }
  #smp <- readLines("~/barrett/sampleIDs.txt")
  new.fn <- cat.path(dirname(fn),rmv.ext(basename(fn),F),suf="SnpMatrix",ext="RData")
  #oo <- file.exists(new.fn); prv(new.fn,oo,restore)
  if(file.exists(new.fn) & restore) {
    s1 <- reader(new.fn); cat("read",new.fn,"\n")
  } else {
    smp <- readLines("~/barrett/sampleIDs.txt")
    s1 <- read.impute(fn,rownames=smp)
    if(any(nchar(colnames(s1))>127)) { colnames(s1) <- substr(colnames(s1),1,127) }
    if(save.snpmat) {
      X <- s1
      save(X,file=new.fn)
      rm(X)
    }
  }
  ph <- get.pheno(s1,sample.info,verbose=FALSE)
  sdat <- as(s1,"SnpMatrix")
  if(any(nchar(colnames(sdat))>127)) { colnames(sdat) <- substr(colnames(sdat),1,127) }
  if(!is.null(samp.subset)) {
    iii <- which(rownames(sdat) %in% samp.subset)
    ph <- ph[iii] ; sdat <- sdat[iii,]
    #print(Dim(sdat))
  } else {
    # print("no subsetting")
  }
  estz1 <- col.summary(sdat[ph==0,])
  estz2 <- col.summary(sdat[ph==1,])
  estz <- cbind(estz1,estz2)
  colnames(estz) <- c(paste(colnames(estz1),"1",sep="."),paste(colnames(estz2),"2",sep="."))
  return(estz)
}



# using a snpmatrix and phenotype vector, run GWAS analysis and add effect direction info
annotated.snp.analysis <- function(snpmat,pheno) {
  ph <- pheno
  estz2 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(snpmat,"SnpMatrix"),uncertain = TRUE)
  res <- estimates.to.results(estz2)
  # do caseway and majmin, adjust directions for ORs
  res <- add.dirs.to.result(rez,snpmat,ph)
  return(res)
}

## add case direction, major minor allele spec, and standardize OR directions for SnpStats analysis results 
add.dirs.to.result <- function(result,s1,pheno) {
  ph <- pheno
  case.way <- caseway(s1,pheno=ph) 
  maj.min <- majmin(s1,pheno=ph) 
  cn <- colnames(s1)
  anno <- cbind(paste(unlist(maj.min)),paste(unlist(case.way)))
  rownames(anno) <- cn
  result[["majmin"]] <- anno[match(rownames(result),rownames(anno)),1]
  result[["caseway"]] <- anno[match(rownames(result),rownames(anno)),2]
  result[["interp"]] <- rep("???",nrow(result))
  result[["interp"]][with(result,caseway=="CasesRef+" & majmin=="major")] <- "T1D Major+"
  result[["interp"]][with(result,caseway=="CasesRef+" & majmin=="minor")] <- "T1D Minor+"
  result[["interp"]][with(result,caseway=="CasesRef-" & majmin=="minor")] <- "T1D Major+"
  result[["interp"]][with(result,caseway=="CasesRef-" & majmin=="major")] <- "T1D Minor+"
  result[["OR_orig"]] <- result$OR
  result[["OR"]][result$OR>1 & result$interp=="T1D Major+"] <- 1/(result[["OR"]][result$OR>1 & result$interp=="T1D Major+"])
  result[["OR"]][result$OR<1 & result$interp=="T1D Minor+"] <- 1/(result[["OR"]][result$OR<1 & result$interp=="T1D Minor+"])
  return(result)
}


do.a.qq <- function(p.values,suf="") {
  ofn <- cat.path("/chiswick/data/ncooper/barrett","chiqq",ext="pdf",suf=suf)  
  pdf(ofn) #cat.path("/chiswick/data/ncooper/barrett","chiqq",ext="pdf",suf=suf))
  #prv(p.meta)
  #p.meta <- narm(p.meta)
  txt <- paste("Imputed GWAS (ChiSq) [",toheader(paste(gsub("."," ",suf,fixed=T),"regions")),"]",sep="")
  LL <- length(which(!is.na(qchisq(p.values,1))))
  xx <- qchisq(1-((1:LL)/LL),1); yy <- qchisq(1-sort(p.values),1)
  cond <- (is.na(xx) | is.na(yy) | !is.finite(yy) | !is.finite(xx))
  xx <- xx[!cond]; yy <- yy[!cond]
  three4 <- function(x) { mean(c(min(x,na.rm=T), rep(max(x,na.rm=T),3))) }
  plot(xx,yy,type="l",main=txt,xlab="expected",ylab="observed",xlim=c(0,20))
  text(three4(xx),three4(yy),paste("slope =",round(coefficients(lm(yy~xx))[2],3)))
  lines(x=qchisq(1-((1:LL)/LL),1),y=qchisq(1-((1:LL)/LL),1),col="red",lty="dotted")
  dev.off()
  cat("wrote QQ-plot to:",ofn,"\n")
}


### 1000 genomes extract to aSnpMatrix functions ##

# FOR PHASE 3 in LEG, HAP format:
## internal ##
is.valid.rn <- function(x) {
  return( if(!is.null(x)) { (!all(paste(x)==paste(1:length(x)))) & !anyDuplicated(x) } else { F } )
}

## internal ##
get.vcf.header.from.files <- function(leg,samp,input.dir,work.dir,just.count=FALSE) {
  leg <- cat.path(input.dir,leg)
  samp <- cat.path(input.dir,samp)
  smps <- reader(samp)
  try.rn <- rownames(smps)
  if(!is.valid.rn(try.rn)) { rn <- paste(smps[[1]]); if(!is.valid.rn(rn)) { stop("couldn't find valid sample ids in 'samp'") } } else { rn <- try.rn } 
  tmp.fn <- "tmphead.txt"
  system(paste0("zcat ",leg," | head -1 > ",tmp.fn))
  if(!file.exists(tmp.fn)) { stop("expected temporary file containing legend headings did not exist") }
  leg.head <- strsplit(readLines(tmp.fn)," ",fixed=T)[[1]]
  if(just.count) { return(length(leg.head)) }
  unlink(tmp.fn)
  header <- paste(c(leg.head,rn),collapse=" ")
  return(header)
}

## convert SNPTEST/IMPUTE leg/hap/sample format to a vcf file (requires bash) ##
convert.leg.hap.samp.to.vcf <- function(leg="1000GP_Phase3_chr22.legend.gz",hap="1000GP_Phase3_chr22.hap.gz",
                                        samp="1000GP_Phase3.sample", out.fn="TESTmyChipChr22TEST", input.dir="~/barrett/1000GP_Phase3/", work.dir="~/PLAY/") 
{
  reqd.cmds <- c("sed","zcat","paste","cat","tail","head")
  if(any(!check.linux.install(reqd.cmds)))  { stop("some commands missing of: ",reqd.cmds,) }
  leg <- cat.path(input.dir,leg)
  hap <- cat.path(input.dir,hap)
  samp <- cat.path(input.dir,samp)
  out.fn1 <- cat.path(work.dir,out.fn,ext="tmp")
  out.fn <- cat.path(work.dir,out.fn,ext="vcf")
  out.tmp <- cat.path(work.dir,"myTEmpSDF",ext="txt")
  smps <- reader(samp)
  header <- get.vcf.header.from.files(leg,samp,input.dir,work.dir)
  cmd <- paste0("paste -d ' ' <(zcat ",leg," | tail -n +2) <(zcat ",hap,
                " | sed 's/\\(.\\) \\(.\\)/\\1\\2/g' | sed 's/00/0|0/g' | sed 's/11/1|1/g' ",
                "| sed 's/01/0|1/g' | sed 's/10/1|0/g') > ",out.fn1,"\n")
  writeLines(header,con=out.tmp)
  cmd <- paste0("bash -c \"",cmd,"\"")
  cmd2 <- paste0("cat ",out.tmp," ",out.fn1," > ",out.fn,"\n")
  cat(cmd)
  system(cmd)
  if(!file.exists(out.fn1)) { stop("bash call to 'sed' and 'zcat' failed") }
  system(cmd2)
  unlink(out.tmp) # clean up files created in the working directory
  unlink(out.fn1)
  return(out.fn)
}


# ## code to get from the crappy vcf format outputted by my R function to
# # a real vcf format recognisable by the nice perl conversion script on the IMPUTE2 website
# 
# # create some dummy files of the correct length using R, to be added into the file to comply with vcf spec
# fn <- ("testChr22.vcf")
# nn <- file.nrow(fn)
# writeLines(paste(rep(22,nn)),con="chr.txt")
# writeLines(paste(rep("PASS",nn)),con="filt.txt")
# writeLines(paste(rep("GT",nn)),con="format.txt")
# writeLines(paste(rep(".",nn)),con="dots.txt"
#            
# ## replace '|' with '/'
# "sed -i 's/|/\//g' testChr22.vcf"
#            
# # paste in all the necessary columns
# "paste -d ' ' chr.txt <(cut -f 2 -d ' ' testChr22.vcf) <(cut -f 1 -d ' ' testChr22.vcf) <(cut -f 3-4 -d ' ' testChr22.vcf) dots.txt filt.txt dots.txt format.txt <(cut -f 12- -d ' ' testChr22.vcf) > big.vcf"
#            
# # change to tab delimited
# "sed -i 's/ /\t/g' big.vcf"
#            
# # convert vcf to .gen format
# "./vcf2impute_gen.pl -vcf big.vcf -gen realbig.gen"

           
# part 1 to extract 1000 genomes data
#e.g my.args <- list(leg="1000GP_Phase3_chr22.legend.gz",hap="1000GP_Phase3_chr22.hap.gz",
#              samp="1000GP_Phase3.sample", out.fn="MAINmyChipChr22MAIN", input.dir="~/barrett/1000GP_Phase3/", #work.dir="~/PLAY/",snp.info=snp.info,chr=22) 
#list.to.env(my.args)

# .smt = snpMatrixText
## extract aSnpMatrix from 1000 genomes SNPTEST/IMPUTE leg/hap/sample format ##
extract.leg.hap.samp.1000g.p1 <- function(leg="1000GP_Phase3_chr22.legend.gz",hap="1000GP_Phase3_chr22.hap.gz",
                                          samp="1000GP_Phase3.sample", out.fn="myChipChr22", input.dir="~/barrett/1000GP_Phase3/", work.dir="~/PLAY/",snp.info=snp.info,chr=22) 
{
  chr.pos.fn <- cat.path(work.dir,"gwasSnpPos",suf=chr,ext="txt")
  writeLines(paste(start(chrSel(snp.info,chr))),con=chr.pos.fn)
  reqd.cmds <- c("sed","zcat","paste","cat","tail","head","awk")
  if(any(!check.linux.install(reqd.cmds)))  { stop("some commands missing of: ",reqd.cmds,) }
  leg <- cat.path(input.dir,leg)
  hap <- cat.path(input.dir,hap)
  samp <- cat.path(input.dir,samp)
  out.fn1 <- cat.path(work.dir,out.fn,ext="tmp")
  out.fn <- cat.path(work.dir,out.fn,ext="smt")
  out.tmp <- cat.path(work.dir,"myTEmpSDF",suf=chr,ext="txt")
  smps <- reader(samp)
  header <- get.vcf.header.from.files(leg,samp,input.dir,work.dir)
  cnt <- get.vcf.header.from.files(leg,samp,input.dir,work.dir,just.count=T) # get size of legend
  cmd <- paste0("/usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ",chr.pos.fn," <(",
                "paste -d ' ' <(zcat ",leg," | tail -n +2) <(zcat ",hap,
                " | sed 's/\\(.\\) \\(.\\)/\\1\\2/g' | sed 's/00/3/g' | sed 's/11/1/g' ",
                "| sed 's/01/2/g' | sed 's/10/2/g')) > ",out.fn1)
  writeLines(header,con=out.tmp)
  #cmd <- paste0("bash -c $(",cmd,")","\n")
  cat("\n",cmd,"\n\nCopy and paste this command ^ into your terminal and do not press enter until it has completed")
  return(list(leg=leg,hap=hap,samp=samp,out.fn=out.fn,out.fn1=out.fn1,out.tmp=out.tmp,cnt=cnt,chr=chr,work.dir=work.dir))
}

# part 2 to extract 1000 genomes data
extract.leg.hap.samp.1000g.p2 <- function(args.list=NULL) {
  if(is.null(args.list)) { stop("args.list should be a list outputted by 'extract.leg.hap.samp.1000g.p1'") }
  list.to.env(args.list)
  
  #system(cmd)
  if(!file.exists(out.fn1)) { stop("bash call to 'sed' and 'zcat' failed") }
  cmd2 <- paste0("cat ",out.tmp," ",out.fn1," > ",out.fn,"\n")
  system(cmd2)
  
  dir.working <- work.dir
  
  ii <- read.table(out.fn,header=T)   # 1min30s???
  
  ## NEEDS TESTING FROM ABOUT HERE ###
  support <- (ii[,1:cnt])
  support[["index"]] <- 1:nrow(support)
  support[["CHROM"]] <- rep(chr,nrow(support))
  sup <- df.to.ranged(support,chr="CHROM",start="position",end="position")
  #sup[["ref"]] <- substr(sup[["ref"]],1,20)
  ic <- snp.info
  sup2 <- subsetByOverlaps(as(sup,"GRanges"),as(ic,"GRanges"))
  
  ind <- narm(match(mcols(sup2)[["id"]],ii$ID))
  
  unid.tg <- paste(chrm(sup2),start(sup2),sep="-")
  unid.ic <- paste(chrm(ic),start(ic),sep="-")
  
  ind <- match(unid.tg,unid.ic)
  sup2 <- sup2[!duplicated(ind),]
  ind <- ind[!duplicated(ind)]
  ichip.ids <- rownames(ic)[ind]
  rownames(sup2) <- ichip.ids
  
  sup2 <- toGenomeOrder(sup2)
  dat <- ii[mcols(sup2)[["index"]],-c(1:9)]
  rownames(dat) <- rownames(sup2)
  dat <- t(dat)
  TG <- df.to.SnpMatrix(dat)
  
  si <- reader("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData")
  si[["affected"]] <- rep(1,nrow(si))
  
  # ref and alt switched because we coded 0|0 as 3, as 0 is the ref allele, for SnpMatrix, dat is # of ref alleles
  mcols(sup2)[["allele.1"]] <- mcols(sup2)[["a1"]]  # alt 
  mcols(sup2)[["allele.2"]] <- mcols(sup2)[["a0"]]  # ref
  
  asm <- annot.sep.support(TG, sup2, si)
  
  # now save 'asm', this is your annotated 1000 genomes aSnpMatrix object
  # that you will need to align your dataset to.
  unlink(out.tmp) # clean up files created in the working directory
  unlink(out.fn1)
  
  #save(asm,file="thousandGenomesSnpMatrix.RData")
  return(asm)
}


check.align.plot <- function(x,y,fn=NULL,...) {
  if(!is(x)[1] %in% c("aSnpMatrix","aXSnpMatrix")) { stop("x must be an aSnpMatrix") } 
  if(!is(y)[1] %in% c("aSnpMatrix","aXSnpMatrix")) { stop("y must be an aSnpMatrix") } 
  x.cs <- col.summary(x)
  y.cs <- col.summary(y)
  if(is.character(fn)) { pdf(fn,...) }
  plot(x.cs[, "RAF"], y.cs[, "RAF"], main = "RAF after switching", 
       xlab = "x", ylab = "y", pch = "+")
  abline(0, 1, col = "red", lty = 1)
  if(is.character(fn)) { dev.off() }
}


# convert a snpmatrix containing uncertain genotypes to a standard SnpMatrix #
convert.snpmat.uncertain.to.normal <- function(X) {
  Y <- matrix(as.raw(apply(pp(X),1,function(x) { which(x==max(x,na.rm=T))[1] })),nrow=nrow(X))
  rownames(Y) <- rownames(X); colnames(Y) <- colnames(X)
  return(Y)
}




one.row <- function(x) { 
  hrs <- days <- mins <- 0
  mn.days <- c(jan=31,feb=28,mar=31,apr=30,may=31,jun=30,jul=31,aug=31,sep=30,oct=31,nov=30,dec=31)
  #Sep  20 01:45 BIG15_16IllA100.sh Sep  20 09:57 BIG15_16IllA100.sh.o5126484
  #1    2   3    4                   5    6   7   8       
  if(x[5]!=x[1]) days <- days + mn.days[tolower(x[1])]
  if(x[6]!=x[2]) days <- days + (as.numeric(x[6])-as.numeric(x[2]))
  ss <- strsplit(paste(x[c(3,7)]),":",fixed=T)
  hrs <- as.numeric(ss[[2]][1]) - as.numeric(ss[[1]][1])
  mins <- as.numeric(ss[[2]][2]) - as.numeric(ss[[1]][2])
  hrs <- hrs + (days*24) + (mins/60)
  return(hrs)
}


nms.x <- function(x) {
  ss <- strsplit(paste(x[4]),"BIG|[_]|Ill|[.]")[[1]]
  ss <- ss[!ss %in% c("","sh","BIG")]
  out <- c(NA,NA,NA,NA)
  out[1:length(ss)] <- ss
  out <- out[1:4]
  return(out)
}

#ii <- reader("timesum.txt")

get.timings <- function(read){
  if(length(Dim(read))==2) { ii <- read } else { stop("invalid 'read'") } 
  if(any(read$Mon=="6Sep")) { read$Mon <- gsub("6Sep","Sep",read$Mon,fixed=T) } 
  if(any(read$Mon=="3Sep")) { read$Mon <- gsub("3Sep","Sep",read$Mon,fixed=T) }
  if(any(read$Mon=="8Sep")) { read$Mon <- gsub("8Sep","Sep",read$Mon,fixed=T) }
  if(any(read$Mon=="6Oct")) { read$Mon <- gsub("6Oct","Oct",read$Mon,fixed=T) }
  if(any(read$Mon=="3Oct")) { read$Mon <- gsub("3Oct","Oct",read$Mon,fixed=T) }
  if(any(read$Mon=="8Oct")) { read$Mon <- gsub("8Oct","Oct",read$Mon,fixed=T) }
  ii <- ii[order(ii$Name),]
  nm <- ii$Name
  pref <- sapply(strsplit(nm,".",fixed=T),head,n=1)
  suf <- get.ext(nm)
#prv(ii,nm)  #,pref,suf)
  pre <- ii[suf=="sh",]
  post <- ii[suf!="sh",][match(pref[suf=="sh"],pref[suf!="sh"]),]
  big <- cbind(pre,post)
  big <- narm(big)

  #orig.big <- big 
  hrs <- apply(big,1,one.row)
  labs <- t(apply(big,1,nms.x))
  #prv(labs,big)
  colnames(labs) <- paste0("nm",1:ncol(labs))
  big[["hrs"]] <- hrs
  big <- cbind(big,labs)
  oo <- sapply(big$Name,readLines,n=1)
  info <- t(sapply(oo, extract.info.from.cmd))
  
#  prv(big,info)
  big <- cbind(big,info)
  return(big)
}



extract.info.from.cmd <- function(cmd) {
  if(is.character(cmd)) { oo <- cmd } else { stop("cmd must be a character string") }
  ss <- strsplit(oo," ")[[1]]
  ssi <- which(ss=="-int") ; if(length(ssi)>0) { lr <- ss[ssi+1]; hr <- ss[ssi+2] } else { lr <- hr <- NA } 
  ssi <- which(ss=="-sample_g") ; if(length(ssi)>0) { smp <- ss[ssi+1]; smp <- basename(smp) } else { smp <- NA } 
  ssi <- which(ss=="-exclude_samples_g") ; if(length(ssi)>0) { excl <- ss[ssi+1]; excl <- basename(excl) } else { excl <- NA } 
  ssi <- which(ss=="-g") ; if(length(ssi)>0) { dat <- ss[ssi+1]; dat <- basename(dat) } else { dat <- NA } 
  chr <- tail(strsplit(dat,"-",fixed=T)[[1]],1)
  out <- c(dat=dat,excl=excl,smp=smp,chr=chr,lr=lr,hr=hr)
  return(out)
}



progress.checkinator <- function(dir) {
  ll <- list.files(dir,pattern=".sh.o",full.names=TRUE)
  sh <- rmv.ext(ll,only.known=F)
  oo <- sapply(sh,readLines,n=1)
# prv(oo)
  info <- t(sapply(oo, extract.info.from.cmd))
 
  filz <- lapply(as.list(ll),readLines) 
  lastz <- sapply(filz,function(x) { paste0(tail(x,2),collapse="") })
  class <- rep("INCOMPLETE",length(lastz))
  class[lastz=="Have a nice day!"] <- "DONE"
  class[substr(lastz,1,12)=="MCMC iterati"] <- "IN_PROGRESS"
  class[substr(lastz,1,12)=="There are no"] <- "EMPTY"
  class[substr(lastz,1,12)=="ERROR: There"] <- "EMPTY"
  class[substr(lastz,1,12)=="Your current"] <- "EMPTY"
  outcome <- class
  out <- cbind(basename(ll),info,outcome)
  rownames(out) <- rmv.ext(rmv.ext(basename(ll),F),F)
  colnames(out)[1] <- "log.file"
  return(out)
}



# get lots of useful information from an IMPUTE_OUTPUT folder from log files
# sometimes coln=6, sometimes =7, needs to be used for some reason
summarise.impute.output.dir <- function(dir,coln=6) {
  cur.dir <- getwd()
  setwd(dir); system(paste0("cd ",dir))
  system(paste0("bash -c \" cat <(echo 'Mon Day Time Name') <(ls -lt | cut -f ",coln,"- -d ' ' | sed 's/[0-9][0-9][0-9]\\ //g' | sed 's/.Sep/Sep/g' | sed 's/\\ \\ /\\ /g' | tail -n +3 | grep -v 'complete_flags' ) > timesum.txt \" "))
  ii <- reader("timesum.txt",header=T) #col.names=T,row.names=F,sep=" ")  
  tt <- get.timings(ii) 
  uu <- progress.checkinator(dir)
  outcome <- uu[,8]
  out <- cbind(tt[match(uu[,"log.file"],tt[,"Name.1"]),],outcome)
  setwd(cur.dir)
  return(out)
}


compare.impute.calls <- function(s1,s2, samps.to.compare=all.affy.ids, 
                                   plot.fn=NULL,uncertain=F) {
  if(uncertain) {
    S1 <- convert.snpmat.uncertain.to.normal(s1)
    S2 <- convert.snpmat.uncertain.to.normal(s2)
  } else { S1 <- s1; S2 <- s2 }
  # sync SNPs
  S2a <- S2[,snpsIn(S2,colnames(S1))]
  S1a <- S1[,snpsIn(S1,colnames(S2a))]
  S2a <- S2a[,colnames(S1a)]
  # sync samples
  S1a <- S1a[sampsIn(S1a, samps.to.compare),]
  S2a <- S2a[sampsIn(S2a, samps.to.compare),]
  S2a <- S2a[sampsIn(S1a,rownames(S2a)),]
  S1a <- S1a[sampsIn(S2a,rownames(S1a)),]
  S2a <- S2a[rownames(S1a),]
  tot <- nrow(S1a); nc <- ncol(S1a); dif <- numeric(nc)
  for (cc in 1:nc) { dif[cc] <- length(which(S1a[,cc]!=S2a[,cc])); loop.tracker(cc,nc) }

  print(summary(dif))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.0    32.0   104.0   309.9   336.0  2542.0 
  print(summary(dif/tot))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.000000 0.004734 0.015390 0.045850 0.049710 0.376100 

  print(cor(col.summary(as(S1a,"SnpMatrix"))$MAF,col.summary(as(S2a,"SnpMatrix"))$MAF,use="pairwise.complete"))  #[1] 0.997459
  if(is.character(plot.fn)) {
    pdf(plot.fn)
    plot(col.summary(as(S1a,"SnpMatrix"))$MAF,
       col.summary(as(S2a,"SnpMatrix"))$MAF,pch=".")
    dev.off()
    send.to.pwf(plot.fn)
  }
}


analyse.dual.snpmatrices <- function(fn,sample.info=NULL,add.info=TRUE,add.interp=TRUE,
   dir1,dir2,prog=FALSE,add.strat=FALSE,verbose=TRUE) {
  COMB <- combine.affy.illu(fn,verbose=verbose)
  ph <- get.pheno(COMB,sample.info,verbose=FALSE)
  sdat <- as(COMB,"SnpMatrix")

  if(add.strat) {
    strat <- factor(get.aff.ill.vec(COMB))
    estz2 <- snp.rhs.estimates(ph ~ strata(strat), family="binomial", snp.data=sdat, uncertain=T)
  } else {
    estz2 <- snp.rhs.estimates(ph ~ 1, family="binomial", snp.data=sdat, uncertain=T)
  }
  res <- estimates.to.results(estz2)
  if(add.info) {
    info <- impute.rn.to.info(sdat)
    #prv(res,info)
    res <- cbind(res,info[rownames(res),])
  }
  if(add.interp) {
    # do caseway and majmin, adjust directions for ORs
    res <- add.dirs.to.result(res,sdat,ph)
  }
  if(prog) { cat("*") }
  #prv(res)
  return(res)
}


best.match.to.common.snps <- function(AFF,ILL,verbose=TRUE,print.fails=FALSE,row.names=FALSE) {
  if(length(Dim(AFF))!=2) { stop("AFF must have 2 dimenions") }
  if(length(Dim(ILL))!=2) { stop("ILL must have 2 dimenions") }
  if(row.names) { AFF <- t(AFF); ILL <- t(ILL) }
  ind <- intersect(colnames(AFF),colnames(ILL))
  # established perfect matches above
  # now below find those that match the rs, e.g, ignore latter of rsid:...:...
  caf <- colnames(AFF)[!colnames(AFF) %in% ind]
  cil <- colnames(ILL)[!colnames(ILL) %in% ind]
  caf2 <- sapply(strsplit(caf,":",fixed=T),"[",1)
  cil2 <- sapply(strsplit(cil,":",fixed=T),"[",1)
  ind2 <- intersect(caf2,cil2)
  caf3 <- caf[match(ind2,caf2)]
  cil3 <- cil[match(ind2,cil2)]
  # in each case, select the longer of the 2 names as the 'official' version
  ind3 <- apply(cbind(caf3,cil3),1,function(x) { if(nchar(x[2])>nchar(x[1])) { x[2] } else { x[1] } })
  # change the names in both files to match the set of common longer names
  colnames(AFF)[match(caf3,colnames(AFF))] <- ind3
  colnames(ILL)[match(cil3,colnames(ILL))] <- ind3
  # run the intersection again
  ind <- intersect(colnames(AFF),colnames(ILL))
  #optionally display those failing the matching (unique in one set or the other)
  if(print.fails) {
    caf <- colnames(AFF)[!colnames(AFF) %in% ind]
    cil <- colnames(ILL)[!colnames(ILL) %in% ind]
    if(length(caf)>0 | length(cil)>0) { print(rbind(sort(caf),sort(cil))) }
  }
  if(verbose){ cat("#affy:",ncol(AFF),"; #illu:",ncol(ILL),"; common:",length(ind),"\n") }
  AFF <- AFF[,ind]
  ILL <- ILL[,colnames(AFF)]
  if(row.names) { AFF <- t(AFF); ILL <- t(ILL) }
  return(list(AFF,ILL))
}


get.common.aff.ill.names <- function(chr=22,full.names=FALSE, full.dir1=TRUE,
                        dir1="OUTPUT_AFFY",dir2="OUTPUT_ILLU") {
  ii <- list.files(dir1,pattern=paste0("impute-",chr,"-"),full.names=full.names)
  af.fn <- ii[c(grep("SnpMatrix",ii))]
  ii <- list.files(dir2,pattern=paste0("impute-",chr,"-"),full.names=full.names)
  il.fn <- ii[c(grep("SnpMatrix",ii))] 
  if(full.dir1) { oo <- c(il.fn,af.fn) } else { oo <- c(af.fn,il.fn) }
  oo <- oo[duplicated(basename(oo))]
  return(oo)
}

combine.affy.illu <- function(fn,dir1="OUTPUT_AFFY",dir2="OUTPUT_ILLU",verbose=TRUE) {
  AFF <- reader(cat.path(dir1,fn)); if(verbose) { cat("read",cat.path(dir1,fn),"\n")}
  ILL <- reader(cat.path(dir2,fn)); if(verbose) { cat("read",cat.path(dir2,fn),"\n")}
  lister <- best.match.to.common.snps(AFF,ILL,verbose=verbose,print.fails=FALSE) 
  AFF <- lister[[1]]; ILL <- lister[[2]] # should have equal colnames
  if(any(colnames(AFF)!=colnames(ILL))) { stop("mismatching column names in AFF and ILL even after running 'best.match.to.common.snps()'") }
#  if(verbose){ cat("#affy:",ncol(AFF),"; #illu:",ncol(ILL),"; common:",length(ind),"\n") }
#  AFF <- AFF[,ind]
#  ILL <- ILL[,colnames(AFF)]
  COMB <- rbind(AFF,ILL)
  if(any(nchar(colnames(COMB))>127)) { colnames(COMB) <- substr(colnames(COMB),1,127) }
  return(COMB)
}

get.aff.ill.vec <- function(X,null.val=NA) {
  if(!exists("all.ill.ids") | !exists("all.affy.ids")) { stop("this function requires both global measures all.affy.ids and all.ill.ids to be in the current environment") }
  x <- rownames(X)
  grp <- rep(null.val,length(x)) 
  grp[x %in% all.ill.ids] <- 2
  grp[x %in% all.affy.ids] <- 1
  return(grp)
}


snp.test.chr <- function(chrz=21:22,base.dir="~/barrett", prog.dir="~/barrett/SNPTEST", out.dir="~/barrett/SNPTEST/SNPTEST_OUTPUT/",samp.fn="barrett.sample", snp.test.cmd="~/snptest_v2.5.1_linux_x86_64_static/snptest_v2.5.1", excl.fn=NULL,out.fn="snptestResultsChr", imp.dir="~/barrett/OUTPUT/", log.pref="SNPtst", bayesian=FALSE,method=if(bayesian) { "score" } else { "em" }, covariates=FALSE, max.conc=10, stagger=30,grid.name="eightcpu") {		
 all.L1000 <- numeric(22)

 dirb <- prog.dir

 unlink(list.files(cat.path(dirb,"complete_flags"),full.names=T))
 for (cchr in chrz) {
  Header(cchr)
  # frequentist analysis by SNPTEST of 1 chr's impute output files #
  fnz <- get.impute.chr.fn(imp.dir,chr=cchr,full.names=T)

  fn.s <- cat.path(dirb, samp.fn)
  lf <- length(fnz)
  system(paste("cd", base.dir)); setwd(base.dir)
  kk <- proc.time(); 
  cmd.list <- vector("list",lf)
  for (cc in 1:lf) {
    cmd <- character(3)
    fn.o <- fnz[cc]
    fn.g <- cat.path("",rmv.ext(fn.o,F),ext="gen")
    fn.o1 <- cat.path(out.dir,basename(rmv.ext(fn.o,F)),ext="stf")
    fn.o2 <- cat.path(out.dir,basename(rmv.ext(fn.o,F)),ext="stb")
    cmd[1] <- paste0("mv ",fn.o," ",fn.g,"\n")
    met.txt <- paste0(if(bayesian) { " -bayesian 1" } else { " -frequentist add" }," -method ",method)
    cmd[2] <- paste0(snp.test.cmd," -data ",fn.g," ",fn.s,"  -o ",fn.o1,met.txt," -pheno phenotype ",if(is.character(excl.fn)) { paste("-exclude_samples",excl.fn) } else { "" },if(covariates) {" -cov_all_discrete "} else {""}."\n") 
    cmd[3] <- paste0("mv ",fn.g," ",fn.o,"\n")
    cmd.list[[cc]] <- cmd
    #loop.tracker(cc,lf,st.time=kk)
  }

  all.cmds <- sapply(cmd.list,function(x){ paste(x,collapse="; ") })
  bash.qsub(all.cmds,dir=dirb,grid.name=grid.name,logpref=paste0(log.pref,"_",cchr,"_"),stagger=stagger,max.conc= max.conc) # in prog

  freq.files <- cat.path(out.dir,basename(rmv.ext(fnz,F)),ext="stf") 
  freq.list <- lapply(freq.files,read.table,header=T)
  freq.results <- do.call("rbind",args=freq.list)
  xx <- freq.results$frequentist_add_pvalue

  lambda.nxt <- round(median(p.to.Z(xx)^2,na.rm=T)/.454,3)
  lambda.1000 <- lambda_nm(lambda.nxt,nr=9110,mr=6158)
  cat("Lambda 1000 for chr",cchr,":",lambda.1000,"\n")
  all.L1000[cchr] <- lambda.1000

  save(freq.results,all.L1000,file=cat.path(dirb,fn=out.fn,suf=cchr,ext="RData"))
  unlink(list.files(cat.path(dirb,"complete_flags"),full.names=T))
 }

}



mclapply.retries <- function(X,...,mc.cores=1,fail.after=10,detect.failure=is.null,VERBOSE=TRUE) {
  n <- length(X)
  nparts <- rep(FALSE,n); counter <- 0
  result <- mclapply(X,...,mc.cores=mc.cores)
  fails <- sapply(result, detect.failure)
  nparts[!fails] <- TRUE
  while(!all(nparts)) {
    counter <- counter + 1;
    if(counter>fail.after) { stop("too many failures, stopping") }
    if(VERBOSE ) { cat("running retry #",counter,"for",length(which(!nparts)),"failing jobs with",mc.cores,"cores\n") }
    REZ <- mclapply(X[which(!nparts)],..., mc.cores= mc.cores) 
    prv(REZ)
    result[which(!nparts)] <- REZ
    fails <- sapply(result, detect.failure)
    nparts[!fails] <- TRUE
    if(VERBOSE ) { cat("Job failures for list elements: ",paste(which(!nparts),collapse=","),"\n") }
    if(mc.cores>1) { mc.cores <- mc.cores-1 }
  }
  return(result)
}

