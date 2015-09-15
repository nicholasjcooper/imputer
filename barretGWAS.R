library(annotSnpStats)
library(humarray)
library(gtools)
library(GenomicRanges)
library(reader)
library(snpStatsWriter)
source("~/github/plumbCNV/SnpMatrixList.R")
source("~/github/imputer/barretFunctions.R")

dat.dir <- ("/dunwich/scratch/chrisw/T1DGC/snpStats-genotypes")
raw.dir <- "/ipswich/data/T1DGC/R-objects/genotypes/"
home.dir <- ("/chiswick/data/ncooper/barrett/")

regenerate.snp.sample.lists <- FALSE

cohorts <- c("C58","WTCCC","T1DGC")
sml <- vector("list",length(cohorts)); names(sml) <- cohorts
ll <- list.files(dat.dir,pattern="C58"); oo <- mixedorder(gsub("C58.","",ll)); sml[[1]] <- as.list(cat.path(dat.dir,ll[oo]))
ll <- list.files(dat.dir,pattern="WTCCC"); oo <- mixedorder(gsub("WTCCC.","",ll)); sml[[2]] <- as.list(cat.path(dat.dir,ll[oo]))
ll <- list.files(dat.dir,pattern="T1DGC"); oo <- mixedorder(gsub("T1DGC.","",ll)); sml[[3]] <- as.list(cat.path(dat.dir,ll[oo]))
ll <- list.files(raw.dir,pattern="WT.C58.snp"); oo <- mixedorder(gsub("WT.C58.snp.","",ll)); smlB <- as.list(cat.path(raw.dir,ll[oo]))[1:22]


SML <- as.list(list.files("~/barrett/IMPUTE_INPUT",pattern="snpmat",full.names=TRUE)); oo <- mixedorder(gsub("snpmat","",SML)); SML <- as.list(SML[oo])

if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) } #get sample.info, snp.info, snp.info2



## Obtain SNP and SAMPLE lists ##

if(regenerate.snp.sample.lists) {
  ### SNPS ###
  cs <- list.colsummary(SML)
  cs5 <- smlapply(SML,function(X) { col.summary(X[1:5,]) },n.cores=3)
  mini <- do.call("rbind",args=cs5)
  # all AFFY SNP IDS  
  all.affy <- rownames(mini)[mini$Calls>0]
  # COMMON AFFY/ILL SNP IDS  
  intersec.snps <- rownames(cs)[which(cs$Call.rate>.8)]
  # ILL ONLY IDS   
  ill.only <- rownames(mini)[mini$Calls==0]
  # ILL plus COMMON IDS 
  all.ill <- ill.plus <- c(ill.only,intersec.snps)
  # AFFY ONLY IDS  
  aff.only <- all.affy[-which(all.affy %in% intersec.snps)]
  #prv()
  ### SAMPLES ###
  (load("/ipswich/data/T1DGC/R-objects/genotypes/WT.C58.snp.chr22.RData"))
  (load("/ipswich/data/T1DGC/R-objects/genotypes/sanger.snp.chr22.RData"))
  ILLU <- sanger.22.data
  AFFY <- C58.22.data
  wtccc.sample.support <- reader("/ipswich/data/T1DGC/R-objects/support/WTCCC-sample-support.RData")
  # Ill ids for samples with both datasets   
  ill.int.ids <- wtccc.sample.support$sanger.id[narm(match(rownames(AFFY),wtccc.sample.support$affy.id))]
  # Affy ids for samples with both datasets     # these guys coded with different pheno
  aff.int.ids <- wtccc.sample.support$affy.id[match(ill.int.ids,wtccc.sample.support$sanger.id)]
  # Affy/Ill lookup table to match these ids  
  aff.ill.lookup <- cbind(ill.int.ids,aff.int.ids)
  # ids for all Ills, including dups  all.ill.ids <-
  rs <- list.rowsummary(SML)
  # ids for Affy only ids, 
  ill.all.ids <- rownames(rs)[which(rs$Call.rate>.55)]
  ill.only.ids <- ill.all.ids[-which(ill.all.ids %in% ill.int.ids)]
  aff.only.ids <- rownames(rs)[which(rs$Call.rate<.55)]
  # all sample ids, final
  all.ids <- rownamesL(SML)
  # all sample ids, incl. dups
  all.ids.dups <- c(all.ids,aff.int.ids)
  # check counts make sense
  prv(ill.int.ids,aff.int.ids,ill.all.ids,ill.only.ids,aff.only.ids,all.ids)
  prv(all.ids.dups,all.affy,intersec.snps,ill.only,all.ill,ill.plus,aff.only)
  save(ill.int.ids,aff.int.ids,ill.all.ids,ill.only.ids,aff.only.ids,all.ids, all.ids.dups,
    all.affy,intersec.snps,ill.only,all.ill,ill.plus,aff.only,file="/chiswick/data/ncooper/barrett/IDLISTS.RData")
} else {
  (load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))
}

## CHECK ALIGNMENT HAS WORKED ##
if(F) {
  ## add sources to sample.info ##
  sample.info[["cohort"]] <- rep(NA,nrow(sample.info))
  sample.info[ill.only.ids,"cohort"] <- "Illumina"
  sample.info[aff.only.ids,"cohort"] <- "Affymetrix"
  sample.info[ill.int.ids,"cohort"] <- "Both"
  
  SML <- as.list(list.files("~/barrett/IMPUTE_INPUT",pattern="snpmat",full.names=TRUE)); oo <- mixedorder(gsub("snpmat","",SML)); SML <- as.list(SML[oo])
  (load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))
  
  tst <- reader(SML[[22]])
  dat <- tst[sampsIn(tst,all.ids),snpsIn(tst,c(intersec.snps))]
  cohort.alignment.check(dat,sample.info=sample.info)
}


### GET the replicated Affy data for replicate samples (only have illumina in Chris' datasets) ###

# extract data for each Chr for SNPs aff.only, group aff.int.ids == > ill.int.ids[lookup]

for (cc in 1:length(smlB)) {
  require(chopsticks)
  nxt <- reader(smlB[[cc]])
  nxt.sm <- as(nxt,"SnpMatrix")
  sm <- nxt.sm[aff.int.ids[narm(match(rownames(nxt.sm),aff.int.ids))],aff.only[narm(match(colnames(nxt.sm),aff.only))]]
  save(sm,file=cat.path(home.dir,"AffyReplicatesChr",suf=cc,ext="RData"))
  loop.tracker(cc,22)
}

# Do alignment of the new Affy data to the current set
# insert into existing snpmatrices to fill in gaps (snpmatXX.RData)
for (cc in 1:length(smlB)) {
  sm <- reader(cat.path(home.dir,"AffyReplicatesChr",suf=cc,ext="RData"))
  SML.nxt <- reader(SML[[cc]])
  ill.idz <- aff.ill.lookup[match(rownames(sm),aff.ill.lookup[,2]),1]
  ii <- which(ill.idz %in% rownames(SML.nxt))
  ind <- match(ill.idz[ii],rownames(SML.nxt))
  jj <- which(colnames(sm) %in% colnames(SML.nxt))
  ind2 <- match(colnames(sm)[jj],colnames(SML.nxt))
  # create alignment
  ind.af <- sampsIn(SML.nxt,aff.only.ids)
  illu.asnp <- SML.nxt[ind.af,ind2]
  affy.mat <- sm[ii,jj]
  # find illumina snpStatsAnnot info, also get affy info and make into aSnpMatrix objects
  affy.asnp <- annot.sep.support(affy.mat, snp.info2 ,sample.info)
  affy.asnp <- affy.asnp[,colnames(illu.asnp)]
  # illu.asnp <- annot.sep.support(illu.mat, snp.info ,sample.info)
  affy.align <- align.alleles(affy.asnp, illu.asnp, mafdiff=0.05)
  sm2 <- as(affy.align,"SnpMatrix") 
  # align.alleles
  # insert aligned values
  SML.nxt[ind,ind2] <- sm2[,]
  dat <- SML.nxt
  save(dat,file=cat.path(home.dir,"AllPlusReplicatesChr",suf=cc,ext="RData"))
  loop.tracker(cc,22)
}

### Analysis Plan ###

# create rules : impute ==>
# (i) ill.only for aff.only.ids, using aff.int.ids, ill.int.ids
# (ii) aff.only for ill.only).ids, using ill.int.ids, aff.int.ids
###############################
# ANALYSE
###############################
# (i) ill.only for all.ill.ids, 
# (ii) aff.only for aff.only.ids, 
# (iii) ill.only for aff.only.ids
# (iv) aff.only for ill.only.ids
# -------------------------------
# (v) intersec.snps for all.ids
# -------------------------------
# using snp.rhs.estimates with imputation.rules passed-in
###############################

## retrieve full chr-wise datasets that have been updated with replicate affy samples ##
SML2 <- as.list(cat.path(home.dir,"AllPlusReplicatesChr",suf=1:22,ext="RData"))

# 1 chromosome at a time, do imputation and analysis
for (dd in 1:22) {
  tst <- reader(SML2[[dd]])
  tst2 <- as(tst,"SnpMatrix")
  # impute from affy to illumina
  pos.to <- snps(tst)$position[snpsIn(tst2,ill.only)]
  pos.fr <- snps(tst)$position[snpsIn(tst2,all.affy)]
  to.impute <- tst2[sampsIn(tst2,ill.int.ids),snpsIn(tst2,ill.only)]
  impute.from <- tst2[sampsIn(tst2,ill.int.ids),snpsIn(tst2,all.affy)]
  imp1 <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
  # impute from illumina to affy
  pos.to <- snps(tst)$position[snpsIn(tst2,aff.only)]
  pos.fr <- snps(tst)$position[snpsIn(tst2,all.ill)]
  to.impute <- tst2[sampsIn(tst2,ill.int.ids),snpsIn(tst2,aff.only)]
  impute.from <- tst2[sampsIn(tst2,ill.int.ids),snpsIn(tst2,all.ill)]
  imp2 <- snp.imputation(impute.from, to.impute, pos.fr, pos.to)
  
  if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) } #get sample.info, snp.info, snp.info2
  
  # run analyses separately on 
  # -------------------------------
  ## affy samples imputed to illumina SNPs
  dat <- tst[sampsIn(tst,aff.only.ids),snpsIn(tst,c(ill.only,all.affy))] # add all affy snps just for imputation
  ph <- get.pheno(dat,sample.info)
  estz1 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(dat,"SnpMatrix"),rules = imp1, uncertain = TRUE)
  ##estz1 <- estz1[ill.only]
  ## illumina samples imputed to affy SNPs
  dat <- tst[sampsIn(tst,ill.only.ids),snpsIn(tst,c(aff.only,all.ill))] # add all ill snps just for imputation
  ph <- get.pheno(dat,sample.info)
  estz2 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(dat,"SnpMatrix"),rules = imp2, uncertain = TRUE)
  ##estz2 <- estz2[aff.only]
  ## affy samples real data for affy SNPs
  dat <- tst[sampsIn(tst,c(aff.only.ids,ill.int.ids)),snpsIn(tst,c(aff.only))]
  ph <- get.pheno(dat,sample.info)
  estz3 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(dat,"SnpMatrix"), uncertain = FALSE)
  ## illumina samples real data for illumina SNPs
  dat <- tst[sampsIn(tst,c(ill.only.ids,ill.int.ids)),snpsIn(tst,c(ill.only))]
  ph <- get.pheno(dat,sample.info)
  estz4 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(dat,"SnpMatrix"), uncertain = FALSE)
  ## intersecting SNPs for which all samples have data
  dat <- tst[sampsIn(tst,all.ids),snpsIn(tst,c(intersec.snps))]
  ph <- get.pheno(dat,sample.info)
  estz5 <- snp.rhs.estimates(ph ~ 1, family = "binomial", snp.data=as(dat,"SnpMatrix"), uncertain = FALSE)
  
  estz1a <- as(as(estz1,"GlmTests"),"GlmTestsScore")[snpsIn(estz1,ill.only)]
  estz2a <- as(as(estz2,"GlmTests"),"GlmTestsScore")[snpsIn(tst,aff.only)]
  estz3a <- as(as(estz3,"GlmTests"),"GlmTestsScore")
  estz4a <- as(as(estz4,"GlmTests"),"GlmTestsScore")
  
  ### COMBINE SEPARATE ANALYSES WITH META-ANALYSIS ###
  #illu.snps.results <- pool(estz1a,estz4a)
  illu.snps.results <- my.pool(estz1,estz4,method="sample.size")
  #affy.snps.results <- pool(estz2a,estz3a)
  affy.snps.results <- my.pool(estz2,estz3,method="sample.size")
  joint.snps.results <- as(as(estz5,"GlmTests"),"GlmTestsScore")
  ####################################################
  
  save(estz1,estz2,estz3,estz4,estz1a,estz2a,estz3a,estz4a,estz5,
       illu.snps.results,affy.snps.results, joint.snps.results, 
       file=cat.path(home.dir,"finalResults5NEW",suf=dd,ext="RData"))
  loop.tracker(dd,22)
}


##### COMBINE RESULTS INTO 1 TABLE ####

nchr <- 22

results1 <- cat.path(paste0(home.dir,"finalResults"),"finalResults5NEW",suf=1:nchr,ext="RData")

All.illu <- All.affy <- All.joint <- vector("list",nchr)

for (dd in 1:nchr) { 
  rr <- reader(results1[dd]) 
  All.illu[[dd]] <- rr$illu.snps.results
  All.affy[[dd]] <- rr$affy.snps.results
  All.joint[[dd]] <- estimates.to.results(rr$estz5)
  loop.tracker(dd,nchr)
}

All.illu.res <- as.data.frame(do.call("rbind",args=All.illu)[,c(1,3,5,4)])
All.affy.res <- as.data.frame(do.call("rbind",args=All.affy)[,c(1,3,5,4)])
All.joint.res <- do.call("rbind",args=All.joint)[,c(1,2,4,3)]

All.illu.res[,4] <- "illu"
All.affy.res[,4] <- "affy"
All.joint.res[,4] <- "joint"

colnames(All.illu.res) <- colnames(All.affy.res) <-  colnames(All.joint.res) <-  c("OR","SE","p","source")

all.res <- rbind(All.illu.res, All.affy.res, All.joint.res)


## ANNOTATE RESULTS WITH DIRECTIONS, MAJOR/MINOR ALLELEs, ETC ##
if(F) {
  # dir should be ~/barrett/
  X <- reader(SML[[22]])
  XX <- X[,1:1000]
  #res <- do.call("rbind",args=ii)
  #res <- do.one(XX)
  ii <- match(rownames(samples(X)),rownames(sample.info))
  #rn <- rownamesL(SML,list=F)
  
  ph <- sample.info$phenotype[ii]-1 # need this pheno before running 'caseway', used as global
  if(F) {
    res <- smlapply(SML,do.one,list,n.cores=3,sample.info=sample.info) #real one
    result <- do.call("rbind",args=res)
    cn <- colnamesL(SML,list=F)
    rownames(result) <- cn
  } else {
    cn <- colnamesL(SML,list=F)
    result <- all.res ## add to existing result NOT generated using 'do.one' ##
  }
  case.way <- smlapply(rev(SML),caseway,list,n.cores=1,pheno=ph) 
  case.way <- rev(case.way[[1]])
  maj.min <- smlapply(rev(SML),majmin,list,n.cores=1,pheno=ph) 
  maj.min <- rev(maj.min[[1]])
  A1 <- smlapply(SML,function(x) { (snps(x)$allele.1)},list,n.cores=3)
  A2 <- smlapply(SML,function(x) { (snps(x)$allele.2)},list,n.cores=3)
  if(T){
    ## add to existing result NOT generated using 'do.one' ##
    anno <- cbind(paste(unlist(maj.min)),paste(unlist(case.way)),paste(unlist(A1[[1]])),paste(unlist(A2[[1]])))
    rownames(anno) <- cn
    result[["majmin"]] <- anno[match(rownames(result),rownames(anno)),1]
    result[["caseway"]] <- anno[match(rownames(result),rownames(anno)),2]
    result[["A1"]] <- anno[match(rownames(result),rownames(anno)),3]
    result[["A2"]] <- anno[match(rownames(result),rownames(anno)),4]
    result[["interp"]] <- rep("???",nrow(result))
    result[["interp"]][with(result,caseway=="CasesRef+" & majmin=="major")] <- "T1D Major+"
    result[["interp"]][with(result,caseway=="CasesRef+" & majmin=="minor")] <- "T1D Minor+"
    result[["interp"]][with(result,caseway=="CasesRef-" & majmin=="minor")] <- "T1D Major+"
    result[["interp"]][with(result,caseway=="CasesRef-" & majmin=="major")] <- "T1D Minor+"
    result[["OR_orig"]] <- result$OR
    result[["OR"]][result$OR>1 & result$interp=="T1D Major+"] <- 1/(result[["OR"]][result$OR>1 & result$interp=="T1D Major+"])
    result[["OR"]][result$OR<1 & result$interp=="T1D Minor+"] <- 1/(result[["OR"]][result$OR<1 & result$interp=="T1D Minor+"])
  } else {
    result[["majmin"]] <- unlist(maj.min)
    result[["case.way"]] <- unlist(case.way)
    
    result[["alt"]] <- unlist(A1[[1]])
    result[["ref"]] <- unlist(A2[[1]])
  }
  save(result,file="/chiswick/data/ncooper/barrett/gwasResultsAnnotJUN4.RData")
  save(res,file="/chiswick/data/ncooper/barrett/gwasResultsProperJUN4.RData")
} else {
  load("/chiswick/data/ncooper/barrett/gwasResultsAnnotJUN4.RData")
}


