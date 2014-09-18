
t1.dr <- "/chiswick/data/ncooper/imputation/T1D/" 
gr.dr <- "/chiswick/data/ncooper/imputation/GRAVES/"
ms.dr <- "/chiswick/data/ncooper/imputation/MS/"
cl.dr <- "/chiswick/data/ncooper/imputation/COELIAC/"
ja.dr <- "/chiswick/data/ncooper/imputation/JIA/"
ra.dr <- "/chiswick/data/ncooper/imputation/RA/"

#dir1 <- "/chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK"
#dir1 <- "/ipswich/data/Immunochip/PLINK/distribution/" # coeliac, graves

dir1 <- "/ipswich/data/chrisw/MS-sawcer/MS_ICdata"  # MS - check binding between cases and specific controls
#dir1 <- "/chiswick/data/store/immunochip/JIA/plink_dataset_files/"
#dir1 <- "/ipswich/data/chrisw/IChip/RA-control-data"   # RA - international cases and controls in same file


my.dr <- gr.dr

setwd(my.dr)

source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
library(NCmisc)
library(reader)
library(snpStats)
library(devtools)
library(annotSnpStats)
#install_github("annotSnpStats",username="chr1swallace")

#fn <- "finn-preqc"
#fn <- "graves-preqc"
fn <- character()
fn[1] <- "Cases_UKIC"
fn[2] <- "Controls_UKIC"


if(F) {
  fn <- character(6)
  ###RA ONES###
#
fn[1] <- "ES_RA_immunochip_dataset/iChip_RACI_PhaseII_ES_QCgp"
#:
fn[2] <- "NL_RA_immunochip_dataset/iChip_RACI_PhaseII_NL_QCgp"
#:
fn[3] <- "SE-E_RA_immunochip_dataset/iChip_RACI_PhaseII_SE-E_QCgp"
#:
fn[4] <- "SE-U_RA_immunochip_dataset/iChip_RACI_PhaseII_SE-U_QCgp"
#:
fn[5] <- "UK_RA_immunochip_dataset/iChip_RACI_PhaseII_UK_QCgp"
#:
fn[6] <- "US_RA_immunochip_dataset/iChip_RACI_PhaseII_US_QCgp"

#fn <- "Immunochip_JIA_PhaseIInew_QCgp_All_SNPQC_UKonly" #jia pt1
#fn <- "Immunochip_JIA_PhaseIInew_QCgp_All_SNPQC_UKonly_polygo" #a subset of jia... not used ?
}


print(load(cat.path(getwd(),"raAnnotSnpStats",ext="RData")))


## combine by groups and split by chromosome
chrz <- rev(c(1:23,25))
for (cc in 1:length(chrz)) {
  nsets <- length(myMat); each.chr <- vector("list",nsets)
  for (nn in 1:nsets) {
    each.chr[[nn]] <- get.chr(myMat[[nn]],chr=chrz[cc])
    loop.tracker((nsets*(cc-1))+nn,nsets*24)
  }
  #ra <- do.call("rbind2",args=each.chr)
  save(each.chr,file=cat.path(ra.dr,"ra.chr",suf=chrz[cc],ext="RData"))
}


if(F) {
  myMat <- vector("list",6)
  for (cc in 1:6) {
    bed <- paste0(dir1,"/",fn[cc],".bed"); bim <- paste0(dir1,"/",fn[cc],".bim"); fam <- paste0(dir1,"/",fn[cc],".fam")
    myMat[[cc]] <- annot.plink(read.plink(bed, bim, fam))  
  }

  myMat <- sync.asnp.mats(myMat)

  names(myMat) <- fn

  save(myMat,file=cat.path(getwd(),"raAnnotSnpStats",ext="RData"))

  jia <- annot.plink(myMat)
  save(jia,file=cat.path(getwd(),"jiaAnnotSnpStats",ext="RData"))
  
  
  graves <- annot.plink(myMat)
  save(graves,file=cat.path(getwd(),"gravesAnnotSnpStats",ext="RData"))

  coeliac <- annot.plink(myMat)
  save(coeliac,file=cat.path(getwd(),"coeliacAnnotSnpStats",ext="RData"))
}

if(F) {
  
  head(myMat$genotypes)
  tail(myMat$genotypes)
  
  sMat <- myMat$genotypes
  map <- myMat$map
  snp.info <- make.snp.info(map=myMat$map)
  sample.info <- make.sample.info(dir=getwd(),id.list=rownames(myMat$fam),
                                  phenotype=column.salvage(myMat$fam,"phenotype","affected"))
  sel.22 <- which(colnames(sMat) %in% rownames(select.autosomes(snp.info)))
  rs <- row.summary(sMat[,sel.22])
  cs <- col.summary(sMat)
  
  #save(cs,rs,file="summaries.finns.RData")
  #save(cs,rs,cs2,file="summaries.finns.RData")
  #save(cs,rs,cs2,snp.info,sample.info,file="summaries.finns.RData")
  
  # assume chr ["X", "Y", "XY", "MT"] == [23, 24, 25, 26]
  #load("summaries.finns.RData")
  
  head(cs)
  head(rs)
  
  call.rate.summary(cs)
  call.rate.summary(rs)
  
  het.density.plots(rs,het.lo=.185,het.hi=.235)
  excl.samp.sel <- hz.vs.callrate.plots(rs,callrate.samp.thr=.97,het.lo=.185,het.hi=.235,excl=T) 
  excl.samp <- rownames(rs)[excl.samp.sel]
  cs2 <- col.summary(sMat[-excl.samp.sel,])
  writeLines(excl.samp,"callrateHzFailingSamples.txt")
  
  
  system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-duplicate-samples.txt ~/Finnish/")
  system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-monomorph-snps.txt  ~/Finnish/")
  system("cp /chiswick/data/Immunochip/project-added-by-Neil-2013-08-27/PLINK/finn-postqc-sexcheck-samples.txt  ~/Finnish/")
  
  mendelz <- read.table(cat.path(dir1,"finn-d6-mendel.imendel"),header=T,comment="")
  ME <- mendelz$N
  hist(ME)
  plot(density(ME[ME<200]))
  
  length(which(ME>200))/length(ME)
  
  hwe.density.plots(cs2)
  hwe.vs.callrate.plots(cs2)
  
  
  sample.info[["call.rate"]] <- rs$Call.rate
  snp.info[["call.rate"]] <- cs$Call.rate
  draw.density.plots(cs,snp.info=snp.info,sample.info=sample.info)
  
}
