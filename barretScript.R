library(annotSnpStats)
library(humarray)
library(gtools)
library(GenomicRanges)
library(reader)
library(snpStatsWriter)
source("~/github/plumbCNV/SnpMatrixList.R")

dat.dir <- ("/dunwich/scratch/chrisw/T1DGC/snpStats-genotypes")
home.dir <- ("/chiswick/data/ncooper/barrett/")
setwd(home.dir)

cohorts <- c("C58","WTCCC","T1DGC")
sml <- vector("list",length(cohorts)); names(sml) <- cohorts
ll <- list.files(dat.dir,pattern="C58"); oo <- mixedorder(gsub("C58.","",ll)); sml[[1]] <- as.list(cat.path(dat.dir,ll[oo]))
ll <- list.files(dat.dir,pattern="WTCCC"); oo <- mixedorder(gsub("WTCCC.","",ll)); sml[[2]] <- as.list(cat.path(dat.dir,ll[oo]))
ll <- list.files(dat.dir,pattern="T1DGC"); oo <- mixedorder(gsub("T1DGC.","",ll)); sml[[3]] <- as.list(cat.path(dat.dir,ll[oo]))
#ll <- list.files(dat.dir,pattern="US"); oo <- mixedorder(gsub("US.","",ll)); sml4 <- as.list(ll[oo]) # not allowed to use

#lcs  <- list.colsummary(sml[[1]])
#summary(lcs$z.HWE)

## files we need
ndir <- "/chiswick/data/ncooper/imputation/"
#f.geno <- paste0(ndir,diseases, if(is.null(args$arm)) {""} else {"/PQDATA"},"/CHR",args$chr,args$arm,".RData")
f.legend <- mixedsort(list.files(cat.path(home.dir,"1000GP_Phase3"),pattern="legend",full.names=TRUE))
#f.legend <- paste0(ndir,"COMMON/ALL.chr",args$chr,".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz")

init <- FALSE # TRUE to recalculate support files (slow)

if(init) {
  # if we need to regenerate the snp and sample support files.. #
  #### PROCESS SNP SUPPORT #####
  # Part 1 of SNPinfo for C58 and WTCCC
  snpsup <- list.files("/dunwich/scratch/chrisw/T1DGC/R-objects/corrected-support",pattern="support",full.names=TRUE)
  snpsup <- snpsup[grep("unplaced",snpsup,invert=TRUE)]
  snpsup <- mixedsort(snpsup)
  snpsup <- lapply(snpsup,reader)
  snpsup <- do.call("rbind",args=snpsup)
  snpsup[["uid"]] <- paste(snpsup$b36_chr,snpsup$b36_start,sep="_")
  anydups <- which(duplicated(snpsup[,"uid"]))
  if(length(anydups)>0) {
    snpsup <- snpsup[-anydups,]; cat("removed",length(anydups),"duplicates\n")
  }
  snpsup <- snpsup[,-1] # remove build 35 column
  colnames(snpsup) <- gsub("b36_","",colnames(snpsup))
  colnames(snpsup)[1] <- "rs.id"
  snpsup[["allele.1"]] <- sapply(strsplit(snpsup$wtccc.affy_alleles,"/"),"[",1)
  snpsup[["allele.2"]] <- sapply(strsplit(snpsup$wtccc.affy_alleles,"/"),"[",2)
  snp.info2 <- df.to.ranged(snpsup,end="stop")
  snp.info2 <- toGenomeOrder(snp.info2)
  snpsup[["allele.1"]] <- sapply(strsplit(snpsup$ilmn_alleles,"/"),"[",1)
  snpsup[["allele.2"]] <- sapply(strsplit(snpsup$ilmn_alleles,"/"),"[",2)
  snp.info <- df.to.ranged(snpsup,end="stop")
  snp.info <- toGenomeOrder(snp.info)
  ### convert from 36 - 37 ######
  snp.info <- conv.36.37(snp.info)
  snp.info2 <- conv.36.37(snp.info2)
  #### END PROCESS SNP SUPPORT #####
  
  #### PROCESS SAMPLE SUPPORT #####
  samp.sup1 <- reader("/ipswich/data/T1DGC/R-objects/support/T1DGC-sample-support.RData") # 100% of WTCCC ids in col 5 + 96% of T1DGC ids
  samp.sup2 <- reader("/ipswich/data/T1DGC/R-objects/support/WTCCC-sample-support.RData") # 82% of C58 ids in col 2 
  rownames(samp.sup1) <- samp.sup1[,2]
  rownames(samp.sup2) <- samp.sup2[,5]
  smp <- rownamesL(sml[[1]])
  sso <- samp.sup2[1:47,]
  rownames(sso) <- smp[(which(!smp %in% rownames(samp.sup2)))]
  sso[,2:8] <- NA
  samp.sup2 <- rbind(samp.sup2,sso)
  colnames(samp.sup2)[c(5,6,7)] <- c("sample","b58cregion","subject") 
  samp.sup2 <- samp.sup2[,c(5,7,2,3,6,8,4)] # rearrange cols to match other sample support v
  smp <- rownamesL(sml[[3]])
  sso <- samp.sup1[1:1422,]
  rownames(sso) <- smp[(which(!smp %in% rownames(samp.sup1)))]
  sso[,2:8] <- NA
  samp.sup1 <- rbind(samp.sup1,sso)
  samp.sup1 <- samp.sup1[,c(2,1,3,4,8,6,9)] # rearrange cols to match other sample support ^
  sample.info <- rbind(samp.sup1,samp.sup2)
  colnames(sample.info)[3] <- "phenotype"
  #### END PROCESS SAMPLE SUPPORT #####
  save(snp.info,snp.info2,sample.info,file="/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")
} else {
  # otherwise load all support from this file #
  (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData"))
}


chrz <- 1:22 # autosomes

for (cc in 3:7) { #chrz) {
  # process impute file for each chromosome in turn
  source("~/github/imputer/barretChr.R")
}


P1 <- 1; P2 <- 3*10^6; ch <- 1

map <- paste0("/chiswick/data/ncooper/barrett/1000GP_Phase3/genetic_map_chr",ch,"_combined_b37.txt")
legend <- paste0("/chiswick/data/ncooper/barrett/1000GP_Phase3/1000GP_Phase3_chr",ch,".legend.gz")
haps <- paste0("/chiswick/data/ncooper/barrett/1000GP_Phase3/1000GP_Phase3_chr",ch,".hap.gz")
f.in <- paste0("/chiswick/data/ncooper/barrett/IMPUTE_INPUT/impute-",ch)
f.out <- paste0("/chiswick/data/ncooper/imputation/T1D/impute-",ch,"-",P1,"-",P2,".out")
imp <- "/home/chrisw/local/bin/impute2"
com <- paste0(imp," -m ",map," -h ",haps," -l ",legend," -g ",f.in," -o ",f.out," -int ",P1," ",P2)


