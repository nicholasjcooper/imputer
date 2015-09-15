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
setwd(home.dir)

cohorts <- c("C58","WTCCC","T1DGC")
sml <- vector("list",length(cohorts)); names(sml) <- cohorts
ll <- list.files(dat.dir,pattern="C58"); oo <- mixedorder(gsub("C58.","",ll)); sml[[1]] <- as.list(cat.path(dat.dir,ll[oo])) #snp.info2 affy
ll <- list.files(dat.dir,pattern="WTCCC"); oo <- mixedorder(gsub("WTCCC.","",ll)); sml[[2]] <- as.list(cat.path(dat.dir,ll[oo])) #snp.info2 affy
ll <- list.files(dat.dir,pattern="T1DGC"); oo <- mixedorder(gsub("T1DGC.","",ll)); sml[[3]] <- as.list(cat.path(dat.dir,ll[oo])) #snp.info  illumina
ll <- list.files(raw.dir,pattern="WT.C58.snp"); oo <- mixedorder(gsub("WT.C58.snp.","",ll)); smlB <- as.list(cat.path(raw.dir,ll[oo]))[1:22] # affy

#ll <- list.files(dat.dir,pattern="US"); oo <- mixedorder(gsub("US.","",ll)); sml4 <- as.list(ll[oo]) # not allowed to use


#lcs  <- list.colsummary(sml[[1]])
#summary(lcs$z.HWE)

## files we need
#f.geno <- paste0(ndir,diseases, if(is.null(args$arm)) {""} else {"/PQDATA"},"/CHR",args$chr,args$arm,".RData")
f.legend <- mixedsort(list.files(cat.path(home.dir,"1000GP_Phase3"),pattern="legend",full.names=TRUE))
#f.legend <- paste0(ndir,"COMMON/ALL.chr",args$chr,".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz")

init <- FALSE # TRUE to recalculate support files (slow)
#proc.impute.files <- FALSE # whether to regenerate the raw impute files (slow)
init.cmds <- FALSE
write.imputes <- FALSE
write.snpmatrix <- TRUE
if(write.imputes | write.snpmatrix) { proc.impute.files <- TRUE } else { proc.impute.files <- FALSE }
base.dir <- "/home/ncooper/barrett/"  #"/chiswick/data/ncooper/barrett/"
ndir <- cat.path(base.dir,"IMPUTE_INPUT")


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
  snp.info <- conv.36.37(snp.info)  # illumina
  snp.info2 <- conv.36.37(snp.info2)  # affy
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
  ##### FIX PHENOS WITH INFILE SUPPORT #####
  (load("/dunwich/scratch/chrisw/T1DGC/snpStats-genotypes/C58.22.cleaned-genotypes.RData"))
  ll <- list.files(dat.dir,pattern="T1DGC",full.names=T)
  (load(ll[1]))
  ll <- list.files(dat.dir,pattern="WTCCC",full.names=T)
  (load(ll[1]))
  wtcn <- names(WTCCC.cohort)
  wtcc <- as.character(WTCCC.cohort)
  wtcc[wtcc!="T1D"] <- 0
  wtcc[wtcc=="T1D"] <- 1
  tgcn <- names(T1DGC.cohort)
  tgcc <- as.character(T1DGC.cohort)
  tgcc[tgcc!="case"] <- 0
  tgcc[tgcc=="case"] <- 1
  c58n <- rownames(C58.22.data)
  c58c <- rep(0,length(c58n))
  all.n <- c(c58n,wtcn,tgcn)
  all.c <- c(c58c,wtcc,tgcc)
  sample.info$phenotype[is.na(sample.info$phenotype)] <- 1
  sample.info[which(sample.info$phenotype==1 & all.c[ii]==1),"phenotype"] <- NA # mismatch between files
  #### END PROCESS SAMPLE SUPPORT #####
  save(snp.info,snp.info2,sample.info,file="/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")
} else {
  # otherwise load all support from this file #
  cat("loading support files ... ")
  if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) }
  cat("done\n")
}


chrz <- 1:22 # autosomes

if(proc.impute.files) {
  for (cc in rev(1:22)) { #chrz) {
    # process impute file for each chromosome in turn
    source("~/github/imputer/barretChr.R")
  }
}



if(init.cmds) {
  snp.info3 <- rbind(snp.info,snp.info2)
  snp.info3 <- toGenomeOrder(snp.info3)
  ci <- chrInfo(snp.info3) # get starts and ends of each chromosome
  ii <- vector("list",length(chrz))
  ci <- ci[,1:2]-ci[,3]
  for (jj in chrz) {
    ii[[jj]] <- seq(from = ci[jj,1], to=ci[jj,2], by=1*10^6)
  }
  cat("Generating IMPUTE2 commands for whole genome\n")
  com <- character() # will have one line for each IMPUTE2 command
  ff <- 0
  for (jj in chrz) {
    for (kk in 1:(length(ii[[jj]])-1)) {
      ff <- ff+1
      P1 <- ii[[jj]][kk]; P2 <- ii[[jj]][kk+1]; ch <- jj
      map <- paste0(base.dir,"1000GP_Phase3/genetic_map_chr",ch,"_combined_b37.txt")
      legend <- paste0(base.dir,"1000GP_Phase3/1000GP_Phase3_chr",ch,".legend.gz")
      haps <- paste0(base.dir,"1000GP_Phase3/1000GP_Phase3_chr",ch,".hap.gz")
      f.in <- paste0(base.dir,"IMPUTE_INPUT/impute-",ch)
      f.out <- paste0(base.dir,"OUTPUT/impute-",ch,"-",P1,"-",P2,".out")
      imp <- "/home/chrisw/local/bin/impute2"
      filt.rules <- " -filt_rules_l 'EUR<0.01'"
      com[ff] <- paste0(imp," -m ",map," -h ",haps," -l ",legend," -g ",f.in," -o ",f.out," -int ",P1," ",P2,filt.rules)
      writeLines(com[ff],con=cat.path(home.dir,pref="IMPUTE_CMDS/",fn="imp",suf=ff,ext="sh"))
    }
    loop.tracker(jj,length(chrz))
  }
  
  cat(ff,"files/commands produced\n")
}

#writeLines(paste(com),con=paste0(base.dir,"commands2870.txt"))
com <- readLines(paste0(base.dir,"commands2870.txt"))
com <- paste(com,"-pgs -pgs_miss")
#

source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
#bash.qsub(com[1:1000],dir=paste0(base.dir,"IMPUTE_OUTPUT"),grid.name="bigmem",logpref="BIG1_",stagger=30,max.conc=100) # in prog
#bash.qsub(com[1001:1500],dir=paste0(base.dir,"IMPUTE_OUTPUT2"),grid.name="fourcpu",logpref="FOUR1_",stagger=30,max.conc=50) # begun
#bash.qsub(com[1501:2000],dir=paste0(base.dir,"IMPUTE_OUTPUT2"),grid.name="fourcpu",logpref="FOUR2_",stagger=30,max.conc=50)
#bash.qsub(com[2001:2870],dir=paste0(base.dir,"IMPUTE_OUTPUT"),grid.name="bigmem",logpref="BIG2_",stagger=30,max.conc=100)
#use grid.name="fourcpu -l h_data=6500000000" to set a ram limit..


#for i in {4278282..4278306}
#do
#qdel "$i"
#done

#bash.qsub(com[1041:1500], dir=paste0(base.dir,"IMPUTE_OUTPUT4"),grid.name="bigmem",logpref="BIGF_",stagger=180,max.conc=50)
#bash.qsub(com[1141:1500], dir=paste0(base.dir,"IMPUTE_OUTPUT4"),grid.name="bigmem",logpref="BIGG_",stagger=60,max.conc=50)

#bash.qsub(com[1500:2000], dir=paste0(base.dir,"IMPUTE_OUTPUT4"),grid.name="bigmem",logpref="BIGI_",stagger=60,max.conc=50)

#bash.qsub(com[2000:2500], dir=paste0(base.dir,"IMPUTE_OUTPUT7"),grid.name="bigmem",logpref="BIGK_",stagger=60,max.conc=50)


#### NEW LOT ####
#bash.qsub(com[2490:2870][-1:-72], dir=paste0(base.dir,"IMPUTE_OUTPUT"),grid.name="bigmem",logpref="BIGA2_",stagger=60,max.conc=50)
#bash.qsub(com[2000:2490], dir=paste0(base.dir,"IMPUTE_OUTPUT2"),grid.name="bigmem",logpref="BIGB_",stagger=60,max.conc=50)

#bash.qsub(com[1500:2000], dir=paste0(base.dir,"IMPUTE_OUTPUT3"),grid.name="bigmem",logpref="BIGC_",stagger=60,max.conc=50)
#bash.qsub(com[1000:1500], dir=paste0(base.dir,"IMPUTE_OUTPUT4"),grid.name="bigmem",logpref="BIGD_",stagger=600,max.conc=50)

#bash.qsub(com[c(1000:1500)[c(349:501,164,195,198,199,335,336,337,338,340,343,345,346)]],
#  dir=paste0(base.dir,"IMPUTE_OUTPUT4b"),grid.name="bigmem",logpref="BIGDb_",stagger=60,max.conc=50)
#bash.qsub(com[c(1500:2000)[c(411:501,204,213,241,251,267)]],
#  dir=paste0(base.dir,"IMPUTE_OUTPUT3b"),grid.name="bigmem",logpref="BIGCb_",stagger=60,max.conc=50)

## in progress
#bash.qsub(com[500:1000], dir=paste0(base.dir,"IMPUTE_OUTPUT5"),grid.name="bigmem",logpref="BIGE_",stagger=60,max.conc=50)


#bash.qsub(com[1:500], dir=paste0(base.dir,"IMPUTE_OUTPUT6"),grid.name="bigmem",logpref="BIGF_",stagger=180,max.conc=50)

#bash.qsub(com[c(456,470,472,474,475,481)], dir=paste0(base.dir,"IMPUTE_OUTPUT6"),grid.name="bigmem",logpref="BIGZ_",stagger=10,max.conc=50)
if(F) {
# those on chr7 that failed
ccz <- c(1282,1284,1297,1300,1301,1307,1308,1309,1310,1311,1312,1313,1315,1316,1318,1321,1322,1323,1324,1325,1326,1327,1328,1330,1331,1332,1333,1338,1340,1341,1343,1346,1347)
bash.qsub(com[ccz],
  dir=paste0(base.dir,"IMPUTE_OUTPUT4b"),grid.name="bigmem",logpref="BIGY_",stagger=60,max.conc=50)


ddz <- c(1845,1853,1854,1855,1859,1861,1866,1867,1868,1869,1870,1871,1873,1878,1879,1880,1881,1883,1884,1885,1886,1887,1888,1889,1891,1892,1893,1894,1895,1896,1897,1898,1899,1900,1901,1902,1903,1904,1905,1906,1907,1908,1909)
bash.qsub(com[ddz],
          dir=paste0(base.dir,"IMPUTE_OUTPUT3b"),grid.name="bigmem",logpref="BIGX_",stagger=60,max.conc=50)


bash.qsub(com[ww],dir=paste0(base.dir,"IMPUTE_OUTPUTX"),grid.name="bigmem",
          logpref="BIGW_",stagger=60,max.conc=96)
}