

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

#"impute-11-45000001-46000001.out"
#"impute-7-99000001-100000001.out"



### run gwas on impute files for each Chr in turn ##
for (dd in c(5:1)) {
  cat("Chr",dd," in progress ..")
  fnz <- get.impute.chr.fn("OUTPUT",dd)
  bad <- bad2 <- c("notarealid")
  for (ee in 1:length(bad)) {
    if(length(grep(bad[ee],fnz))>0) { fnz <- fnz[-grep(bad[ee],fnz)] }
  }
  for (ee in 1:length(bad2)) {
    if(length(grep(bad2[ee],fnz))>0) { fnz <- fnz[-grep(bad2[ee],fnz)] }
  }
  cat(".")
  res.list <- lapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE)
  cat(".")
  SNPMAT <- do.call("rbind",args=res.list)
  cat(".")
  save(SNPMAT,res.list,file=cat.path("/chiswick/data/ncooper/barrett/",pref="chr",fn=dd,suf="_GWASresults",ext=".RData"))
  cat(" done\n")
}


bash.qsub(com[c(456,470,472,474,475,481)], dir=paste0(base.dir,"IMPUTE_OUTPUT6"),grid.name="bigmem",logpref="BIGZ_",stagger=10,max.conc=50)

com2 <- NULL;
for (ee in 1:length(bad2)) {
  if(length(grep(bad2[ee],com))>0) { com2 <- c(com2,com[grep(bad2[ee],com)]) }
}

com3 <- NULL;
for (ee in 1:length(bad)) {
  if(length(grep(bad[ee],com))>0) { com3 <- c(com3,com[grep(bad[ee],com)]) }
}

which(com %in% com2)
which(com %in% com2)

list.files(cat.path(base.dir,"IMPUTE_OUTPUT6/check_flags"))

iii <- paste0("ls ~/barrett/OUTPUT/*",bad2[-1],"*")
for (jjj in 1:length(iii)) {
  system(iii[jjj])
}

iii <- paste0("ls ~/barrett/OUTPUT/*",bad,"*")
for (jjj in 1:length(iii)) {
  system(iii[jjj])
}

##
fnz <- get.impute.chr.fn("OUTPUT",22)

fn <- fnz[1]
rez <- analyse.impute.file(fn,sample.info)





fnzC[(which(!fnzC %in% fnz))]






# test which segments of 1MB are missing throughout the genome

mis.list <- vector("list",22)
for (dd in 1:22) {
  fnz <- gsub(".out","",gsub(paste0("/home/ncooper/barrett/OUTPUT/impute-",dd,"-"),"",get.impute.chr.fn("~/barrett/OUTPUT",dd),fixed=T),fixed=T)
  enz <- 1+(10^6*(1:(get.chr.lens()[dd] %/% 10^6)))
  stz <- c(1,enz[-length(enz)])
  fnzC <- paste0(stz,"-",enz)
  cat("missing:",length(which(!fnzC %in% fnz))," ; present: ",length(which(fnzC %in% fnz)),"\n")
  mis.list[[dd]] <- paste0("impute-",dd,"-",fnzC[(which(!fnzC %in% fnz))])
}



missos <- unlist(mis.list)

ww <- NULL; for (ee in 1:length(missos)) { ww <- c(ww,grep(missos[ee],com)) }



bad <- paste0("impute-11-",
  c("34000001-35000001","42000001-43000001","43000001-44000001","44000001-45000001",
"45000001-46000001",
"48000001-49000001",
"50000001-51000001","55000001-56000001","56000001-57000001","57000001-58000001","58000001-59000001","59000001-60000001",
"60000001-61000001","62000001-63000001","67000001-68000001","68000001-69000001","69000001-70000001",
"70000001-71000001","72000001-73000001","73000001-74000001","74000001-75000001","75000001-76000001",
"76000001-77000001","77000001-78000001","78000001-79000001","80000001-81000001","81000001-82000001",
"82000001-83000001","83000001-84000001","84000001-85000001","85000001-86000001","86000001-87000001",
"87000001-88000001","88000001-89000001","89000001-90000001","90000001-91000001","91000001-92000001",
"92000001-93000001","93000001-94000001","94000001-95000001","95000001-96000001","96000001-97000001",
"97000001-98000001","98000001-99000001")
  ,".out")

  bad2 <- paste0("impute-7-",c("100000001-101000001","101000001-102000001","102000001-103000001","103000001-104000001","108000001-109000001",
"110000001-111000001","111000001-112000001","113000001-114000001","116000001-117000001","117000001-118000001",
"52000001-53000001","54000001-55000001","67000001-68000001","70000001-71000001","71000001-72000001",
"77000001-78000001","78000001-79000001","79000001-80000001","80000001-81000001","81000001-82000001",
"82000001-83000001","83000001-84000001","85000001-86000001","86000001-87000001","88000001-89000001",
"91000001-92000001","92000001-93000001","93000001-94000001","94000001-95000001","95000001-96000001",
"96000001-97000001","97000001-98000001","98000001-99000001",
"99000001-100000001")   ,".out")


# read in and save impute2 information files (rather than the GWAS files)

for (dd in c(22:1)) {
  cat("Chr",dd," in progress ..")
  fnz <- get.impute.chr.fn("OUTPUT",dd,info.only=TRUE,by.sample=FALSE)

  res.list <- lapply(as.list(fnz),reader)
  SNPINFO <- do.call("rbind",args=res.list)
  save(SNPINFO,res.list,file=cat.path("/chiswick/data/ncooper/barrett/",pref="chr",fn=dd,suf="_GWASinfo",ext=".RData"))
  cat(" done\n")
}

smp <- rownamesL(SML)


for (dd in c(22:1)) {
  cat("Chr",dd," in progress ..")
  fnz <- get.impute.chr.fn("OUTPUT",dd,info.only=TRUE,by.sample=TRUE)
  res.list <- lapply(as.list(fnz),reader)
  res.list <- lapply(res.list,function(x) { x[["sampleid"]] <- smp; return(x) })
  lenz <- sapply(res.list,nrow)
  rngz <- gsub("OUTPUT/impute-","",sapply(strsplit(fnz,".",fixed=TRUE),"[",1))
  lab.col <- rep(rngz,lenz)
  SAMPINFO <- do.call("rbind",args=res.list)
  SAMPINFO[["range"]] <- lab.col
  save(SAMPINFO,res.list,file=cat.path("/chiswick/data/ncooper/barrett/",pref="chr",fn=dd,suf="_GWASinfo_by_sample",ext=".RData"))
  cat(" done\n")
}


(load("chr22_GWASinfo_by_sample.RData"))
(load("chr22_GWASinfo.RData"))
(load("chr22_GWASresults.RData"))


save(SML,SML2,outl,outl2,file="/chiswick/data/ncooper/barrett/nickWorkingThu.RData")

### run for PCs derived for intersecting SNPs in barrett
setwd("~/barrett")
lambda.1000.pc4 <- numeric(22)
res.list.pc4 <- vector("list",22)
for (j in 21:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- lapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,nPCs=4)
  res.list.pc4[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.pc4[[j]]$p)^2,na.rm=T)/.454,3)
  lambda.1000 <- lambda_nm(lambda.nxt,nr=9110,mr=6158)
  lambda.1000.pc4[j] <-  lambda.1000
  save(res.list.pc4,lambda.1000.pc4,file="/chiswick/data/ncooper/barrett/L1000Calcs_pc4.RData")
  cat("L1000 for chr",j,"=",lambda.1000,"\n")
}



## get lamba1000s for PCAs 1-4 GWAS analysis ##
(load("/chiswick/data/ncooper/barrett/L1000Calcs_pc4.RData"))
lambda.1000.pc4

## get lamba1000s for initial naive GWAS analysis ##
for (j in 1:22) {
  (load(cat.path("/chiswick/data/ncooper/barrett/",pref="chr",fn=j,suf="_GWASresults",ext=".RData")))
  cat("L1000 for chr",j,"=",lambda_nm( (median(p.to.Z(SNPMAT$p),na.rm=T)^2)/.454,nr=9110,mr=6158),"\n")
}




lambdas.naive <- c(2.167,2.047,2.041,2.141,2.133,2.312,2.196,1.869,1.979,2.219,1.851,1.939,2.135,2.199,2.202,2.591,2.614,2.513,2.736,2.935,2.319,2.638)

lambdas.pc4  <-  c(2.150,2.033,2.027,2.122,2.115,2.290,2.174,1.853,1.968,2.208,1.836,1.923,2.126,2.181,2.187,2.567,2.589,2.493,2.699,2.903,2.303,2.610)

### run for PCs derived for intersecting SNPs in 1000 genomes, loadings applied to aligned barrett

if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) }

setwd("~/barrett")
lambda.1000.pc4g <- numeric(22)
res.list.pc4g <- vector("list",22)
for (j in 2:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- lapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,nPCs=4,PC.source.fn="/chiswick/data/ncooper/barrett/PCsFor1000gPCAIntersectSnpsForT1Dand1KG.RData")
  res.list.pc4g[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.pc4g[[j]]$p)^2,na.rm=T)/.454,3)
  lambda.1000 <- lambda_nm(lambda.nxt,nr=9110,mr=6158)
  lambda.1000.pc4g[j] <-  lambda.1000
  save(res.list.pc4g,lambda.1000.pc4g,file="/chiswick/data/ncooper/barrett/L1000Calcs_pc4g.RData")
  cat("L1000 for chr",j,"=",lambda.1000,"\n")
}


## COMPARE ORIGINAL VS IMPUTED DATA
fj <- cat.path(paste0("~/barrett/","IMPUTE_INPUT"),fn="snpmat",suf=22,ext="RData") # SnpMatrix of data just before 'write.impute'
FJ <- reader(fj) # read it in
jf <- reader("~/barrett/OUTPUT/impute-22-24000001-25000001SnpMatrix.RData") # SnpMatrix of imputed data
cjf <- colnames(jf) 
cjf <- sapply(strsplit(cjf,":",fixed=T),"[",1)  # extract just the rs id from the imputed colnames
colnames(jf) <- cjf
FFF <- (FJ[,colnames(jf[,colnames(jf) %in% colnames(FJ)])]) # get just the SNPs that are in the original dataset and in imputed dataset
GGG <- (jf[,colnames(jf) %in% colnames(FJ)])   # get just the SNPs that are in the original dataset
cs1 <- col.summary(FFF)
cs2 <- col.summary(GGG)
HHH <- convert.snpmat.uncertain.to.normal(GGG) # convert uncertain genotypes to regular
fff <- FFF[aff.only.ids,]  # get affy barrett samples
hhh <- HHH[aff.only.ids,]  # get matching imputed data
CS1 <- col.summary(fff)
hhh <- as(hhh,"SnpMatrix")
CS2 <- col.summary(hhh)
hhh2 <- HHH[ill.only.ids,]  # get illumina original samples
hhh2 <- as(hhh2,"SnpMatrix")
fff2 <- FFF[ill.only.ids,]  # get matching imputed data
CS12 <- col.summary(fff2)
CS22 <- col.summary(hhh2)
hhh <- hhh[,CS1$Calls>5000,drop=F] # keep only high number of call SNPs [affy snps]
fff <- fff[,CS1$Calls>5000,drop=F] # keep only affy SNPs
hhh2 <- hhh2[,CS12$Calls>5000,drop=F] # keep only high number of call SNPs [illu snps]
fff2 <- fff2[,CS12$Calls>5000,drop=F] # keep only illu SNPs
#affy performs better than illumina at imputation
tt <- diag(cor(SnpMatrix.to.data.frame(fff[,]),SnpMatrix.to.data.frame(hhh[,]),use="pairwise.complete"))
summary> summary(tt)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7412  0.9620  0.9887  0.9674  0.9932  0.9997 
tt2 <- diag(cor(SnpMatrix.to.data.frame(fff2[,]),SnpMatrix.to.data.frame(hhh2[,]),use="pairwise.complete"))
summary(tt2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2106  0.5827  0.7294  0.6932  0.8437  0.9940 

csa <- col.summary(fff)
csb <- col.summary(hhh)
cor(csa$RAF,csb$RAF,use="pairwise.complete") # = .999 so no alignment issue.

csa2 <- col.summary(fff2)
csb2 <- col.summary(hhh2)
cor(csa2$RAF,csb2$RAF,use="pairwise.complete") # = .988 so no alignment issue.




### run separately for affy and illumina ids
if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) }
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))

all.affy.ids <- unique(c(aff.int.ids, aff.only.ids))
all.ill.ids <- unique(c(ill.int.ids, ill.only.ids))

setwd("~/barrett")
lambda.1000.affy <- lambda.1000.ill <- lambda.1000 <- numeric(22)
res.list.affy <- res.list.ill <- res.list <- vector("list",22)
for (j in 22:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.affy.ids,mc.cores=20)
  res.list.affy[[j]] <- do.call("rbind",args=res.list.nxt)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.ill.ids,mc.cores=20)
  res.list.ill[[j]] <- do.call("rbind",args=res.list.nxt)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,mc.cores=20)
  res.list[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.affy[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=4830,mr=1930)
  lambda.1000.affy[j] <- l1000
  cat("L1000 for Affy chr",j,"=",l1000,"\n")
  lambda.nxt <- round(median(p.to.Z(res.list.ill[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999,mr=3983)
  lambda.1000.ill[j] <- l1000
  cat("L1000 for Illumina chr",j,"=",l1000,"\n")
  lambda.nxt <- round(median(p.to.Z(res.list[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999+4830,mr=3983+1930)
  lambda.1000[j] <- l1000
  cat("L1000 for BOTH chr",j,"=",l1000,"\n")
  save(res.list.affy,res.list.ill,res.list,lambda.1000.affy,lambda.1000.ill,lambda.1000,file="/chiswick/data/ncooper/barrett/L1000Calcs_illVSaffy.RData")
}



#### analyse the affy ppl with some PCAs
library(reader)
library(humarray)
source("~/github/plumbCNV/SnpMatrixList.R")
source("~/github/imputer/barretFunctions.R")

if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) }
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))

all.affy.ids <- unique(c(aff.int.ids, aff.only.ids))
all.ill.ids <- unique(c(ill.int.ids, ill.only.ids))

setwd("~/barrett")
lambda.1000.affy.pc4 <- lambda.1000.affy.pc19  <- numeric(22)
res.list.affy.pc4 <- res.list.affy.pc19 <- vector("list",22)
for (j in 12:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.affy.ids,mc.cores=10,nPCs=4,
    PC.source.fn="~/barrett/AFFY/PCA/result.quick.affy.RData")
  res.list.affy.pc4[[j]] <- do.call("rbind",args=res.list.nxt)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.affy.ids,mc.cores=10,nPCs=19,
    PC.source.fn="~/barrett/AFFY/PCA/result.quick.affy.RData")
  res.list.affy.pc19[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.affy.pc4[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=4826,mr=1930)
  lambda.1000.affy.pc4[j] <- l1000
  cat("L1000 for Affy PC4 chr",j,"=",l1000,"\n")
  lambda.nxt <- round(median(p.to.Z(res.list.affy.pc19[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=4826,mr=1930)
  lambda.1000.affy.pc19[j] <- l1000
  cat("L1000 for Affy PC19 chr",j,"=",l1000,"\n") 
  save(res.list.affy.pc4,res.list.affy.pc19,lambda.1000.affy.pc4,lambda.1000.affy.pc19,
      file="/chiswick/data/ncooper/barrett/L1000Calcs_affyPC4vs19.RData")
}



#### analyse the ill ppl with some PCAs
if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) }
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))

all.affy.ids <- unique(c(aff.int.ids, aff.only.ids))
all.ill.ids <- unique(c(ill.int.ids, ill.only.ids))

setwd("~/barrett")
lambda.1000.ill.pc4 <- lambda.1000.ill.pc18 <- lambda.1000.ill.pc0 <- numeric(22)
res.list.ill.pc4 <- res.list.ill.pc18 <- res.list.ill.pc0 <- vector("list",22)
for (j in 11:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.ill.ids,mc.cores=20,nPCs=4,
    PC.source.fn="~/barrett/ILL/PCA/result.quick.ill.RData")
  res.list.ill.pc4[[j]] <- do.call("rbind",args=res.list.nxt)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.ill.ids,mc.cores=20,nPCs=18,
    PC.source.fn="~/barrett/ILL/PCA/result.quick.ill.RData")
  res.list.ill.pc18[[j]] <- do.call("rbind",args=res.list.nxt)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=sample.info,add.interp=TRUE,samp.subset=all.ill.ids,mc.cores=20,nPCs=0,
    PC.source.fn="~/barrett/ILL/PCA/result.quick.ill.RData")
  res.list.ill.pc0[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.ill.pc4[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999, mr=3983)
  lambda.1000.ill.pc4[j] <- l1000
  cat("L1000 for Illu PC4 chr",j,"=",l1000,"\n")
  lambda.nxt <- round(median(p.to.Z(res.list.ill.pc18[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999, mr=3983)
  lambda.1000.ill.pc18[j] <- l1000
  cat("L1000 for Illu PC18 chr",j,"=",l1000,"\n") 
  lambda.nxt <- round(median(p.to.Z(res.list.ill.pc0[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999, mr=3983)
  lambda.1000.ill.pc18[j] <- l1000
  cat("L1000 for Illu PC0 chr",j,"=",l1000,"\n") 
  save(res.list.ill.pc4,res.list.ill.pc18,res.list.ill.pc0,
    lambda.1000.ill.pc4,lambda.1000.ill.pc18,lambda.1000.ill.pc0,
      file="/chiswick/data/ncooper/barrett/L1000Calcs_illPC4vs18b.RData")
}



#### Run control versus Control ######

all.cont <- rownames(sample.info[sample.info$phenotype==1,])
# make new sample info where 'phenotype' contains affy vs ill codes instead
si2 <- sample.info[all.cont,]
si2$phenotype <- NA
si2$phenotype[rownames(si2) %in% all.affy.ids] <- 1
si2$phenotype[rownames(si2) %in% all.ill.ids] <- 2
si2[["phenocode"]] <- c("affy","illu")[si2$phenotype]
save(si2,file="/chiswick/data/ncooper/barrett/sampleInfoCtrlVsCtrl.RData")

# run a new loop

lambda.1000.ctrl <- numeric(22)
res.list.ctrl <- vector("list",22)
for (j in 22:1) {
  fnz <- get.impute.chr.fn("OUTPUT",j)
  res.list.nxt <- mclapply(as.list(fnz),analyse.impute.file,sample.info=si2,add.interp=TRUE,samp.subset=all.cont,mc.cores=20)
  res.list.ctrl[[j]] <- do.call("rbind",args=res.list.nxt)
  lambda.nxt <- round(median(p.to.Z(res.list.ctrl[[j]]$p)^2,na.rm=T)/.454,3)
  l1000 <- lambda_nm(lambda.nxt,nr=3999, mr=4830)
  lambda.1000.ctrl[j] <- l1000
  cat("L1000 for Affy vs Illu CTRLS only, chr",j,"=",l1000,"\n")
  save(res.list.ctrl,lambda.1000.ctrl,file="/chiswick/data/ncooper/barrett/L1000Calcs_CTRLAffVsIll.RData")
}



