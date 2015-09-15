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



### PCA the 1000genomes ppl for all the GWAS snps ###



### MISSION102 ###

### Calc PCs for the 1000genomes ppl for the common GWAS snps ###
# check alignment of 1000g against our T1D set 
SML <- as.list(mixedsort(list.files("~/PLAY/aSnpMat",pattern="RData",full.names=T)))
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))
aff.ill.lookup <- cbind(ill.int.ids,aff.int.ids)

ready.to.impute.sml.dir <- "~/barrett/IMPUTE_INPUT"
set.dir <- "~/PLAY/"

new.dir <- "/chiswick/data/ncooper/barrett/AllPlusReplicates/"
SML2 <- as.list(cat.path(new.dir,"AllPlusReplicatesChr",suf=1:22,ext="RData"))

sample.info.1000 <- reader("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData")

for (cc in 1:length(SML)) {
  #sm <- reader(cat.path(ready.to.impute.sml.dir,"snpmat",suf=cc,ext="RData")) # barret newest data
  sm <- reader(SML2[[cc]])
  SML.nxt <- reader(SML[[cc]]) # 1000 genomes
  
  # get union list
  union.snps <- colnames(sm)[colnames(sm) %in% colnames(SML.nxt)]
  # take only the dual set ids from barrett, just for simplicity
  thousm <- SML.nxt[,union.snps] # only europeans used for align
  sm <- sm[,union.snps]
  # select same set in each
  # check alignment/align, should be quite similar
    # create alignment
  euros <- rownames(sample.info.1000)[sample.info.1000$ancestry %in% c("CEU","GBR","FIN","ESN")]
  ## align the europeans to our dataset
###pdf("/chiswick/data/ncooper/barrett/workingAlign.pdf")
  barr.align <- align.alleles( sm, thousm[sampsIn(thousm,euros),],mafdiff=0.05,do.plot=T)
###dev.off()
  #align.1000.all <- align.alleles(SML.nxt, sm, mafdiff=0.05,do.plot=F)
  #dd <- dups(align.1000.all,align.1000.eur)
  #align.1000.all.b <- align.alleles(align.1000.all, align.1000.eur, mafdiff=0.05,do.plot=F, known.dups=dd)
  save(barr.align,file=cat.path(set.dir,"barrettOnlyAlignedtoTG_Chr",suf=cc,ext="RData"))
  check.align.plot(barr.align,thousm,fn=cat.path("/chiswick/data/ncooper/","checkMe",suf=cc,ext="pdf"))
  loop.tracker(cc,22)
}
# plot 1000genomes ppl against the 1324 to see whether they match europeans

### MISSION101 ###

## retrieve full chr-wise datasets that have been updated with replicate affy samples ##
new.dir <- "/chiswick/data/ncooper/barrett/AllPlusReplicates/"
SML2 <- as.list(cat.path(new.dir,"AllPlusReplicatesChr",suf=1:22,ext="RData"))

new.dir <- "~/barrett"
set.dir <- "~/PLAY/"

all.dirs <- cat.path(set.dir,"barrettAlignedtoTG_Chr",suf=1:22,ext="RData")
for (cc in 1:length(all.dirs)) {
  (load(all.dirs[cc]))
  save(barr.align,file=cat.path(set.dir,"barrettOnlyAlignedtoTG_Chr",suf=cc,ext="RData"))
  loop.tracker(cc,length(all.dirs))
}
SML2 <- as.list(cat.path(set.dir,"barrettOnlyAlignedtoTG_Chr",suf=1:22,ext="RData"))

rs <- row.summary(barr.align)
# get samples typed on both chips
full.samp <- rownames(rs)[(which(rs$Call.rate>.99))]
# combine into large 1324 x 980711 SnpMatrix
big.snpMatrix <- sampSel(SML2,samples=full.samp)
save(big.snpMatrix,file="/chiswick/data/ncooper/barrett/allCompleteSamples1324.RData")
#big.snpMatrix <- reader("/chiswick/data/ncooper/barrett/allCompleteSamples1324.RData")

# LD-prune SNPs one chr at a time
chrz <- seqnames(snp.info[colnames(big.snpMatrix),])
chr.list <- unique(chrz)
init <- T
if(init) {
  res <- vector("list",length(chr.list))
  for (dd in 1:length(chr.list)) {
    Header(paste(dd))
    ii <- big.snpMatrix[,chrz==chr.list[dd]]
    res[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res,file="/chiswick/data/ncooper/barrett/prunedLists.RData")
  }
} else {
  res <- reader("/chiswick/data/ncooper/barrett/prunedLists.RData")
}

prune.snps <- unlist(lapply(res,names))
prune.snps <- unique(prune.snps)

intersect.snpMatrix <- big.snpMatrix
TGpruned <- intersect.snpMatrix[,prune.snps]

save(TGpruned,file="/chiswick/data/ncooper/barrett/TGpruned1324.RData")
#TGpruned <- reader("/chiswick/data/ncooper/barrett/TGpruned.RData")

num.to.use <- length(prune.snps)
TG.MR <- randomize.missing2(TGpruned,verbose=TRUE)
save(TG.MR,file="/chiswick/data/ncooper/barrett/TG.MR1324.RData")
#TG.MR <- reader("/chiswick/data/ncooper/barrett/TG.MR1324.RData")

big1000.GWA <- bigSnpMatrix(TG.MR,cat.path(new.dir,"big1000gA"),n.cores=8)
mySM <- describe(big1000.GWA)
save(mySM,file="/chiswick/data/ncooper/barrett/big1000gA1324.RData")

#big1000.GWA <- get.big.matrix(reader("/chiswick/data/ncooper/barrett/big1000gA1324.RData"))
rm.a <- colmean(big1000.GWA)

result.quick.A <- big.PCA(big.t(big1000.GWA),return.loadings=TRUE,center=rm.a)
save(result.quick.A,file="/chiswick/data/ncooper/barrett/result.quick.A1324.RData")
#result.quick.A <- reader("/chiswick/data/ncooper/barrett/result.quick.A1324.RData")

# now we have loadings, generated on the 1324 samples with all the SNPs.
# now impute the missing SNPS for everyone else and
# then apply the loadings to generate PCs for everyone
# IMPUTE SNPs

# APPLY PCs to this (barrett) dataset #
top5pc.A <- rev(order((1*abs(result.quick.A$loadings[,1]))+(1.00*abs(result.quick.A$loadings[,2]))))[1:num.to.use]
top5pc.An <- colnames(big1000.GWA)[top5pc.A]

mean.crct.t1d2 <- (big1000.GWA[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.GWA)))

n.pc <- 10
pca.pred.a.t1d2 <- as.data.frame(matrix(nrow=nrow(mean.crct.t1d2),ncol=n.pc))
colnames(pca.pred.a.t1d2) <- paste0("PC",1:n.pc)
for (cc in 1:n.pc) {
  pca.pred.a.t1d2[,cc] <-  mean.crct.t1d2 %*% result.quick.A$loadings[top5pc.An,cc]
}
rownames(pca.pred.a.t1d2) <- rownames(pca.pred.a.t1d2[,cc])

# APPLY to 1000 genomes (alternate) dataset #

TG.pruned <- snpSel(SML,snps=prune.snps)
TG.pruned2 <- snpSel(SML,snps=intersec.snps)

save(TG.pruned,file="/chiswick/data/ncooper/barrett/TGpruned1000g.RData")
save(TG.pruned2,file="/chiswick/data/ncooper/barrett/TGpruned1000gIntersect.RData")

#TG.pruned <- reader("/chiswick/data/ncooper/barrett/TGpruned1000g.RData")

chrz <- seqnames(snp.info[colnames(TG.pruned2),])
chr.list <- unique(chrz)
init <- T
if(init) {
  res <- vector("list",length(chr.list))
  for (dd in 1:length(chr.list)) {
    Header(paste(dd))
    ii <- TG.pruned2[,chrz==chr.list[dd]]
    res[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res,file="/chiswick/data/ncooper/barrett/prunedLists1000gI.RData")
  }
} else {
  res <- reader("/chiswick/data/ncooper/barrett/prunedLists1000gI.RData")
}
prune.snps <- unlist(lapply(res,names))
prune.snps <- unique(prune.snps)

num.to.use <- length(prune.snps)

TG.MR1000 <- randomize.missing2(TG.pruned,verbose=TRUE)
save(TG.MR1000,file="/chiswick/data/ncooper/barrett/TG.MR1000g.RData")
#TG.MR1000 <- reader("/chiswick/data/ncooper/barrett/TG.MR1000g.RData")
TG.MR1000i <- randomize.missing2(TG.pruned2[,prune.snps],verbose=TRUE)
save(TG.MR1000i,file="/chiswick/data/ncooper/barrett/TG.MR1000i.RData")

big1000.tg <- bigSnpMatrix(TG.MR1000,cat.path(new.dir,"big1000gTG"),n.cores=6)
big1000.tgi <- bigSnpMatrix(TG.MR1000i,cat.path(new.dir,"big1000gTGi"),n.cores=6)

mySM1 <- describe(big1000.tg)
save(mySM1,file="/chiswick/data/ncooper/barrett/big1000gTG.RData")
mySM2 <- describe(big1000.tgi)
save(mySM2,file="/chiswick/data/ncooper/barrett/big1000gTGi.RData")

rm.a <- colmean(big1000.tgi)

result.quick.A <- big.PCA(big.t(big1000.tgi),return.loadings=TRUE,center=rm.a)
save(result.quick.A,file="/chiswick/data/ncooper/barrett/result.quick.tgi.RData")

all(colnames(big1000.tg)==colnames(big1000.GWA))
mean.crct.tg <- (big1000.tg[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.tg)))
mean.crct.tgi <- (big1000.tgi[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.tgi)))

n.pc <- 10

pca.pred.a.tgi <- as.data.frame(matrix(nrow=nrow(mean.crct.tgi),ncol=n.pc))
colnames(pca.pred.a.tgi) <- paste0("PC",1:n.pc)
for (cc in 1:n.pc) {
  pca.pred.a.tgi[,cc] <-  mean.crct.tgi %*% result.quick.A$loadings[top5pc.An,cc]
}
rownames(pca.pred.a.tgi) <- rownames(pca.pred.a.tgi[,cc])
# DONE
t1d.snpmat <- snpSel(SML2,snps=top5pc.An)  #"~/PLAY/barrettOnlyAlignedtoTG_Chr1.RData"
T1d.MR <- randomize.missing2(t1d.snpmat,verbose=TRUE)
save(T1d.MR,file="/chiswick/data/ncooper/barrett/T1d.MR.RData")
T1d.MR.i <- bigSnpMatrix(T1d.MR,cat.path(new.dir,"T1d.MRi"),n.cores=4)
mySM3 <- describe(T1d.MR.i)
save(mySM3,file="/chiswick/data/ncooper/barrett/bigT1d.MRi.RData")

mean.crct.t1di <- (T1d.MR.i[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(T1d.MR.i)))
pca.pred.a.t1di <- as.data.frame(matrix(nrow=nrow(mean.crct.t1di),ncol=n.pc))
colnames(pca.pred.a.t1di) <- paste0("PC",1:n.pc)
for (cc in 1:n.pc) {
  pca.pred.a.t1di[,cc] <-  mean.crct.t1di %*% result.quick.A$loadings[top5pc.An,cc]
}
rownames(pca.pred.a.t1di) <- rownames(T1d.MR.i)

# Contains PCs created using the common affy/illumina SNPs, these were derived using the 1000
# genomes samples, and are for each T1D sample (which were aligned to 1000g), and for the 1000 genomes 
# samples in a second object with an extra ancestry column for plotting, etc.
PCs <- pca.pred.a.t1di
save(PCs,pca.pred.a.tgi,DESCR,file="/chiswick/data/ncooper/barrett/PCsFor1000gPCAIntersectSnpsForT1Dand1KG.RData")

pdf("/chiswick/data/ncooper/pca.plot.3.pdf"); 
with(pca.pred.a.t1di,plot(PC3,PC2)); 
dev.off()

## 1000 genomes plot
pca.pred.a.tgi[["anc"]] <- sample.info.1000$anc[match(rownames(pca.pred.a.tgi),rownames(sample.info.1000))]
pdf("/chiswick/data/ncooper/pca.plot.4.pdf"); 
with(pca.pred.a.tgi,plot(PC1,PC2,col=as.numeric(factor(anc)))); 
dev.off()


### NEW MISSION - AFFY PCA ###
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))
all.affy.ids <- unique(c(aff.int.ids, aff.only.ids))
dir.out <- "~/barrett/AFFY"
fnZ <- cat.path(paste0("~/barrett/","IMPUTE_INPUT"),fn="snpmat",suf=1:22,ext="RData") # SnpMatrix of data just before 'write.impute'
fnZ.out <- cat.path(dir.out,rmv.ext(basename(fnZ)),pref="AFFY",ext=".RData")
for (cc in 1:length(fnZ)) {
  X <- reader(fnZ[cc])
  snpmat <- X[sampsIn(X,all.affy.ids),snpsIn(X,all.affy)]
  save(snpmat,file=fnZ.out[cc])
  loop.tracker(cc,length(fnZ))
}

#TG.pruned2 <- snpSel(SML,snps=intersec.snps)
#save(TG.pruned2,file="/chiswick/data/ncooper/barrett/TGpruned1000gIntersect.RData")

#chrz <- seqnames(snp.info[colnames(TG.pruned2),])
init <- T
if(init) {
  res <- vector("list",22)
  for (dd in rev(1:22)) {
    Header(paste(dd))
    ii <- reader(fnZ.out[dd])
    res[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res,file="/chiswick/data/ncooper/barrett/prunedListsAffy.RData")
  }
} else {
  res <- reader("/chiswick/data/ncooper/barrett/prunedListsAffy.RData")
}
prune.snps <- unlist(lapply(res,names))
prune.snps <- unique(prune.snps)

num.to.use <- length(prune.snps)

new.dir <- "~/barrett/AFFY/PCA"

affy.all <- snpSel(as.list(fnZ.out),snps=prune.snps)

affy.MR <- randomize.missing2(affy.all,verbose=TRUE)
save(affy.MR,file="~/barrett/AFFY/PCA/affy.MR.RData")
#affy.MR <- reader("/chiswick/data/ncooper/barrett/AFFY/PCA/affy.MR.RData")

big.affy.MR <- bigSnpMatrix(affy.MR,cat.path(new.dir,"bigaffy.MR"),n.cores=6)

myAffy <- describe(big.affy.MR)
save(myAffy,file="~/barrett/AFFY/PCA/myAffy.RData")

rm.a <- colmean(big.affy.MR)

result.quick.A.affy <- big.PCA(big.t(big.affy.MR),return.loadings=TRUE,center=rm.a)
save(result.quick.A.affy,file="~/barrett/AFFY/PCA/result.quick.affy.RData")

###################
### NEW MISSION - ILLU PCA ###
(load("/chiswick/data/ncooper/barrett/IDLISTS.RData"))
all.ill.ids <- unique(c(ill.int.ids, ill.only.ids))
dir.out <- "~/barrett/ILL"
fnZ <- cat.path(paste0("~/barrett/","IMPUTE_INPUT"),fn="snpmat",suf=1:22,ext="RData") # SnpMatrix of data just before 'write.impute'
fnZ.out <- cat.path(dir.out,rmv.ext(basename(fnZ)),pref="ILL",ext=".RData")
for (cc in 1:length(fnZ)) {
  X <- reader(fnZ[cc])
  snpmat <- X[sampsIn(X,all.ill.ids),snpsIn(X,all.ill)]
  save(snpmat,file=fnZ.out[cc])
  loop.tracker(cc,length(fnZ))
}

#TG.pruned2 <- snpSel(SML,snps=intersec.snps)
#save(TG.pruned2,file="/chiswick/data/ncooper/barrett/TGpruned1000gIntersect.RData")

#chrz <- seqnames(snp.info[colnames(TG.pruned2),])
init <- T
if(init) {
  res <- vector("list",22)
  for (dd in 1:22) {
    Header(paste(dd))
    ii <- reader(fnZ.out[dd])
    res[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res,file="/chiswick/data/ncooper/barrett/prunedListsIll.RData")
  }
} else {
  res <- reader("/chiswick/data/ncooper/barrett/prunedListsIll.RData")
}
prune.snps <- unlist(lapply(res,names))
prune.snps <- unique(prune.snps)

num.to.use <- length(prune.snps)

new.dir <- "~/barrett/ILL/PCA"

ill.all <- snpSel(as.list(fnZ.out),snps=prune.snps)

ill.MR <- randomize.missing2(ill.all,verbose=TRUE)
save(ill.MR,file="~/barrett/ILL/PCA/ill.MR.RData")
#ill.MR <- reader("/chiswick/data/ncooper/barrett/ILL/PCA/ill.MR.RData")

big.ill.MR <- bigSnpMatrix(ill.MR,cat.path(new.dir,"bigill.MR"),n.cores=6)

myIll <- describe(big.ill.MR)
save(myIll,file="~/barrett/ILL/PCA/myIll.RData")

rm.a <- colmean(big.ill.MR)

result.quick.A.ill <- big.PCA(big.t(big.ill.MR),return.loadings=TRUE,center=rm.a)
save(result.quick.A.ill,file="~/barrett/ILL/PCA/result.quick.ill.RData")
#############
### MISSION103 ###

#ALTERNATIVE USING COMMON SNPs between PLATFORMS
intersect.snpMatrix <- snpSel(SML2,snps=intersec.snps)
save(intersect.snpMatrix,file="/chiswick/data/ncooper/barrett/allSamplesIntersectingSNPs.RData")
#intersect.snpMatrix <- reader("/chiswick/data/ncooper/barrett/allSamplesIntersectingSNPs.RData")

# LD-prune SNPs one chr at a time
chrz <- seqnames(snp.info[colnames(intersect.snpMatrix),])
chr.list <- unique(chrz)
init <- T
if(init) {
  res2 <- vector("list",length(chr.list))
  for (dd in 1:length(chr.list)) {
    Header(paste(dd))
    ii <- intersect.snpMatrix[,chrz==chr.list[dd]]
    res2[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res2,file="/chiswick/data/ncooper/barrett/prunedListsCommonSnps.RData")
  }
} else {
  res2 <- reader("/chiswick/data/ncooper/barrett/prunedListsCommonSnps.RData")
}


prune.snps2 <- unlist(lapply(res2,names))
prune.snps2 <- unique(prune.snps2)

#mini.prune <- sample(prune.snps,round(38000*.75))

TGpruned <- intersect.snpMatrix[,prune.snps2]

save(TGpruned,file="/chiswick/data/ncooper/barrett/TGpruned.RData")
#TGpruned <- reader("/chiswick/data/ncooper/barrett/TGpruned.RData")

num.to.use <- length(prune.snps2)

TG.MR <- randomize.missing2(TGpruned,verbose=TRUE)
save(TG.MR,file="/chiswick/data/ncooper/barrett/TG.MR.RData")
#TG.MR <- reader("/chiswick/data/ncooper/barrett/TG.MR.RData")

big1000.GWA <- bigSnpMatrix(TG.MR,cat.path(new.dir,"big1000gA"))
inter.big <- describe(big1000.GWA)
save(inter.big,file="/chiswick/data/ncooper/barrett/big1000gA.RData")

big1000.GWA <- get.big.matrix(reader("/chiswick/data/ncooper/barrett/big1000gA.RData"))
rm.a <- colmean(big1000.GWA)

result.quick.A <- big.PCA(big.t(big1000.GWA),return.loadings=TRUE,center=rm.a)
save(result.quick.A,file="/chiswick/data/ncooper/barrett/result.quick.A.RData")
result.quick.A <- reader("/chiswick/data/ncooper/barrett/result.quick.A.RData")


### Get the top 'n' loading SNP list ###
# divide 1 by 2 of vec: estimate.eig.vpcs(result.quick.A$Evalues,M=big.t(big1000.A))$variance.pcs[1:2]
top5pc.A <- rev(order((1*abs(result.quick.A$loadings[,1]))+(1.00*abs(result.quick.A$loadings[,2]))))[1:num.to.use]
top5pc.An <- colnames(big1000.GWA)[top5pc.A]
#top5pc.An <- top5pc.An[!top5pc.An %in% snp.excl]
#save(top5pc.A,top5pc.An,result.quick.A, sample.info,pruned.1000g.MR,SM, file="pruned.list.PCA.1000g.RData")
#load("pruned.list.PCA.1000g.RData")
#rm.a <- colMeans(big1000.A[,top5pc.An])
############################################################# 

pca1.pred.a <- scale(big1000.GWA[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.An,1]
pca2.pred.a <- scale(big1000.GWA[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.An,2]
pca.bsm.a.pred.1 <-   (BSM.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(BSM.A))) %*% result.quick.A$loadings[top5pc.An,1]
pca.bsm.a.pred.2 <-   (BSM.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(BSM.A))) %*% result.quick.A$loadings[top5pc.An,2]

######## MAKE DATA FRAME OF T1D 10 PCs #########

mean.crct.t1d <- (big1000.GWA[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.GWA)))

n.pc <- 10
pca.pred.a.t1d <- as.data.frame(matrix(nrow=nrow(mean.crct.t1d),ncol=n.pc))
colnames(pca.pred.a.t1d) <- paste0("PC",1:n.pc)
for (cc in 1:n.pc) {
  pca.pred.a.t1d[,cc] <-  mean.crct.t1d %*% result.quick.A$loadings[top5pc.An,cc]
}
rownames(pca.pred.a.t1d) <- rownames(pca.pred.a.t1d[,cc])

# correlations are all 1 - same dataset, ths is good!
for (jj in 1:10) { print(cor(result.quick.A$PCs[,jj],pca.pred.a.t1d[,jj])) }
############################################################# 

########

# set to bigMatrix type, with missing replaced, transformations
new.dir <- "~/barrett/";# setwd(new.dir) # 5 cores max, else memory gets a bit overloaded
bigSnpMat <- bigSnpMatrix(big.snpMatrix,filename=cat.path(new.dir,"tempMatrix"), n.cores=5, limit.ram=F, mn.zero=T, sd.hwe=T, replace.missing=T) 
  

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


