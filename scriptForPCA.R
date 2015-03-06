#pca 1000 genomes for gwas

init <- TRUE # FALSE 
setwd("/chiswick/data/ncooper/imputation/THOUSAND/")

library(reader); #source("~/github/iChip/iFunctions.R"); source("~/github/plumbCNV/FunctionsCNVAnalysis.R")
library(humarray); library(bigpca)
source("~/github/plumbCNV/SnpMatrixList.R")
source("~/github/imputer/imputeFunctions.R"); library(gtools)

sml <- mixedsort(list.files("/chiswick/data/ncooper/imputation/THOUSAND/TG12M.files/"))
sml <- cat.path("/chiswick/data/ncooper/imputation/THOUSAND/TG12M.files/",sml)
sml <- as.list(sml)

if(init) {
  res <- reader("/chiswick/data/ncooper/imputation/THOUSAND/prunedLists.RData")
  #res <- vector("list",39)
  for (dd in 1:1) { #length(res)) {
    Header(paste(dd))
    ii <- get.SnpMatrix.in.file(sml[[dd]])
    res[[dd]] <- ld.prune.big(ii,thresh=.1,n.cores=16)
    save(res,file="/chiswick/data/ncooper/imputation/THOUSAND/prunedLists.RData")
  }
} else {
  res <- reader("/chiswick/data/ncooper/imputation/THOUSAND/prunedLists.RData")
}

prune.snps <- unlist(lapply(res,names))
prune.snps <- unique(prune.snps)
if(init) { writeLines(prune.snps,con="/chiswick/data/ncooper/imputation/THOUSAND/noLDsnpList1000g.txt") }

## get cases which only have 660K (~half) SNPs, and keep only the SNPs that are in both arrays for PCA ##
ms2.dir <- "/chiswick/data/ncooper/imputation/MS/WTCCC2/Cases_UKN/"
casel <- mixedsort(list.files(ms2.dir,pattern="RData"))
casel <- cat.path(ms2.dir,casel)
case.msml <- as.list(casel) # only 379 samples
cn <- colnamesL(case.msml,list=F)
prune.snps <- prune.snps[prune.snps %in% cn]
if(init) { writeLines(prune.snps,con="/chiswick/data/ncooper/imputation/THOUSAND/noLDsnpList1000g660.txt") }
# CN <- colnames(TG.aligned.pca)
# TG.aligned.pca <- TG.aligned.pca[,which(CN %in% cn)] # now have manageable 1092 x 76962 matrix for PCA

## get the MS 1958BC controls to align to ##
msdir <- "/chiswick/data/ncooper/imputation/MS/WTCCC2/Controls_Illu58C/"
msml <- mixedsort(list.files(msdir,pattern="RData"))
msml <- cat.path(msdir,msml)
msml <- as.list(msml)
sub.msml <- sampSel(msml,samples=sample(2930,1000))

TGpruned <- snpSel(sml,prune.snps)
which(!colnames(TGpruned) %in% colnames(sub.msml))
# [1] 61901 # if this is here, need to remove using line below
# TGpruned <- TGpruned[,-61901]

# remove TG uncertain genotypes
sn <- (snps(TGpruned))
ja <- paste(sn$allele.1,sn$allele.2,sep="")
# remove uncertain genotypes #
uncertain <- ja %in% c("AT","TA","CG","GC")
TGpruned <- TGpruned[,-which(uncertain)]

msmlDat <- sub.msml[,colnames(TGpruned)]

# remove MS uncertain genotypes
sn <- (snps(msmlDat))
ja <- paste(sn$allele.1,sn$allele.2,sep="")
# remove uncertain genotypes #
uncertain <- ja %in% c("AT","TA","CG","GC")
msmlDat <- msmlDat[,-which(uncertain)]

TGpruned <- TGpruned[,colnames(msmlDat)]




sample.info <- reader("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData")
anc <- factor(sample.info$ancestry[match(rownames(TGpruned),rownames(sample.info))])
TG.gbr <- TGpruned[anc=="GBR",]
pdf("chrisplotMS.pdf") ; TG.gbr.aligned <- align.alleles(TG.gbr,msmlDat,mafdiff=.01); dev.off()

tg.raf <- col.summary(TG.gbr.aligned)$RAF
ms.raf <- col.summary(msmlDat)$RAF
ms.raf[is.na(ms.raf)] <- 0
tg.raf[is.na(tg.raf)] <- 0
bad.match <- which(abs(tg.raf-ms.raf)>0.1)
TGpruned2 <- TGpruned[,-bad.match]
msmlDat2 <- msmlDat[,-bad.match]
TG.gbr2 <- TGpruned2[anc=="GBR",]
pdf("chrisplotMS3.pdf") ; TG.gbr.aligned2 <- align.alleles(TG.gbr2,msmlDat2,mafdiff=.01); dev.off()


pdf("chrisplotMS4.pdf") ; TG.aligned3 <- align.alleles(TGpruned2,TG.gbr.aligned2,mafdiff=.01); dev.off()



### RUN PCA on 1000 genomes GWAS and check it looks right ###
# replace missing and create bigSnpMatrix #
pruned.TG.MR <- randomize.missing2(TG.aligned3,verbose=TRUE)
big1000.GWA <- bigSnpMatrix(pruned.TG.MR,"big1000gA")
rm.a <- colmean(big1000.GWA)
result.quick.A <- big.PCA(big.t(big1000.GWA),return.loadings=TRUE,center=rm.a)
#anc <- factor(sample.info$ancestry[match(rownames(big1000.GWA),rownames(sample.info))])
pdf("PCA1000.Gwas.pdf")
 plot(result.quick.A$PCs[,"PC1"],result.quick.A$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
 legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4)
dev.off()
############################################################# 

snp.excl <- readLines("/chiswick/data/ncooper/imputation/MS/WTCCC2/snpsPreExcludedAll.txt")

### Get the top 'n' loading SNP list ###
# divide 1 by 2 of vec: estimate.eig.vpcs(result.quick.A$Evalues,M=big.t(big1000.GWA))$variance.pcs[1:2]
top5pc.A <- rev(order((1*abs(result.quick.A$loadings[,1]))+(1.00*abs(result.quick.A$loadings[,2]))))[1:31578]
top5pc.An <- colnames(big1000.GWA)[top5pc.A]
top5pc.An <- top5pc.An[!top5pc.An %in% snp.excl]
save(top5pc.A, top5pc.An, result.quick.A, anc,sample.info,pruned.TG.MR, file="pruned.list.PCA.1000g.gwas.RData")
load("pruned.list.PCA.1000g.gwas.RData")
rm.a <- colMeans(big1000.GWA[,top5pc.An])
############################################################# 






############################################################# 


#### GET GWAS DATA FOR ALL THE MS GROUPS! ######

gwas.sets <- c("Controls_Illu58C","Controls_IlluNBS","Cases_UKC","Cases_UKN","Cases_UKP","Cases_UKW") #[-1]
n.pc <- 10
rm.a <- colmean(big1000.GWA) # get 1000 genomes means for each SNP
top5pc.A <- match(top5pc.An,colnames(big1000.GWA)) # index instead of name
dir.th <- "/chiswick/data/ncooper/imputation/THOUSAND"
align <- TRUE


for (dd in 1:length(gwas.sets)) {
  cat("Loading SnpMatrixList for ",gwas.sets[dd],"...")
  sml <- get.sml(gwas.sets[dd])
  SM <- snpSel(sml,snps=top5pc.An)
  cat("done\n")
  bad.ids.a <- rmv.common.missing(X=SM,Y=pruned.TG.MR)
  top5pc.Bn <- top5pc.An[!clean.snp.ids(top5pc.An) %in% clean.snp.ids(bad.ids.a)]
  top5pc.B <- match(top5pc.Bn,colnames(big1000.GWA)) # index instead of name
  SM <- SM[,top5pc.Bn]
  ## align the alleles to the original T1D dataset ###
  if(align) {
    cat("Aligning aSnpMatrix to the 1000 genomes GBR dataset...")
    baseTG <- TG.gbr2[,top5pc.Bn]
    pdf(cat.path(dir.th,"AlignPlot",suf=gwas.sets[dd],ext="pdf")) 
    SM2 <- align.alleles(SM,baseTG,mafdiff=.05); dev.off()
    cat("done\n")
  } else { SM2 <- SM }
  cat("creating bigSnpMatrix...")
  BSM.A <- bigSnpMatrix(SM2,paste0("bigPCA",gwas.sets[dd]),ref.data=pruned.TG.MR,n.cores=10)
  cat("done\n")
  mean.crct.dis <- (BSM.A[,top5pc.Bn]-rep(rm.a[top5pc.B],each=nrow(BSM.A))) # mean corrected SNP matrix per disease
  cat("calculating ",n.pc," principal components...")
  pca.pred.a.dis <- as.data.frame(matrix(nrow=nrow(mean.crct.dis),ncol=n.pc))
  colnames(pca.pred.a.dis) <- paste0("PC",1:n.pc)
  for (cc in 1:n.pc) {
    pca.pred.a.dis[,cc] <-  mean.crct.dis %*% result.quick.A$loadings[top5pc.Bn,cc]
  }
  rownames(pca.pred.a.dis) <- rownames(pca.pred.a.dis[,cc])
  PCs <- pca.pred.a.dis
  cat("done\n")
  ofn <- cat.path(dir.th,"PCs",suf=gwas.sets[dd],ext="RData")
  save(PCs,file=ofn)
  cat("wrote file: ",ofn,"\n")
  align <- TRUE

if(dd==1) {
  ######## MAKE DATA FRAME OF 1000 GENOMES top GWAS 10 PCs #########
  mean.crct.tg <- (big1000.GWA[,top5pc.Bn]-rep(rm.a[top5pc.Bn],each=nrow(big1000.GWA)))

  n.pc <- 10
  pca.pred.a.tg <- as.data.frame(matrix(nrow=nrow(mean.crct.tg),ncol=n.pc+1))
  colnames(pca.pred.a.tg) <- c(paste0("PC",1:n.pc),"ancestry")
  for (cc in 1:n.pc) {
    pca.pred.a.tg[,cc] <-  mean.crct.tg %*% result.quick.A$loadings[top5pc.Bn,cc]
    corz <- cor(pca.pred.a.tg[,cc],result.quick.A$PCs[,cc],use="pairwise.complete")
    cat("correlation of PC",cc,"with projection was: ",corz,"\n")
  }
  rownames(pca.pred.a.tg) <- rownames(pca.pred.a.tg[,cc])
  pca.pred.a.tg[["ancestry"]] <- sample.info$ancestry[match(rownames(pca.pred.a.tg),rownames(sample.info))]
}

  pdf(cat.path(dir.th,"PCA",suf=gwas.sets[dd],ext="pdf")); 
   xl <- range(c(PCs[,1],pca.pred.a.tg[,1])[-1]); yl <- range(c(PCs[,2],pca.pred.a.tg[,2])[-1])
   plot(pca.pred.a.tg[,1],pca.pred.a.tg[,2],col=get.distinct.cols(14)[as.numeric(anc)] ,ylim=yl, xlim=xl)
   points(PCs[,1],PCs[,2],col="black",pch="+")
   legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4)
  dev.off()
}


## AT a loss as to why this is not quite working. thinking about reimporting the 1000 genomes data, as perhaps the colnames got mixed up at some point??? ###














###########
#  ICHIP  #
###########


# to PCA 1000 genomes for available iChip Snps
# assume already extracted 1000 genomes for ichip and aligned using annotSnpStats to our dataset (T1D-ichip)
(load("/chiswick/data/ncooper/imputation/THOUSAND/aligned1000g.RData"))
newsml <- SnpMatrix.to.sml.by.chr(asm.1000g.aligned)
ld.pruned <- smlapply(newsml,ld.prune.chr,thresh=0.1,n.cores=12) # this step is slow (several hours?)
pruned.1000g <- NULL
pruned.1000g <- asm.1000g.aligned[,names(unlist(ld.pruned))]
save(ld.pruned,pruned.1000g,file="pruned.list1000g.RData")
#save(ld.pruned,file="pruned.list_gwas_1000g.RData")

lsml <- sampSel(sml,samples=sample(15000,1000),dir=dd)
sample.info <- reader("/chiswick/data/ncooper/iChipData/sample.info.RData")
asmT1D <- annot.sep.support(snpMat=lsml,snp.info=ichip.local, sample.info=sample.info)
(load("/chiswick/data/ncooper/imputation/THOUSAND/thousandGenomesAndIChip.RData"))
colnames(asm.1000g) <- clean.snp.ids(colnames(asm.1000g))
colnames(asmT1D) <- clean.snp.ids(colnames(asmT1D))
comn <- colnames(asm.1000g)[(colnames(asm.1000g) %in% colnames(asmT1D))]
ii <- match(comn,colnames(asmT1D))
asmT1D <- asmT1D[,ii]
asm.1000g <- asm.1000g[,comn]
asm.1000gbr <- asm.1000g[which(asm.1000g@samples$ancestry %in% c("GBR","CEU","FIN")),]

pdf("chrisplot.pdf") ; asm.1000gbr.aligned <- align.alleles(asm.1000gbr,asmT1D,mafdiff=.05); dev.off()
pdf("chrisplot.pdf") ; asm.1000.aligned <- align.alleles(asm.1000g,asm.1000gbr.aligned,mafdiff=.05); dev.off()
save(asm.1000.aligned,"1000gAlignedToT1D.RData")
pruned.1000g <- asm.1000.aligned[,names(unlist(ld.pruned))]

### Align 1000 genomes to T1D (as T1D is biggest) ###
# start here not to be very slow #
(load("pruned.list1000g.RData"))
sn <- (snps(pruned.1000g))
ja <- paste(sn$allele.1,sn$allele.2,sep="")
# remove uncertain genotypes #
uncertain <- ja %in% c("AT","TA","CG","GC")
pruned.1000g <- pruned.1000g[,-which(uncertain)]
(load("/chiswick/data/ncooper/imputation/THOUSAND/sample.info.1000g.RData"))
# replace missing and create bigSnpMatrix #
pruned.1000g.MR <- randomize.missing2(pruned.1000g,verbose=TRUE)
big1000.A <- bigSnpMatrix(pruned.1000g.MR,"big1000gA")
############################################################# 

### RUN PCA on 1000 genomes and check it looks right ###
rm.a <- colmean(big1000.A)
result.quick.A <- big.PCA(big.t(big1000.A),return.loadings=TRUE,center=rm.a)
anc <- factor(sample.info$ancestry[match(rownames(result.quick.A$PCs),rownames(sample.info))])
pdf("PCA1000.A.pdf"); plot(result.quick.A$PCs[,"PC1"],result.quick.A$PCs[,"PC2"],col=get.distinct.cols(14)[as.numeric(anc)],ylim=c(-.05,.1))
legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()
############################################################# 


#### find the SNPs excluded in any of the datasets ####
snp.excl.files <- c("/chiswick/data/ncooper/imputation/RA/snpsExcludedAll.txt",
"/chiswick/data/ncooper/imputation/JIA/snpsExcluded2.txt",
"/chiswick/data/ncooper/imputation/GRAVES/snpsExcluded.txt",
"/chiswick/data/ncooper/imputation/COELIAC/snpsExcluded.txt",
"/chiswick/data/ncooper/imputation/MS/iCHIP/snpsExcludedComb.txt",
"/chiswick/data/ncooper/imputation/T1D/QC/snpsExcludedPREV.txt")

snp.excl <- unique(unlist(sapply(snp.excl.files,readLines)))
############################


### Get the top 'n' loading SNP list ###
# divide 1 by 2 of vec: estimate.eig.vpcs(result.quick.A$Evalues,M=big.t(big1000.A))$variance.pcs[1:2]
top5pc.A <- rev(order((1*abs(result.quick.A$loadings[,1]))+(1.00*abs(result.quick.A$loadings[,2]))))[1:11578]
top5pc.An <- colnames(big1000.A)[top5pc.A]
top5pc.An <- top5pc.An[!top5pc.An %in% snp.excl]
save(top5pc.A,top5pc.An,result.quick.A, anc,sample.info,pruned.1000g.MR,SM, file="pruned.list.PCA.1000g.RData")
load("pruned.list.PCA.1000g.RData")
rm.a <- colMeans(big1000.A[,top5pc.An])
############################################################# 


### GET all T1D data ###
dd <- ("/chiswick/data/ncooper/imputation/T1D/PQDATA/")
sml <- list.files(dd)
sml <- as.list(sml)
sml <- sml[c(18:19,26:40,1:17,20:24)]
SM <- snpSel(sml,snps=top5pc.An,dir=dd)
bad.ids.a <- rmv.common.missing(X=SM,Y=pruned.1000g)
top5pc.An <- top5pc.An[!clean.snp.ids(top5pc.An) %in% clean.snp.ids(bad.ids.a)]
SM <- SM[,top5pc.An]
SM.big <- SM; SM <- SM[sample(nrow(SM),200),]
BSM.A <- bigSnpMatrix(SM,"bigPCAa",ref.data=pruned.1000g.MR)

######## MAKE DATA FRAME OF T1D 10 PCs #########

mean.crct.t1d <- (BSM.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(BSM.A)))

n.pc <- 10
pca.pred.a.t1d <- as.data.frame(matrix(nrow=nrow(mean.crct.t1d),ncol=n.pc))
colnames(pca.pred.a.t1d) <- paste0("PC",1:n.pc0)
for (cc in 1:n.pc) {
  pca.pred.a.t1d[,cc] <-  mean.crct.t1d %*% result.quick.A$loadings[top5pc.An,cc]
}
rownames(pca.pred.a.t1d) <- rownames(pca.pred.a.t1d[,cc])

############################################################# 


######## MAKE DATA FRAME OF 1000 GENOMES top 10 PCs #########
mean.crct.tg <- (big1000.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.A)))

n.pc <- 10
pca.pred.a.tg <- as.data.frame(matrix(nrow=nrow(mean.crct.tg),ncol=n.pc+1))
colnames(pca.pred.a.tg) <- c(paste0("PC",1:n.pc),"ancestry")
for (cc in 1:n.pc) {
  pca.pred.a.tg[,cc] <-  mean.crct.tg %*% result.quick.A$loadings[top5pc.An,cc]
  corz <- cor(pca.pred.a.tg[,cc],result.quick.A$PCs[,cc],use="pairwise.complete")
  cat("correlation of PC",cc,"with projection was: ",corz,"\n")
}
rownames(pca.pred.a.tg) <- rownames(pca.pred.a.tg[,cc])
pca.pred.a.tg[["ancestry"]] <- sample.info$ancestry[match(rownames(pca.pred.a.tg),rownames(sample.info))]

#pca2.pred.a <-   (big1000.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(big1000.A))) %*% result.quick.A$loadings[top5pc.An,2]

save(pca.pred.a.tg,pca.pred.a.t1d,file="/chiswick/data/ncooper/imputation/THOUSAND/AncestryPCsT1DandTG.RData")
############################################################# 


######## MAKE PLOT OF 1000 GENOMES top 10 PCs vs T1D #########
pca1.pred.a <- scale(big1000.A[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.An,1]
pca2.pred.a <- scale(big1000.A[,top5pc.An],center=T,scale=F) %*% result.quick.A$loadings[top5pc.An,2]
pca.bsm.a.pred.1 <-   (BSM.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(BSM.A))) %*% result.quick.A$loadings[top5pc.An,1]
pca.bsm.a.pred.2 <-   (BSM.A[,top5pc.An]-rep(rm.a[top5pc.An],each=nrow(BSM.A))) %*% result.quick.A$loadings[top5pc.An,2]


pdf("PCAT1D.A5.pdf"); 
xl <- range(c(pca.pred.a.t1d[,1],pca.pred.a.tg[,1])); yl <- range(c(pca.pred.a.t1d[,2],pca.pred.a.tg[,2]))
plot(pca.pred.a.tg[,1],pca.pred.a.tg[,2],col=get.distinct.cols(14)[as.numeric(anc)] ,ylim=yl, xlim=xl)
points(pca.pred.a.t1d[,1],pca.pred.a.t1d[,2],col="black",pch="+")
legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()


pdf("PCAT1D.A.pdf"); 
xl <- range(c(pca.bsm.a.pred.1,pca1.pred.a)); yl <- range(c(pca.bsm.a.pred.2,pca2.pred.a))
plot(pca1.pred.a,pca2.pred.a,col=get.distinct.cols(14)[as.numeric(anc)] ,ylim=yl, xlim=xl)
points(pca.bsm.a.pred.1,pca.bsm.a.pred.2,col="black",pch="+")
legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()



######## MAKE PLOT OF ALIGNMENT #########

pdf("chrisplot.pdf") ; asm.1000g.aligned2 <- align.alleles(asm.1000g,asmT1D,mafdiff=.05); dev.off()


#### GET DATA FOR ALL THE DISEASES! ######

ichip.sets <- c("T1D","COELIAC","GRAVES","JIA","RA1","RA2","RA3","RA4","RA5","RA6","IC.CASE","IC.CTRL") #[-1]
n.pc <- 10
rm.a <- colmean(big1000.A) # get 1000 genomes means for each SNP
top5pc.A <- match(top5pc.An,colnames(big1000.A)) # index instead of name
dir.th <- "/chiswick/data/ncooper/imputation/THOUSAND"
align <- FALSE

for (dd in 2:length(ichip.sets)) {
  cat("Loading SnpMatrixList for ",ichip.sets[dd],"...")
  sml <- get.sml(ichip.sets[dd])
  SM <- snpSel(sml,snps=top5pc.An)
  cat("done\n")
  bad.ids.a <- rmv.common.missing(X=SM,Y=pruned.1000g.MR)
  top5pc.Bn <- top5pc.An[!clean.snp.ids(top5pc.An) %in% clean.snp.ids(bad.ids.a)]
  top5pc.B <- match(top5pc.Bn,colnames(big1000.A)) # index instead of name
  SM <- SM[,top5pc.Bn]
  ## align the alleles to the original T1D dataset ###
  if(align) {
    cat("Aligning aSnpMatrix to the T1D dataset...")
    baseT1D <- asmT1D[,top5pc.Bn]
    pdf(cat.path(dir.th,"AlignPlot",suf=ichip.sets[dd],ext="pdf")) 
    SM2 <- align.alleles(SM,baseT1D,mafdiff=.05); dev.off()
    cat("done\n")
  } else { SM2 <- SM }
  cat("creating bigSnpMatrix...")
  BSM.A <- bigSnpMatrix(SM2,paste0("bigPCA",ichip.sets[dd]),ref.data=pruned.1000g.MR,n.cores=10)
  cat("done\n")
  mean.crct.dis <- (BSM.A[,top5pc.Bn]-rep(rm.a[top5pc.B],each=nrow(BSM.A))) # mean corrected SNP matrix per disease
  cat("calculating ",n.pc," principal components...")
  pca.pred.a.dis <- as.data.frame(matrix(nrow=nrow(mean.crct.dis),ncol=n.pc))
  colnames(pca.pred.a.dis) <- paste0("PC",1:n.pc)
  for (cc in 1:n.pc) {
    pca.pred.a.dis[,cc] <-  mean.crct.dis %*% result.quick.A$loadings[top5pc.Bn,cc]
  }
  rownames(pca.pred.a.dis) <- rownames(pca.pred.a.dis[,cc])
  PCs <- pca.pred.a.dis
  cat("done\n")
  ofn <- cat.path(dir.th,"PCs",suf=ichip.sets[dd],ext="RData")
  save(PCs,file=ofn)
  cat("wrote file: ",ofn,"\n")
  align <- TRUE
  pdf(cat.path(dir.th,"PCA",suf=ichip.sets[dd],ext="pdf")); 
  xl <- range(c(PCs[,1],pca.pred.a.tg[,1])); yl <- range(c(PCs[,2],pca.pred.a.tg[,2]))
  plot(pca.pred.a.tg[,1],pca.pred.a.tg[,2],col=get.distinct.cols(14)[as.numeric(anc)] ,ylim=yl, xlim=xl)
  points(PCs[,1],PCs[,2],col="black",pch="+")
  legend("top",legend=paste(unique(anc)),col=get.distinct.cols(14)[as.numeric(unique(anc))],pch=19,ncol=4); dev.off()


}



