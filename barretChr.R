

Header(paste("Chr",cc))


## files we will make
#f.plot <- paste0(ndir,diseases,"/allelematch-",chrz[cc],args$arm,".png")

f.impute <- cat.path(ndir,paste0("impute-",chrz[cc]))
f.snps <- cat.path(ndir,paste0("snps-dropped-for-impute-",chrz[cc],".csv.gz"))

## load IMPUTE 1000 genomes Phase 3 legend file
head(legend <- read.table(f.legend[cc],header=TRUE,as.is=TRUE))
LEG <- legend
LEG[["chr"]] <- cc
rownames(LEG) <- LEG[,1]
LEG <- df.to.ranged(LEG,start="position",end="position")
#################


XX <- vector("list",length(cohorts))
for(i in seq_along(cohorts)) {
  message(cohorts[i])
  #snp.info2 is for AFFY, snp.info is for illumina
  XX[[i]] <- annot.sep.support(get.SnpMatrix.in.file(sml[[i]][[cc]]),if(i!=3) { snp.info2 } else { snp.info },sample.info)
  # must read snp.info, sample.info then annot.sep.support
  #See /dunwich/scratch/chrisw/T1DGC/R-objects/corrected-support/ for the SNP support objects
  #See /ipswich/data/T1DGC/R-objects/support/ for the sample support objects
}

cr.thr <- c(.95,.95,.98)

for(i in seq_along(cohorts)) {
  csumm <- col.summary(XX[[i]])
  XX[[i]] <- XX[[i]][, csumm[,"Call.rate"] > cr.thr[i] & abs(csumm[,"z.HWE"])<5 ]
}


if(FALSE) {
  i <- 1
  # fn <- cat.path(getwd(),"callRateDistributions",suf=i,ext="pdf")
  ii <- doSnpQC(sml[[i]],sample.info=sample.info,snp.info=snp.info)
  hwe.density.plots(Z.hwe=ii$SNP.INFO$Z.hwe,hwe.thr=5,dir=getwd(),fn=cat.path("","hwe",suf=i,ext="pdf"))
  hwe.vs.callrate.plots(call.rate=ii$SNP.INFO$call.rate,Z.hwe=ii$SNP.INFO$Z.hwe,
                      callrate.snp.thr=.95,hwe.thr=5,dir=getwd())
  jj <- doSampQC(sml[[i]],het.lo=.29,het.hi=.34)
  draw.density.plots(fn,jj$SAMPLE.INFO,ii$SNP.INFO)

}


## check SNP overlap
for(i in 1:(length(cohorts)-1)) {
  for(j in (i+1):length(cohorts)) {
    ci <- make.names(colnames(XX[[i]]))
    cj <- make.names(colnames(XX[[j]]))
    cat(cohorts[i],cohorts[j],length(ci),length(cj),length(intersect(ci,cj)),"\n",sep="\t")
  }
}

## align and merge
pdf(cat.path(home.dir,"alignplot",suf=cc,ext="pdf"))

par(mfrow=c(2,2))
#for(j in 1:(length(cohorts)-1)) {
j <- 1
  for(i in (j+1):length(cohorts)) {
    message()
    message(cohorts[[i]], " ", cohorts[[j]]) ## in each case 1 was 'j'
    m <- match(colnames(XX[[i]]),colnames(XX[[j]])) ## in each case 1 was 'j'
    d1 <- XX[[i]][,!is.na(m)]
    d2 <- XX[[j]][,m[!is.na(m)]]  ## in each case 1 was 'j'
    # print(col.summary(d1)["rs12251307",])
    # print(col.summary(d2)["rs12251307",])
    XX[[i]][,!is.na(m)] <- align.alleles(d1, d2, mafdiff=0.05)  # was 'sw'
    # print(col.summary(sw)["rs12251307",])
    title(sub=paste(cohorts[[i]], "-", cohorts[[j]]))  ## in each case 1 was 'j'
  }
#}
dev.off()

cat("combining 3 cohorts to one file\n")
YY <- do.call("sync.asnp.mats", args=XX)
X <- gt.for.impute <- do.call("rbind3",args=YY)
cat("new dimension for chr ",cc," is",dim(X)[1],",",dim(X)[2],"\n")

alleles <- snps(X)
alleles$RAF <- col.summary(X)[,"RAF"]
si <- snp.info.from.annot(X)
leg.uid <- paste(chrm(LEG),start(LEG),sep="_")
si.uid <- paste(chrm(si),start(si),sep="_")
## match against legend
summary(m <- match(si.uid,leg.uid))
message("Attempted to match:\t",length(m))
message("Found:\t",sum(!is.na(m)))
message("Missing:\t",sum(is.na(m)))
index.disease <- which(!is.na(m))
index.legend <- m[!is.na(m)]

alleles <- alleles[index.disease,]

Header("HERE")
print(length(index.disease))
print(Dim(alleles))



## compare alleles
x.alleles <- paste(alleles$allele.1,alleles$allele.2,sep="/") #[index.disease]
y.alleles <-  paste(legend$a0,legend$a1,sep="/")[index.legend]
message("allele codes before switching")

print(tt <- as.matrix(table(disease=x.alleles,legend=y.alleles)))
alleles$legend.alleles <- alleles$legend.id <- NA
alleles$legend.alleles <- y.alleles   #alleles$legend.alleles[index.disease] <- y.alleles
alleles$legend.id <- legend$id[index.legend]  #alleles$legend.id[index.disease] <- legend$id[index.legend]


## genotype classes
message("SNP classes")
table(sw.class <- g.class(x.alleles,y.alleles))

## unambiguous switches
sw <- sw.class %in% c("rev","revcomp")
any.comp <- any(sw.class %in% c("comp","revcomp"))
any.ambig <- any(sw.class=="ambig")

## any impossible?
any.impossible <- any(sw.class=="impossible")
n.impossible <- 0
if(any.impossible) {
  wh <- which(sw.class=="impossible")
  n.impossible <- length(wh)
  message(n.impossible, " impossible matches found, they will be removed:")
  #print(cbind(snp=alleles$rs.id[index.disease][wh],disease=x.alleles[wh],legend=y.alleles[wh]))
  prv(cbind(snp=alleles$rs.id[wh],disease=x.alleles[wh],legend=y.alleles[wh]))
  sw[wh] <- NA
}

if(any.comp & any.ambig) { 
  # there are reverse complements in the distinguishable cases
  message("We observed both strand switches and allele reversals. Resolving ambiguous SNPs using MAF")
  ind <- which(sw.class=="ambig")
  ind.ok <- which(!(sw.class %in% c("ambig","impossible")))
  message(length(ind)," SNPs have alleles not completely resolvable without strand information")
  #x.raf <- ifelse(sw,1-alleles[index.disease,"RAF"],alleles[index.disease,"RAF"])
  x.raf <- ifelse(sw,1-alleles[,"RAF"],alleles[,"RAF"])
  y.raf <- legend[index.legend,grep("EUR",colnames(legend),ignore.case=T)[1]] #"eur.aaf"]
  
  ## generate mahalanobis distances from "good" points
  S <- cov(x.raf[ind.ok],y.raf[ind.ok])
  d <- (x.raf[ind] - y.raf[ind])^2/S
  d.sw <- (1 - x.raf[ind] - y.raf[ind])^2/S
  d.limit <- min(0.1, 1.1 * max( (x.raf[ind.ok] - y.raf[ind.ok])^2/S )) # set a sensible upper limit of 0.1 in case of occasional oddities
  sw[ind][ d.sw < d ] <- TRUE
  sw[ind][ d.sw < d.limit & d < d.limit & x.raf[ind]>0.40 & x.raf[ind]<0.6 ] <- NA
}


## snp list
#alleles$sw[index.disease] <- sw
#alleles$sw.class[index.disease] <- sw.class
alleles$sw <- sw
alleles$sw.class <- sw.class

## do switch and write impute
print(Dim(alleles))
print(Dim(X))
X <- gt.for.impute[,rownames(alleles)]
print(Dim(alleles))
print(Dim(X))
X <- X[,!is.na(sw)]
X <- switch.alleles(X,sw[!is.na(sw)])
print(Dim(alleles))
print(Dim(X))

alleles <- alleles[!is.na(sw),]
table(alleles$sw.class,rownames(alleles) %in% colnames(X))

print(Dim(alleles))
print(Dim(X))

## duplicates
drop <- which(duplicated(alleles$position))
if(length(drop)) {
  cat(length(drop),"duplicated SNPs being dropped\n")
  alleles <- alleles[-drop,]
  X <- X[,-drop]
}

csumm <- col.summary(X)
drop <- which(csumm[,"Call.rate"]<0.05)
if(length(drop)) {
  cat(length(drop),"SNPs with >95% missingness being dropped\n")
  alleles <- alleles[-drop,]
  X <- X[,-drop]
}


drop <- which(!substr(alleles$legend.alleles,1,1) %in% c("A","C","G","T"))
if(length(drop)) {
  cat(length(drop),"SNPs with invalid NON-REF allele codes being dropped\n")
  alleles <- alleles[-drop,]
  X <- X[,-drop]
}

drop <- which(!substr(alleles$legend.alleles,3,3) %in% c("A","C","G","T"))
if(length(drop)) {
  cat(length(drop),"SNPs with invalid REF allele codes being dropped\n")
  alleles <- alleles[-drop,]
  X <- X[,-drop]
}

drop <- which(substr(alleles$legend.alleles,3,3)==substr(alleles$legend.alleles,1,1))
if(length(drop)) {
  cat(length(drop),"SNPs with impossible allele codes being dropped\n")
  alleles <- alleles[-drop,]
  X <- X[,-drop]
}

if(write.imputes) {
  #prv(f.snps,f.impute)
  write.table(alleles,file=gzfile(f.snps),quote=FALSE,sep="\t",row.names=FALSE)
  #write.impute(X, a1=substr(alleles$legend.alleles,1,1),a2=substr(alleles$legend.alleles,3,3),bp=alleles$position, fileroot=f.impute,
  #rs.id=alleles$legend.id, snp.id=alleles$Name)
  cat("\nwrote summary file for chromosome",cc,"to ",f.snps,"\n")
  write.impute(X, a1=substr(alleles$legend.alleles,1,1),a2=substr(alleles$legend.alleles,3,3),bp=alleles$position, pedfile=f.impute, snp.id=rownames(alleles))
  cat("wrote impute file for chromosome",cc,"to ",f.impute,"\n")
}

if(write.snpmatrix) {
  ofn <- cat.path(paste0(base.dir,"IMPUTE_INPUT"),fn="snpmat",suf=cc,ext="RData")
  save(X,file=ofn)
  cat("wrote aSnpMatrix object for chromosome",cc,"to",ofn,"\n")
 # save(X,file=cat.path(paste0("/chiswick/data/ncooper/barrett","IMPUTE_INPUT"),fn="snpmat",suf=cc,ext="RData"))
}
