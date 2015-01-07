ichip.dir <- "/chiswick/data/ncooper/iChipData/"
impute.dir <- "/chiswick/data/ncooper/imputation/common/"
data.dir <- paste0(impute.dir,"ALL.integrated_phase1_SHAPEIT_16-06-14.nosing/")
temp.txt <- "nextCHR.txt"

source("~/github/imputer/imputeFunctions.R")


options(ucsc="hg19")

#update allele info in ichip object
if(F) {
  
  rm(all.support)
  ichip37 <- chip.support(build=37)
  ichip36 <- chip.support(build=36)
  
  
print(load(cat.path(impute.dir,"all.rr.RData")))

#######

print(load("/chiswick/data/ncooper/iChipData/ImmunoChip_ChipInfo.RData"))
table(paste0(A1(ichip37[[19]]),A2(ichip37[[19]])),exclude=NULL)

bad.na1 <- which(is.na(A1(sup.for.impute)))
bad.na2 <- which(is.na(A2(sup.for.impute)))
any.bad <- unique(c(bad.na1,bad.na2))
if(length(any.bad)>0) {
  cat("missing alleles\n")
  print(sup.for.impute[any.bad,])
}
equal.ones <- which(A1(sup.for.impute)==A2(sup.for.impute))
if(length(equal.ones)>0) {
  cat("equal alleles\n")
  print(sup.for.impute[equal.ones,])
}

## REPORT ON ALL DUPLICATED POSITIONS ##
for(chrnum in c(1,2,5:13,16:19)) { show((sup[chr(sup)==chrnum,])[dup.pairs(start(sup[chr(sup)==chrnum])),]) }
## fix allele codes 
}

#mkdir imputation
#cd imputation
#wget https://mathgen.stats.ox.ac.uk/impute/ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz
#tar -zxvf ALL.integrated_phase1_SHAPEIT_16-06-14.nosing.tgz 

take.cls <- c("id","position","a0","a1","type","eur.aaf","eur.maf") # only take these cols from legend files

## EXCLUSION FILES ##
anc.excl <- reader("/chiswick/data/ncooper/iChipData/ancestryExclSamps.txt")
snp.excl <- reader("/chiswick/data/ncooper/iChipData/snpsExcluded.txt")
samp.excl <- unique(c(anc.excl,reader("/chiswick/data/ncooper/iChipData/sampsExcluded.txt")))

sup <- chip.support()

for (chrnum in 22:1) {

  cat("processing chromosome",chrnum,"...\n")

  ichip.fn <- cat.path(ichip.dir,fn="temp.ichip-data",suf=chrnum,ext="RData")
  next.rn <- cat.path(impute.dir,fn="CHR",suf=chrnum,ext="RData")
  next.fn <- cat.path(data.dir,pref="ALL.chr",fn=chrnum,suf=".integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend",ext="gz")

  posz <- start(sup)[chr(sup)==chrnum]
  rsz <- rs.id(sup)[chr(sup)==chrnum]

  if(!file.exists(next.rn) | T) {
    gt.for.impute <- get.SnpMatrix.in.file(ichip.fn,warn=FALSE)
    colnames(gt.for.impute) <- clean.snp.ids(colnames(gt.for.impute))
    if(chrnum==11){
      rsmtch <- which( colnames(gt.for.impute) %in% "imm_11_2138800")
      if(length(rsmtch)>0) {
        gt.for.impute[,rsmtch] <- rep(as.raw("00"),times=nrow(gt.for.impute)) # "imm_11_2138800" = "rs689"
        rs689.snpmatrix <- get(paste(load(cat.path(ichip.dir,"taqman_rs689.RData")))) #get.taqman.snp.from.text.file()
        indzz <- match(rownames(rs689.snpmatrix),rownames(gt.for.impute))
        gt.for.impute[narm(indzz),rsmtch] <- rs689.snpmatrix[which(!is.na(indzz)),1]
      } else {
        warning("did not find SNP rs689 in chromosome 11")
      }
    }
    system(paste("zcat",next.fn," > ",temp.txt))
    jj <- proc.time() ; rr <- reader(temp.txt); kk <- proc.time() ; print(kk[3]-jj[3])
    rr <- rr[,take.cls]
#    rids <- rr$id[(which((rr$position %in% posz | rr$id %in% rsz) & rr$type!="INDEL"))]
#    rpos <- rr$position[(which((rr$position %in% posz | rr$id %in% rsz) & rr$type!="INDEL"))]
    qids <- which((rr$position %in% posz) & rr$type!="INDEL")  # same rs-id and same position
    my.rr <- rr[qids,]
    print(length(qids))  #8746
#    indo <- match(rr$position[qids],posz)
    indo <- match(posz,my.rr$position)
    ichip.ids <- rsz[!is.na(indo)]
    dindo <- duplicated(indo)
    if(any(dindo)) {
      dups <- rsz[dup.pairs(indo)]
      warning(length(which(dindo))," duplicate positions were found in chromosome ",chrnum,
              ":\n  ",comma(head(dups,50)))
      #match(rr$position[qids],posz)
      #pdups <- duplicated(posz)
    }
    ichip.ids <- ichip.ids[!rs.to.id(ichip.ids) %in% rs.to.id(snp.excl)]
    #print(ichip.ids[grep("seq",ichip.ids)])
    #print(ichip.ids[grep("seq",rs.to.id(ichip.ids))])
    #print(rs.to.id(ichip.ids[grep("seq",rs.to.id(ichip.ids))]))
#   qids <- which(((rr$position %in% posz) & (rr$id %in% rsz)) & rr$type!="INDEL")  # same rs-id and same position
#   length(qids)  #8696
#   qids <- (which((((rr$position %in% posz) & (!rr$id %in% rsz))) & rr$type!="INDEL"))  # same position but different rs-id
#   length(qids) #1436
  
    indx <- narm(match(rs.to.id(ichip.ids),rs.to.id(colnames(gt.for.impute)))) 
    dd <- duplicated(indx); if(any(dd)) { warning(paste(indx[which(dd)],collapse=",")," were duplicated indices") }
    valid.snps <- colnames(gt.for.impute)[indx]
    dd <- duplicated(valid.snps); if(any(dd)) { stop(paste(valid.snps[which(dd)],collapse=",")," were duplicated SNPs") }
    samps <- (1:nrow(gt.for.impute))
    sei <- narm(match(samp.excl,rownames(gt.for.impute)))
    if(length(sei)>0) { samps <- samps[-sei] }
    gt.for.impute <- gt.for.impute[samps,indx]
    sup.for.impute <- sup[rs.to.id(valid.snps),]
    take.rws <- my.rr$type=="SNP"
    my.rr <- my.rr[take.rws,]
    rr <- my.rr[my.rr$position %in% start(sup.for.impute),]
    save(rr, gt.for.impute,sup.for.impute, file=next.rn)
    prv(rr, gt.for.impute,sup.for.impute)
  } else {
    print(load(next.rn))
    #str <- get.strands(rr, gt.for.impute, sup.for.impute)
  }
  fn <- paste0("CHR",chrnum)
  write.impute(gt.for.impute,pedfile=cat.path(data.dir,fn,ext=".gens"),bp=start(sup.for.impute), 
               a1=A1(sup.for.impute), a2=A2(sup.for.impute), snp.id=rownames(sup.for.impute))
}

