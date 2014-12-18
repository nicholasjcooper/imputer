#!/usr/bin/Rscript


## Paste the following into terminal to run once for each CHR/pq-arm ##
#
# d=/chiswick/data/ncooper/imputation/COMBINED/PQDATA
# for f in $d/*.RData; do
#   /home/ncooper/github/imputer/single-snp-tests.R disease=COMBINED ifile=$f pc=5
# done
##

### Set run options ##
debug.in.r.session <- FALSE # set TRUE if doing test runs in an R-session rather than running from bash cmd-line #
fail.if.any.fail <- FALSE  # set TRUE to exclude any SNP that fails in any cohort, rather than just the current one
do.ancestry.exclude <- TRUE # set TRUE to exclude outliers on principle components for ancestry
exclude.pc34 <- FALSE # set TRUE to exclude outliers on PCs 3,4 in addition to the default PCs 1,2
verbose <- TRUE # set TRUE to print values for manual checking of successful running

library(reader)
library(snpStats)

## get command arguments ##
if(!debug.in.r.session) {
  defaults <- list(disease="COMBINED",ifile="/chiswick/data/ncooper/imputation/COMBINED/PQDATA/CHR1p.RData", pc="fail")
  args <- parse.args(coms=names(defaults),def=defaults,list.out=T,verbose=T)
} else {
  args <- list(disease ="COMBINED",ifile = "/chiswick/data/ncooper/imputation/COMBINED/PQDATA/CHR10p.RData",pc= "5")
  print(args)
}

setwd("/home/chrisw/Projects/FineMapping")

message("!!! single-snp-tests.R")


## files to read
d <- "/chiswick/data/ncooper/imputation"
f.in <- args$ifile
n.pc <- args$pc
# file with QC/ancestry-PC failure codes for each sample, with final index column mapping ids to annotSnpStats rownames #
sample.pc.fail.support <- "/chiswick/data/ncooper/imputation/COMMON/SAMPLES/iChipOverallSampleSupport.RData"
file.with.all.PCs <- "/chiswick/data/ncooper/imputation/THOUSAND/allGroupsPCs.RData"

## files to create
od <- file.path(d,args$disease,"postpca-single-snp-tests")
if(!file.exists(od))
  dir.create(od)
f.out <- file.path(od,sub(".RData",".csv.gz",basename(args$ifile)))

## read data
library(devtools)
load_all("/home/chrisw/RP/annotSnpStats")
(load(f.in))
if(!is.list(X))
  X <- list(X)

## tests
nulls <- sapply(X,is.null)
if(any(nulls)) 
  X <- X[which(!nulls)]
snps <- do.call("rbind",lapply(X,snps))
sX <- lapply(X,sm)


### define main function to run single SNP tests ###
do.test <- function(Y) {
  ss <- vector("list",length(sinfo))
  for(i in seq_along(X)) {
    message(i)
    # make selection vector of controls plus current phenotype that are in the list of sample IDs passing QC #
    use <- which(sinfo[[i]]$phenotype %in% c("CONTROL",Y) & (rownames(sinfo[[i]]) %in% pass.qc))
    if(verbose) { print(length(use)); print(max(use)) }
    
    if(fail.if.any.fail) {
      # exclude any SNP that fails in any cohort:
      filt <- all.filt
    } else {
      ## exclude SNPs that fail for the current phenotype:
      filt <- readLines(cat.path("/chiswick/data/ncooper/imputation/COMMON/SNPS/","snpsExcluded",ext="txt",suf=Y)) 
    }

    if(!all(rownames(sX[[i]])==rownames(sinfo[[i]]))) { stop("rownames in SnpMatrix did not match rownames in sample.info") }
    tmp.X <- sX[[i]][use,] # get aSnpMatrix for current CHR/pq-arm and sample-set
    tmp.X <- tmp.X[,!colnames(tmp.X) %in% filt] # filter any QC/PC-QC snps
    
    ## if only using failures for current cohort, some low call-rate SNPs are getting through
    # so here run a quick exclusion for low call rate SNPS just in case #
    pp <- col.summary(tmp.X)
    fail.cr.snps <- rownames(pp)[(which(pp$Call.rate<.94))]
    tmp.X <- tmp.X[,!colnames(tmp.X) %in% fail.cr.snps] # filter any extra call-rate failing snps
    ######
    
    # set up covariates for each sample #
    tmp.S <- sinfo[[i]][use,]
    if(verbose) {  print(dim(tmp.X));  print(dim(tmp.S)) }  #check that dims make sense
    ph <- tmp.S$affected # phenotype vector
    region <- factor(tmp.S$country) # region covariate (only used for RA)
    country <- tmp.S$country
    PC.DAT <- all.pc[match(rownames(tmp.X),rownames(all.pc)),1:n.pc,drop=F] # PC-components matrix for sample-set and n-PCs
    
    ## run model (with modified specs if current phenotype is RA) ##
    npcs <- paste(paste0("PC",1:n.pc),collapse=" + ") # PC(s) text for formula string
    ctrl <- glm.test.control(maxit = 50, epsilon = 1.e-5, R2Max = 0.999)
    if(Y=="RA") {
      form1 <- paste("ph ~ region +",npcs); prv(form1); form1 <- as.formula(form1)
      tmp.ss <- snp.rhs.tests(formula=form1,snp.data=tmp.X,data=PC.DAT,subset=(country!="UK"),control=ctrl)
    } else {
      form2 <- paste("ph ~",npcs); prv(form2) ;form2 <- as.formula(form2)
      tmp.ss <- snp.rhs.tests(formula=form2,snp.data=tmp.X,data=PC.DAT,subset=(country=="UK"),control=ctrl)
    }
    ss[[i]] <- tmp.ss # store the results for the next set
  }
  if(verbose) { 
    print(lapply(ss,head)) # little preview of p-value/chisq to ensure it worked #
  }
  p <- lapply(ss,p.value) #,1)
  p <- unlist(p)  
  if(length(p)!=ncol(tmp.X)) { stop("number of p-values should match number of columns in data") }
  names(p) <- colnames(tmp.X) # p-values are now connected to snp-names
  return(p)
}
###############


#### MAIN PROGRAM ####

sinfo <- lapply(X,samples) # generate sample.info

# generate country and phenotype columns in sample.info
if(!("country" %in% colnames(sinfo[[1]]))) {
  sinfo <- lapply(sinfo, function(s) {
    s$country <- ifelse(grepl("RA.",s$file,fixed=TRUE),sub("RA.","",s$file),"UK")
    s$phenotype <- ifelse(s$affected==1,"CONTROL",sub("\\..*","",s$file))
    return(s)
  })
}

### need to match support with PCA sample rownames to chris' aSnpMatrix IDs:
# this will allow extraction of a list of QC/PC-QC passing samples, plus PCs matrix for each cohort
si <- reader(sample.pc.fail.support) # load support file with PC-QC/QC flags
all.pc <- reader(file.with.all.PCs) # load PCs for all diseases

xids.of.si <- match(rownames(X[[1]]),si$X.id) # match chris' aSnpMatrix object IDs to custom id column in support file
SI <- si[xids.of.si,] # create new support object with same row order as 'samples(X[[1]])'
rownames(SI) <- rownames(X[[1]]) # replace rownames to be same as X[[1]]

# create vector of sample ids passing regular QC and PC-ancestry QC
# pca.fail <2 only excludes outliers on PC1,PC2; pca.fail==0 would also exclude outliers on PC3,PC4
if(do.ancestry.exclude) {
  if(exclude.pc34) {
    # exclude regular QC failers, and outliers on any of PCs 1-4
    pass.qc <- rownames(SI)[SI$pca.fail==0 & SI$QCfail==0] 
  } else {
    # exclude regular QC failers, and outliers on PCs 1,2
    pass.qc <- rownames(SI)[SI$pca.fail<2 & SI$QCfail==0] 
  }
} else {
  # just regular QC
  pass.qc <- rownames(SI)[SI$QCfail==0]
}

all.pc <- all.pc[rownames(si)[match(rownames(X[[1]]),si$X.id)],] # now we have the matrix of PCs for each sample
prv(all.pc) # uncomment to preview PC data.frame
################

if(fail.if.any.fail) {
  ## create a[n optional] filter for any SNP excluded in ANY disease ##
  all.filt <- NULL
  for (grp in c("COELIAC","T1D","JIA","MS","GRAVES","RA")) {
    all.filt <- unique(c(all.filt,readLines(cat.path("/chiswick/data/ncooper/imputation/COMMON/SNPS/",
                                                     "snpsExcluded",ext="txt",suf=grp))))
  }
}

### Run the main function for each disease sequentially ###
p.celiac <- do.test(Y="COELIAC")
p.t1d <- do.test("T1D")
p.jia <- do.test("JIA")
p.ms <- do.test("MS")
p.atd <- do.test("GRAVES")
p.ra <- do.test("RA")

# create vector of any SNP that appears in ANY disease results set #
all.names <- unique(c(names(p.celiac),names(p.t1d),names(p.jia),names(p.ms),names(p.atd),names(p.ra)))

# cbind the results together, indexing by SNPs named list to ensure rows are in sync #
results <- cbind(snps[match(all.names,rownames(snps)),],
                 p.celiac=p.celiac[all.names],
                 p.t1d=p.t1d[all.names],
                 p.jia=p.jia[all.names],
                 p.ms=p.ms[all.names],
                 p.atd=p.atd[all.names],
                 p.ra=p.ra[all.names])

results <- results[order(results$position),]
                 

## write results as compressed csv files ##
write.table(results,file=gzfile(f.out),sep="\t",quote=FALSE,row.names=FALSE)
## write results as RData binary files ##
ofn <- sub(".csv.gz",".RData",f.out)
save(results,file=ofn)
cat("wrote file: ",ofn,"\n")


## make GG-manhatten plots ##
library(reshape)
head(x <- melt(results[,grep("^p.",colnames(results))],"position"))
x$variable <- toupper(sub("p.","",x$variable))
x$flag <- 0
x$flag[x$value<5e-8] <- 1
x$flag[x$value<1e-100] <- 2
levels(x$flag <- factor(x$flag))
x$value[x$value<1e-100] <- 1e-100
library(ggplot2)
pdf(sub(".csv.gz",".pdf",f.out),height=8,width=8)
  ggplot(x,aes(x=position,y=-log10(value),colour=flag,pch=flag)) + geom_point() + facet_grid(variable ~ .,scales="free_y") + ylab("-log10(p)") + xlab(paste("Chr",args$chr)) + theme(legend.position="none")
dev.off()


##############################



## Chris' old code, not sure what difference:


## sinfo <- lapply(sinfo,function(s) {colnames(s) <- sub("phenotype","y",colnames(s)); return(s)})
## ss.celiac <- mapply(single.snp.tests,snp.data=sX,data=sinfo,
##                     MoreArgs=list(phenotype="affected",
##                       subset=which(y %in% c("CONTROL","COELIAC"))))

## results <- vector("list",length(sX))
## for(i in seq_along(sX))  {
## head(ss.celiac <- single.snp.tests(snp.data=sX[[i]],phenotype=sinfo[[i]]$affected,
##                               subset=which(sinfo[[i]]$y %in% c("CONTROL","COELIAC"))))
## ss.t1d <- single.snp.tests(snp.data=sX[[i]],phenotype="affected",data=sinfo[[i]],
##                               subset=phenotype %in% c("CONTROL","T1D"))
## ss.jia <- single.snp.tests(snp.data=sX[[i]],phenotype="affected",data=sinfo[[i]],
##                               subset=phenotype %in% c("CONTROL","JIA"))
## ss.ms <- single.snp.tests(snp.data=sX[[i]],phenotype="affected",data=sinfo[[i]],
##                               subset=phenotype %in% c("CONTROL","MS"))
## ss.atd <- single.snp.tests(snp.data=sX[[i]],phenotype="affected",data=sinfo[[i]],
##                               subset=phenotype %in% c("CONTROL","GRAVES"))
## ss.ra <- single.snp.tests(snp.data=sX[[i]],phenotype="affected",data=sinfo[[i]],
##                               subset=phenotype %in% c("CONTROL","RA"))
## results[[i]] <- cbind(chr=args$chr,
##                  snps[[i]],
##                  p.celiac=p.value(ss.celiac,1),
##                  p.t1d=p.value(ss.t1d,1),
##                  p.jia=p.value(ss.jia,1),
##                  p.ms=p.value(ss.ms,1),
##                  p.atd=p.value(ss.atd,1),
##                  p.ra=p.value(ss.ra,1))
## }
