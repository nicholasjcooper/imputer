
do.plot <- TRUE # whether to plot lambda1000 for each disease
num.pcs <- 10 # number of PCs to test
## note that several TRUE/FALSE options can be set in the 'single-snp-tests.R' file, namely:
#    debug.in.r.session, fail.if.any.fail, do.ancestry.exclude, exclude.pc34, verbose 

#### SCRIPT TO GENERATE DATA FOR 1 TO 10 PCs ########
library(gtools)
library(reader)
source("/home/ncooper/github/plumbCNV/FunctionsCNVAnalysis.R") # source to get Nick's extra functions loaded

# location of the input annotSnpStats datafiles for all diseases combined, by CHR/pq-arm
dir.raw <- "/chiswick/data/ncooper/imputation/COMBINED/PQDATA/"
# location that single-snp-tests.R will write *.RData files to:
dir.postpc <- "/chiswick/data/ncooper/imputation/COMBINED/postpca-single-snp-tests/"
# location to write temporary qsub/grid log files:
dir.grid <- "/chiswick/data/ncooper/imputation/COMBINED/GRID" 
# location for output of plots/data files
dir.out <- "/chiswick/data/ncooper/imputation/PCcovTests/"

### set up named list to store the results for lambda for each DISEASE x PC-n ###
dis <- paste0("p.",c("celiac","t1d","jia","ms","atd","ra")) # colnames for p.values in 'results' objects
num.pcs <- 10;  pcs <- numeric(num.pcs)
lambda.list <- list(pcs) 
lambda.list <- lambda.list[rep(1,times=length(dis))]
names(lambda.list) <- dis

### Run single-snp-tests.R for a model with each of 1 thru 10 PC covariates ###
for (n.pc in c(1:num.pcs)) {
  Header(paste("PC =",n.pc),h = "@")

  ## delete all the files from the previous run ##
  to.delete <- list.files(dir.postpc,pattern="RData")
  to.delete <- cat.path(dir.postpc,to.delete)
  unlink(to.delete) # remove last lot of RData files before writing new ones
  
  ## generate the bash command to run single-snp-tests.R for each CHR/PQ-arm ##
  fns <- list.files(dir.raw,pattern="RData")
  fns <- mixedsort(fns) # sort names to genome order
  cmds <- paste0("/home/ncooper/github/imputer/single-snp-tests.R disease=COMBINED ",
                 "ifile=",cat.path(dir.raw,fns)," pc=",n.pc)
  
  #### submit run of single-snp-tests.R to the queue in parallel for each chromosome/arm ####
  bash.qsub(cmds, dir=dir.grid, interval=15)  # NB: this function has its code in FunctionsCNVAnalysis.R
  
  # script will save the results summary for each chromosome to these locations:
  datl <- list.files(dir.postpc,pattern="RData")
  datl <- mixedsort(datl)
  
  ## read the result for the p arm of Chromosome 1 ##
  results <- reader(cat.path(dir.postpc,datl[1]))
  
  ## bind the results together for all chromosome/pq-arm's
  RES <- results
  for (cc in 2:length(datl)) { 
    results <- reader(cat.path(dir.postpc,datl[cc]))
    RES <- rbind(RES,results) 
  }
  
  # get James' null snp set
  null_p <- reader("/chiswick/data/ncooper/imputation/COMMON/null_p.RData")
  
  ### calculate inflation and lambda 1000 for each disease ###
  for (dd in 1:length(dis)) { 
    X <- (RES[null_p,dis[dd]])  # use only NULL snps for lambda calc
    #X <- (RES[,dis[dd]])  # uncomment to use ALL snps for lambda calc
    chisq <- qchisq(1-X,1)
    lambda.list[[dd]][n.pc] <- lambda <- median(chisq,na.rm=T)/0.4549364
    cat(dis[dd],", lambda =",lambda,"\n")
    # do qq-plot for next DISEASE x PC-n combination #
    pdf(cat.path(dir.out,paste0("QQtest_PC",n.pc,"_",dis[dd]),ext="pdf"))
      qqplot(sort(X),runif(length(X))); 
    dev.off() 
  }
  ## Save the combined results for 'n' PCs to a binary file ##
  save(RES,file=cat.path(dir.out,"resultFor",suf=paste0(n.pc,"PCs"),ext="RData"))
}

########################################

###### plot of lamdas ######
if(do.plot) {
  lam.range <- c(0.75,1.5) # observed - manually set limits for the plot
  pdf(cat.path(dir.out,"inflationVsLambda1000K.pdf"))
    ## N's for case/controls for lambda 1000 calculations ##
    dis.ncases <- c(6317,6551,1214,4379,2694,1171) # for celiac,t1d,jia,ms,atd,ra respectively
    dis.nctrls <- c(12443,12443,12443,12443,12443,7234)
    plot(1:num.pcs,runif(num.pcs,lam.range[1],lam.range[2]),col="white",ylim=lam.range,
       xlab="number of PCs covaried for",ylab="lambda1000",
       main="inflation for varying levels of PC-correction across 6 AI diseases")
    for (dd in 1:length(dis)) {
      l1000 <- lambda_nm(lambda.list[[dd]],nr=dis.ncases[dd],mr=dis.nctrls[dd])
      lines(1:num.pcs,l1000,col=dd,lwd=2)  
    }
    legend("topright",legend=dis,col=1:length(dis),lwd=2,bty="n")
  dev.off()
}  
############################





#l1000 <- lambda_nm(lambda.list,nr=6317,mr=12443) #coeliac
#l1000 <- lambda_nm(Lnm,nr=6551,mr=12443) #t1d
#l1000 <- lambda_nm(Lnm,nr=1214,mr=12443) #jia
#l1000 <- lambda_nm(Lnm,nr=4379,mr=12443) #ms
##l1000 <- lambda_nm(Lnm,nr=2694,mr=12443) #atd
#l1000 <- lambda_nm(Lnm,nr=11171,mr=7234) #ra

#affected
#phenotype   1     2
#COELIAC     0  6317
#RACTRL   7234     0
#UKCTRL  12443     0
#GRAVES      0  2694
#JIA         0   307
#MS      ?9358  4379
#RA          0 11171
#T1D         0  6551

