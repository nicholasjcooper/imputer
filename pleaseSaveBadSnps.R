#!/usr/bin/Rscript
files <- list.files("/stats/chrisw/FM-ss/impdata",pattern=".RData",full=TRUE)

pvars <- c("p.celiac", "p.t1d", "p.jia", "p.ms", "p.atd", "p.ra")

p0 <- p2 <- matrix(NA,length(files),6,dimnames=list(basename(files),pvars))
outlist <- list()
for(i in seq_along(files)) {
  (load(files[[i]]))
  p0[i,] <- apply(subset(results,type==0,select=pvars),2,min,na.rm=TRUE)
  p2[i,] <- apply(subset(results,type==2,select=pvars),2,min,na.rm=TRUE)

}

# type2  = genotyped, type0 = imputed
p0 <- pmax(p0,p0,1e-50)
p2 <- pmax(p2,p2,1e-50)

library(reshape)
m0 <- melt(p0)
m2 <- melt(p2)
p02 <- merge(m2,m0,by=c("X1","X2"))

library(ggplot2)
ggplot(p02,aes(x=-log10(value.x),y=-log10(value.y))) + geom_abline() + geom_point() + labs(x="genotyped (type 2)",y="imputed (type 0)") + facet_wrap(~X2)

bad.dis <- vector("list",length(pvars))
names(bad.dis) <- pvars

for (dd in seq_along(pvars)) {
  l0 <- (-log10(p0[,pvars[dd]]))
  l2 <- (-log10(p2[,pvars[dd]]))
  difz <- (l0 - l2)
  ratz <- difz/l2
  #fz <- which.outlier(difz,thr=1.5,low=F)
  fz <- which(difz>10 | (ratz>2 & difz>5))

  minz <- (-log10(p0[,pvars[dd]]))[fz]

  thr <- .2  # percentage within the outlying value to tag

  bad.rows <- vector("list",length(fz))
  cat(length(fz),pvars[dd],"\n")
  for(cc in seq_along(fz)) {
    i <- fz[cc]
    (load(files[[i]]))
    P2 <- subset(results,type==0,select=pvars[dd])
    P3 <- -log10(P2)
    bad.rows[[cc]] <- P3[which(P3>.8*minz[1]),,drop=F]
    #badz <- which(P3>(1-thr)*minz[1])
    #bad.snps <- rownames(P3)[badz]
  }
  all.rows <- do.call("rbind",args=bad.rows)
  bad.dis[[dd]] <- if(is.null(all.rows)) { NA } else { all.rows }
}
#bad.dis
sapply(bad.dis,nrow)


