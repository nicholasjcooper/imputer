
www@dan.de

## functions used, now in 'barretFunctions.R'
# is.valid.rn()
# get.vcf.header.from.files()
# convert.leg.hap.samp.to.vcf() # not used here but inspired by this exercise
# extract.leg.hap.samp.1000g.p1
# extract.leg.hap.samp.1000g.p2


if(!exists("sample.info")) { (load("/chiswick/data/ncooper/imputation/COMMON/allsupport.RData")) } 
# ie, get sample.info, snp.info, snp.info2

####### EXTRACT 1000 genomes for any chip, e.g GWAS (set of SNPs that appears in 'snp.info') #######
# FOR PHASE 3 in LEG, HAP format:

# initially interleaved the bash commands with the R commands but actually you can run as:
# 1) ALL R part 1 : extract.leg.hap.samp.1000g.p1()
# 2) ALL Bash : /usr/bin/awk ... etc
# 3) ALL R part 2 : extract.leg.hap.samp.1000g.p2()

# initialise all file name vectors
legz <- paste0("1000GP_Phase3_chr",1:22,".legend.gz")
hapz <- paste0("1000GP_Phase3_chr",1:22,".hap.gz")
out.fnz <- paste0("myChipChr",1:22)

# run part 1 in R to generate the bash code (must paste rather than use 'system' for strange technical reasons)
# also prepares chromosome position files
ll <- vector("list",22)
for (j in 1:22) {
  ll[[j]] <- extract.leg.hap.samp.1000g.p1(leg=legz[j],hap=hapz[j],out.fn=out.fnz[j],snp.info=snp.info,chr=j)
}

## NOW RUN BASH TO CREATE ALMOST READY FILES (no headers) ##

## Run R to add the header, and extract the file as a SnpMatrix and add annotation
# to make an aSnpMatrix for each chromosome
wait(3,"h")
for (j in 3) {
  cat("processing part 2 of extraction for chr: ",j,"..")
  asm <- extract.leg.hap.samp.1000g.p2(ll[[j]])
  save(asm,file=cat.path("~/PLAY","asmchr",suf=j,ext="RData"))
  cat("done\n")
}

## after some cleaning the resulting files are: ~/PLAY/aSnpMat
SML <- as.list(mixedsort(list.files("~/PLAY/aSnpMat",pattern="RData",full.names=T)))

### BASH SCRIPT FOR CHR 1-22 (as generated for my case) ###
# sleep and touch commands were added manually to allow a parallel loop with the R code #

sleep 20s
 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos1.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr1.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr1.tmp 
touch done1.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos2.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr2.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr2.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr2.tmp 
rm done1.log; touch done2.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos3.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr3.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr3.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr3.tmp 
rm done2.log; touch done3.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos4.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr4.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr4.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr4.tmp 
rm done3.log; touch done4.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos5.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr5.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr5.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr5.tmp 
rm done4.log; touch done5.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos6.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr6.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr6.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr6.tmp 
rm done5.log; touch done6.log
sleep 10m

Copy and paste this command ^ into your terminal and do not press enter until it has completed
 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos7.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr7.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr7.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr7.tmp 
rm done6.log; touch done7.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos8.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr8.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr8.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr8.tmp 
rm done7.log; touch done8.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos9.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr9.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr9.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr9.tmp 
rm done8.log; touch done9.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos10.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr10.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr10.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr10.tmp 
rm done9.log; touch done10.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos11.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr11.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr11.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr11.tmp 
rm done10.log; touch done11.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos12.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr12.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr12.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr12.tmp 
rm done11.log; touch done12.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos13.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr13.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr13.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr13.tmp 
rm done12.log; touch done13.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos14.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr14.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr14.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr14.tmp 
rm done13.log; touch done14.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos15.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr15.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr15.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr15.tmp 
rm done14.log; touch done15.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos16.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr16.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr16.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr16.tmp 
rm done15.log; touch done16.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos17.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr17.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr17.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr17.tmp 
rm done16.log; touch done17.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos18.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr18.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr18.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr18.tmp 
rm done17.log; touch done18.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos19.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr19.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr19.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr19.tmp 
rm done18.log; touch done19.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos20.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr20.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr20.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr20.tmp 
rm done19.log; touch done20.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos21.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr21.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr21.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr21.tmp 
rm done20.log; touch done21.log
sleep 10m

 /usr/bin/awk 'NR==FNR{pats[$0]; next} $2 in pats' ~/PLAY/gwasSnpPos22.txt <(paste -d ' ' <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr22.legend.gz | tail -n +2) <(zcat ~/barrett/1000GP_Phase3/1000GP_Phase3_chr22.hap.gz | sed 's/\(.\) \(.\)/\1\2/g' | sed 's/00/3/g' | sed 's/11/1/g' | sed 's/01/2/g' | sed 's/10/2/g')) > ~/PLAY/myChipChr22.tmp 
rm done21.log; touch done22.log
sleep 10m
