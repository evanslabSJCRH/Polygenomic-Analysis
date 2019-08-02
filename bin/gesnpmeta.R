invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig



#drug <- tolower(as.character(gsub("--", "", args.orig[8])))
idx.orig <- gsub("--", "", args.orig[9])
idx.orig <- as.numeric(idx.orig)

probe <- gsub ("--", "", args.orig[11])

cprobe <- gsub("/", "_", probe)

#protocol
idx.orig

if(file.exists(paste("ge_snp_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
#require (EvansData, lib.loc=r.lib)
#require (EvansAnalysis, lib.loc=r.lib)
#require (SJHMGEData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (VGAM, lib.loc=r.lib)
require (lawstat, lib.loc=r.lib)



r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
#require (EvansData, lib.loc=r.lib)
#require (EvansAnalysis, lib.loc=r.lib)
#require (SJHMGEData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (VGAM, lib.loc=r.lib)
require (lawstat, lib.loc=r.lib)




totxv.snp.ge.results <- read.table(paste("totxv_ge_snp/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)
totxvi.snp.ge.results <- read.table(paste("totxvi_ge_snp/output.",cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)
all.snp.ge.results <- read.table(paste("all_ge_snp/output.",cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)

#colnames (totxv.snp.ge.results) <- c("SNP.Probe.Set.ID", "GE.Probe.Set.ID", "totxv.fisher.p", "

head (totxv.snp.ge.results)
head (totxvi.snp.ge.results)

#q(save="no")

snp.ge.result <- cbind (totxv.snp.ge.results, totxvi.snp.ge.results)
snp.ge.result <- cbind (snp.ge.result,all.snp.ge.results)

head(snp.ge.result)
cmh.meta <- rep(NA, times=nrow(snp.ge.result))

for (i in 1:nrow(snp.ge.result)){


  if (!is.na(snp.ge.result[i,9])&!is.na(snp.ge.result[i,2])){

    totxv.snp.matrix <- matrix(c(snp.ge.result[i,6], snp.ge.result[i,7], snp.ge.result[i,8], snp.ge.result[i,9], snp.ge.result[i,10], snp.ge.result[i,11]), nrow=2)
    totxvi.snp.matrix <- matrix(c(snp.ge.result[i,15], snp.ge.result[i,16], snp.ge.result[i,17], snp.ge.result[i,18], snp.ge.result[i,19], snp.ge.result[i,20]), nrow=2)

    snparray <- array(c(totxv.snp.matrix,totxvi.snp.matrix), dim=c(2,3,2))

    if (length (which(apply (rbind (snparray[,,1], snparray[,,2]), 2, function (x){sum(x==0)}) != 4))>1){
      if (length (which (apply(snparray[,,1], 2, function (x){sum(x==0)}) != 2))>1&length (which (apply(snparray[,,2], 2, function (x){sum(x==0)}) != 2))>1){
        if (length (which (apply(snparray[,,1], 1, function (x){sum(x==0)})!= 3 ))>1&length (which (apply(snparray[,,2], 1, function (x){sum(x==0)}) != 3))>1){
          snparray <- snparray [,which(apply (rbind (snparray[,,1], snparray[,,2]), 2, function (x){sum(x==0)}) != 4),]
          if (all(snparray !=0)){
          cmh.meta[i] <-  mantelhaen.test (snparray)$p.value
        }else{
      cmh.meta[i] <- NA
    }
}
      }
    }
  }
}

snp.ge.result <- cbind(snp.ge.result,cmh.meta)
snp.ge.result$meta.stat.b <- apply (snp.ge.result[,c("totxv.lm.stat.b", "totxvi.lm.stat.b")], 1, function (x){sum(x)/sqrt(length(x))})
snp.ge.result$meta.p.b <- 2*(pnorm(abs(snp.ge.result$meta.stat.b),lower.tail=FALSE))


write.table (snp.ge.result, file=paste("ge_snp_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")
