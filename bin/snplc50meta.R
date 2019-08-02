invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))
drug <- as.character(gsub("--", "", args.orig[10])) 
protocol <- tolower(as.character(gsub("--", "", args.orig[11])))

split.num <- 200


if(file.exists(paste("snp_lc50_meta/output.",sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"

require (mvtnorm, lib.loc=r.lib)
require (VGAM, lib.loc=r.lib)
require (lawstat, lib.loc=r.lib)

getwd()

snp.result <- read.delim("snp_result.tsv")


ntests <- nrow(snp.result)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1  
end.row <-  min (block.size*idx.orig, ntests)
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row

snp.result <- snp.result[begin.row:end.row,]




##################

cmh.meta <- rep(NA, times=nrow(snp.result))

for (i in 1:nrow(snp.result)){


  if (!is.na(snp.result[i,"totxvi.fisher.p.b"])&!is.na(snp.result[i,"totxv.fisher.p.b"])){
    
    totxv.snp.matrix <- matrix(c(snp.result[i,"totxv.fisher.01.b"], snp.result[i,"totxv.fisher.03.b"], snp.result[i,"totxv.fisher.11.b"], snp.result[i,"totxv.fisher.13.b"], snp.result[i,"totxv.fisher.21.b"], snp.result[i,"totxv.fisher.23.b"]), nrow=2)
    
    totxvi.snp.matrix <- matrix(c(snp.result[i,"totxvi.fisher.01.b"], snp.result[i,"totxvi.fisher.03.b"], snp.result[i,"totxvi.fisher.11.b"], snp.result[i,"totxvi.fisher.13.b"], snp.result[i,"totxvi.fisher.21.b"], snp.result[i,"totxvi.fisher.23.b"]), nrow=2)

    
    snparray <- array(c(totxv.snp.matrix,totxvi.snp.matrix), dim=c(2,3,2))
    
    if (length (which(apply (rbind (snparray[,,1], snparray[,,2]), 2, function (x){sum(x==0)}) != 4))>1){
      if (length (which (apply(snparray[,,1], 2, function (x){sum(x==0)}) != 2))>1&length (which (apply(snparray[,,2], 2, function (x){sum(x==0)}) != 2))>1){
        if (length (which (apply(snparray[,,1], 1, function (x){sum(x==0)})!= 3 ))>1&length (which (apply(snparray[,,2], 1, function (x){sum(x==0)}) != 3))>1){
          snparray <- snparray [,which(apply (rbind (snparray[,,1], snparray[,,2]), 2, function (x){sum(x==0)}) != 4),]
          cmh.meta[i] <-  mantelhaen.test (snparray)$p.value
        }
      }
    }
  }
}




snp.result <- cbind (snp.result, cmh.meta)
##################

snp.result$p.b <- snp.result$all.snp.lm.p.b
snp.result$stat.b <- snp.result$all.snp.lm.stat.b


write.table (snp.result, file=paste("snp_lc50_meta/output.",sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")



