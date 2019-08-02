invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig

#protocol <- tolower(as.character(gsub("--", "", args.orig[11])))
idx.orig <- gsub("--", "", args.orig[9])
idx.orig <- as.numeric(idx.orig)
#protocol
idx.orig
#q(save="no")
if(file.exists(paste("ge_meth_meta/output.", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)

totxv.meth.results <- read.table(paste("totxv_ge_meth/output.", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)
totxvi.meth.results <- read.table(paste("totxvi_ge_meth/output.", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)
all.meth.results <- read.table(paste("all_ge_meth/output.", sprintf("%03.0f", idx.orig), sep=""), header=TRUE, stringsAsFactors=FALSE)

meth.results <- merge (totxv.meth.results, totxvi.meth.results)
meth.results <- merge(meth.results,all.meth.results)
colnames (meth.results)


meth.results$meta.stat.b <- apply (meth.results[,c("totxv.lm.stat", "totxvi.lm.stat")], 1, function (x){sum(x)/sqrt(length(x))})
meth.results$meta.p.b <- 2*(pnorm(abs(meth.results$meta.stat.b),lower.tail=FALSE))


write.table (meth.results, file=paste("ge_meth_meta/output.", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")
