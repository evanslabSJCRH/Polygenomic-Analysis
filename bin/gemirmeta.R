invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig

protocol <- tolower(as.character(gsub("--", "", args.orig[11])))
mirprobe <- gsub("--", "", args.orig[11])

protocol
mirprobe

if(file.exists(paste("ge_mir_meta/", mirprobe, ".tsv", sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)

totxv.mir.results <- read.table(paste("totxv_ge_mir/", mirprobe, ".tsv", sep=""), header=TRUE, stringsAsFactors=FALSE)
totxvi.mir.results <- read.table(paste("totxvi_ge_mir/", mirprobe, ".tsv", sep=""), header=TRUE, stringsAsFactors=FALSE)
all.mir.results <- read.table(paste("all_ge_mir/", mirprobe, ".tsv", sep=""), header=TRUE, stringsAsFactors=FALSE)

mir.results <- merge (totxv.mir.results, totxvi.mir.results)
mir.results <- merge (mir.results, all.mir.results)

colnames (mir.results)


mir.results$meta.stat.b <- apply (mir.results[,c("totxv.lm.stat", "totxvi.lm.stat")], 1, function (x){sum(x)/sqrt(length(x))})
mir.results$meta.p.b <- 2*(pnorm(abs(mir.results$meta.stat.b),lower.tail=FALSE))
mir.results$p.b <- mir.results$all.lm.p
mir.results$stat.b <- mir.results$all.lm.stat

write.table (mir.results, file=paste("ge_mir_meta/", mirprobe, ".tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")




