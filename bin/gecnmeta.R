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

if(file.exists(paste("ge_cn_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
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




totxv.cn.ge.results <- read.table(paste("totxv_ge_cn/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=FALSE, stringsAsFactors=FALSE)
totxvi.cn.ge.results <- read.table(paste("totxvi_ge_cn/output.",cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=FALSE, stringsAsFactors=FALSE)
all.cn.ge.results <- read.table(paste("all_ge_cn/output.",cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), header=FALSE, stringsAsFactors=FALSE)

#colnames (totxv.snp.ge.results) <- c("SNP.Probe.Set.ID", "GE.Probe.Set.ID", "totxv.fisher.p", "

head (totxv.cn.ge.results)
head (totxvi.cn.ge.results)

#q(save="no")

cn.ge.result <- cbind (totxv.cn.ge.results, totxvi.cn.ge.results[,3:ncol(totxvi.cn.ge.results)])
cn.ge.result <- cbind (cn.ge.result, all.cn.ge.results[,3:ncol(all.cn.ge.results)])

head (cn.ge.result)
#q(save="no")

 colnames (cn.ge.result) <- c("SNPProbeSetID", "GEProbeSetID", "totxv.t.p.b", "totxv.t.stat.b", "totxv.w.p.b", "totxv.w.stat.b","totxv.lm.p.b", "totxv.lm.stat.b","totxvi.t.p.b", "totxvi.t.stat.b", "totxvi.w.p.b", "totxvi.w.stat.b","totxvi.lm.p.b", "totxvi.lm.stat.b","all.t.p.b","all.t.stat.b","all.w.p.b","all.w.stat.b","all.lm.stat.b","all.lm.p.b")



cn.ge.result$meta.stat.b <- apply (cn.ge.result[,c("totxv.lm.stat.b", "totxvi.lm.stat.b")], 1, function (x){sum(x)/sqrt(length(x))})
cn.ge.result$meta.p.b <- 2*(pnorm(abs(cn.ge.result$meta.stat.b),lower.tail=FALSE))


write.table (cn.ge.result, file=paste("ge_cn_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")
