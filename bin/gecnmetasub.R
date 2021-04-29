invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))

args.orig
idx.orig
probe <- gsub ("--", "", args.orig[11])
cprobe <- gsub ("/", "_", probe)


if(file.exists(paste("ge_cn_meta_sub/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

cn.result <- read.delim (paste ("ge_cn_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))


head (cn.result)
#q(save="no")

#cn.result <- subset (cn.result, !is.na(cn.result$meta.p.b) & cn.result$meta.p.b < 0.05)
cn.result <- subset (cn.result, !is.na(cn.result$meta.p.b) & cn.result$all.lm.p.b < 0.05)


dim (cn.result)
write.table (cn.result, file=paste("ge_cn_meta_sub/small.", cprobe, ".", sprintf ("%03.0f", idx.orig), sep=""), col.names=FALSE, row.names=FALSE, sep="\t")
