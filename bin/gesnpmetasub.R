invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))

args.orig
idx.orig
probe <- gsub ("--", "", args.orig[11])
cprobe <- gsub ("/", "_", probe)


if(file.exists(paste("ge_snp_meta_sub/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

snp.result <- read.delim (paste ("ge_snp_meta/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))


head (snp.result)
#q(save="no")

snp.result <- subset (snp.result, !is.na(snp.result$cmh.meta) & snp.result$cmh.meta < 0.05)


dim (snp.result)
write.table (snp.result, file=paste("ge_snp_meta_sub/output.", cprobe, ".", sprintf ("%03.0f", idx.orig), sep=""), col.names=FALSE, row.names=TRUE, sep="\t")
