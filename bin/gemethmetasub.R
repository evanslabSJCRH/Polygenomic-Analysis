invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))

args.orig
idx.orig
#q(save="no")
meth.result <- read.delim (paste ("ge_meth_meta/output.", sprintf("%03.0f", idx.orig), sep=""))


meth.result <- subset (meth.result, !is.na(meth.result$meta.p.b) & meth.result$meta.p.b < 0.05 | meth.result$all.lm.p <0.05)


dim (meth.result)
write.table (meth.result, file=paste("ge_meth_meta_sub/small.", sprintf ("%03.0f", idx.orig), sep=""), row.names=FALSE, sep="\t")
