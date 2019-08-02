invisible(options(echo = TRUE)) #Turn on echoing

meth.results <- read.table ("meth_lc50.tsv", header=TRUE)

head (meth.results)
dim (meth.results)

colnames (meth.results)

meth.results <- subset (meth.results, meth.results$p.b < 0.05)

meth.results <- meth.results[,1]

write.table (meth.results, file="meth_sig_probes.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



