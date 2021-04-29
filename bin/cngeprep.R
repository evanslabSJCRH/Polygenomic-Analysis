invisible(options(echo = TRUE)) #Turn on echoing

ge.results <- read.table ("ge_lc50.tsv", header=TRUE)

head (ge.results)
dim (ge.results)

colnames (ge.results)

meta.alpha <- read.delim ("gecutoff.txt", header=FALSE)
meta.alpha <- unlist (meta.alpha)
names (meta.alpha) <- NULL


ge.results <- subset (ge.results, ge.results$p.b < meta.alpha)

#ge.results <- ge.results[1:10,]#remove me later after testing complete

ge.results <- ge.results[,1]

write.table (ge.results, file="ge_sig_probes_cn.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



