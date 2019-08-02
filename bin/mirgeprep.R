invisible(options(echo = TRUE)) #Turn on echoing

mir.results <- read.table ("mir_lc50.tsv", header=TRUE)


mir.cutoff <- read.delim ("mircutoff.txt", header=FALSE)
mir.cutoff <- unlist (mir.cutoff)
names (mir.cutoff) <- NULL

mir.results <- subset (mir.results, mir.results$p.b < mir.cutoff)

mir.results <- mir.results[,1]


mir.anno <- read.csv ("miranno.csv", as.is=TRUE, stringsAsFactors=FALSE)


mir.result <- merge (mir.results, mir.anno)
mir.result <- subset (mir.result, mir.result$Int == 1)



write.table (mir.results, file="mir_sig_probes.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")



