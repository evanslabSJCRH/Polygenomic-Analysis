invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

drug <- as.character(gsub("--", "", args.orig[9])) 

drug

args.orig


if(file.exists("snp_result.tsv")){
  q(save="no")
}

totxv.result <- read.delim("totxv_snp_lc50.tsv")

rownames (totxv.result) <- totxv.result$ProbeSetID
totxv.result$ProbeSetID <- NULL

colnames (totxv.result) <- paste ("totxv.", colnames (totxv.result), sep="")


totxvi.result <- read.delim("totxvi_snp_lc50.tsv")

rownames (totxvi.result) <- totxvi.result$ProbeSetID
totxvi.result$ProbeSetID <- NULL

colnames (totxvi.result) <- paste ("totxvi.", colnames (totxvi.result), sep="")


all.result <- read.delim("all_snp_lc50.tsv")

rownames (all.result) <- all.result$ProbeSetID
all.result$ProbeSetID <- NULL

colnames (all.result) <- paste ("all.", colnames (all.result), sep="")


probe.int <- intersect (rownames(totxv.result), rownames (totxvi.result))

probe.int <- intersect (probe.int, rownames(all.result))

totxv.result <- totxv.result[probe.int,]
totxvi.result <- totxvi.result[probe.int,]
all.result <- all.result[probe.int,]

snp.result <- cbind (totxv.result, totxvi.result, all.result)

snp.result$ProbeSetID <- rownames (snp.result)

snp.result <- snp.result[,c(ncol(snp.result),1:(ncol(snp.result)-1))]

write.table (snp.result, file="snp_result.tsv", row.names=FALSE, quote=FALSE, sep="\t")


