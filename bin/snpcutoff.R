rm(list=ls())

#if(file.exists("snpcutoff.txt")){
#  q(save="no")
#}


r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
#"C:/Documents and Settings/rautry/My Documents/R/win-library/2.12"

require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)
require (SJHMSNPData, lib.loc=r.lib)
require (amap, lib.loc=r.lib)

wd <- getwd()
wdsplit <- unlist (strsplit (wd, "/"))
drug <- wdsplit[length(wdsplit)]
drug

##############


mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")

rownames(mtt.select) <- mtt.select$MRN
mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))


data ("20111216totxvsnpdatasom")

totxv.snp.data <- snp.data.totxv.som

totxv.snp.pt.int <- intersect (mtt.select.b$MRN, colnames (totxv.snp.data))

totxv.snp.data <- snp.data.totxv.som[,as.character(totxv.snp.pt.int)]

mtt.select.totxv.b <- mtt.select.b[as.character(totxv.snp.pt.int),]

totxv.snp.data[is.na(totxv.snp.data)] <- NA
totxv.snp.data[totxv.snp.data == "NN"] <- NA
totxv.snp.data[totxv.snp.data == "NC"] <- NA
totxv.snp.data[totxv.snp.data == "AA"] <- 0
totxv.snp.data[totxv.snp.data == "AB"] <- 1
totxv.snp.data[totxv.snp.data == "BB"] <- 2

totxv.snp.data <- as.data.frame(totxv.snp.data, as.is=TRUE, stringsAsFactors=FALSE)

for (i in 1:ncol(totxv.snp.data)){
totxv.snp.data[,i] <- as.numeric(totxv.snp.data[,i])
}

totxv.snp.data <- as.matrix (totxv.snp.data)


######################

data ("20111216totxvisnpdatasom")

totxvi.snp.data <- snp.data.totxvi.som

totxvi.snp.pt.int <- intersect (mtt.select.b$MRN, colnames (totxvi.snp.data))

totxvi.snp.data <- snp.data.totxvi.som[,as.character(totxvi.snp.pt.int)]

mtt.select.totxvi.b <- mtt.select.b[as.character(totxvi.snp.pt.int),]

totxvi.snp.data[is.na(totxvi.snp.data)] <- NA
totxvi.snp.data[totxvi.snp.data == "NN"] <- NA
totxvi.snp.data[totxvi.snp.data == "NC"] <- NA
totxvi.snp.data[totxvi.snp.data == "AA"] <- 0
totxvi.snp.data[totxvi.snp.data == "AB"] <- 1
totxvi.snp.data[totxvi.snp.data == "BB"] <- 2

totxvi.snp.data <- as.data.frame(totxvi.snp.data, as.is=TRUE, stringsAsFactors=FALSE)

for (i in 1:ncol(totxvi.snp.data)){
totxvi.snp.data[,i] <- as.numeric(totxvi.snp.data[,i])
}

totxvi.snp.data <- as.matrix (totxvi.snp.data)

############################
snp.result <- read.delim("snp_lc50.tsv", as.is=TRUE, stringsAsFactors=FALSE, header=TRUE)

snp.result <- snp.result[order (snp.result$p.b),]
snp.result <- snp.result[!is.na(snp.result$p.b),]


snp.result$xv.01 <- apply (snp.result[,c("totxv.fisher.01.b", "totxv.fisher.11.b", "totxv.fisher.21.b")], 1, sum)

snp.result$xv.03 <- apply (snp.result[,c("totxv.fisher.03.b", "totxv.fisher.13.b", "totxv.fisher.23.b")], 1, sum)

snp.result$xvi.01 <- apply (snp.result[,c("totxvi.fisher.01.b", "totxvi.fisher.11.b", "totxvi.fisher.21.b")], 1, sum)

snp.result$xvi.03 <- apply (snp.result[,c("totxvi.fisher.03.b", "totxvi.fisher.13.b", "totxvi.fisher.23.b")], 1, sum)

snp.result$min.count <- apply (snp.result[,c("xv.01", "xv.03", "xvi.01", "xvi.03")], 1, min)

snp.result <- snp.result[snp.result$min.count > 0,]

snp.result$xv.01 <- NULL
snp.result$xv.03 <- NULL
snp.result$xvi.01 <- NULL
snp.result$xvi.03 <- NULL
snp.result$min.count <- NULL

snp.result <- snp.result[1:500,]

xv.snp.clust.p <- rep (NA, times=nrow(snp.result))
xvi.snp.clust.p <- rep (NA, times=nrow(snp.result))

#save.image ("trouble.RData")

for (i in 10:nrow(snp.result)){
  snp.probes.int <- snp.result[1:i,"ProbeSetID"]

  xv.snp.data <- totxv.snp.data[snp.probes.int,]

  totxv.snp.hc2 <- hcluster (t(xv.snp.data), method = "bin", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

  totxv.snp.clust.p <- fisher.test(table(cutree(totxv.snp.hc2, k=2), mtt.select.totxv.b$LC50.GROUP))$p.value
  totxv.snp.clust.p

  xv.snp.clust.p[i] <- totxv.snp.clust.p
  
  xvi.snp.data <- totxvi.snp.data[snp.probes.int,]

  totxvi.snp.hc2 <- hcluster (t(xvi.snp.data), method = "bin", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

  totxvi.snp.clust.p <- fisher.test(table(cutree(totxvi.snp.hc2, k=2), mtt.select.totxvi.b$LC50.GROUP))$p.value
  totxvi.snp.clust.p

  xvi.snp.clust.p[i] <- totxvi.snp.clust.p 

cat (i)
cat ("\n")
}

##############




results <- cbind (xv.snp.clust.p, xvi.snp.clust.p, snp.result$p.b)

colnames (results) <- c("xv.clust", "xvi.clust", "p.b")

results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)

results <- results[-1,]

results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2
results$order <- (1:nrow(results))+1


write.csv (results, file="snpcutoffresult.csv")


results <- subset(results, results$order >= 10)
    results$exponent.sum <- floor(-log10(results$xv.clust))+floor(-log10(results$xvi.clust))
results.orig <- results

write.table (results[min(which (results$meta == min (results$meta)))+1,"p.b"], "snpcutoff.txt", col.names=FALSE, row.names=FALSE)

save.image("snpcutoffdata.RData")
columns.int <- c("xv.clust", "xvi.clust", "p.b", "meta")
plims <- -log10(c(0.99, min (unlist(results.orig[,columns.int]))))
plims[2] <- plims[2]+2

cutoff <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]+1==results.orig$order),"p.b"]  
myorder <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]==results.orig$order),"order"]  
myminp <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]==results.orig$order),"meta"]  

pdf(file=paste(drug, "_snpcutoff.pdf", sep=""))
plot (results.orig$order, -log10(results.orig$xv.clust), type="l", ylim=plims, col=2, xlab="Number of Features", ylab="-log10(p-value)")
points (results.orig$order, -log10(results.orig$xvi.clust), type="l", col=3)
points (results.orig$order, -log10(results.orig$meta), type="l", col=4)
points (results.orig$order, -log10(results.orig$p.b), type="l", col=5)

mytext <- c("St. Jude Protocol XV", "St. Jude Protocol XVI", "Meta Clustering", "Alpha value")
mycols <- 2:5

legend(500, plims[2], mytext, col=mycols, pch=16)
arrows(myorder, -log10(myminp)+1.2, myorder, -log10(myminp)+0.2, length=0.1)
text(myorder, -log10(myminp)+1.2, paste("Meta Clustering p=", signif(myminp, 3),"\nAlpha=", signif(cutoff, 3), "\nFeatures=", myorder, sep=""), pos=4)
dev.off()
