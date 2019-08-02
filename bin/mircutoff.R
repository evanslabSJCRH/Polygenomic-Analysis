rm(list=ls())

if(file.exists("mircutoff.txt")){
  q(save="no")
}


r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
#"C:/Documents and Settings/rautry/My Documents/R/win-library/2.12"

require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)
require (amap, lib.loc=r.lib)

wd <- getwd()
wdsplit <- unlist (strsplit (wd, "/"))
drug <- wdsplit[length(wdsplit)]
drug

##############

mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")

totxv.mir.mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")



totxv.mir.mtt.select.b <- subset (totxv.mir.mtt.select.b, totxv.mir.mtt.select.b$LC50.GROUP %in% c(1,3))
rownames(totxv.mir.mtt.select.b) <- totxv.mir.mtt.select.b$MRN

totxv.mir.select.b <- intersect (colnames (mir.data.totxv.som), totxv.mir.mtt.select.b$MRN)

mir.data.totxv.som <- t(scale(t(mir.data.totxv.som)))

totxv.mirdata.b <- mir.data.totxv.som[,as.character(totxv.mir.select.b)]
totxv.mir.mtt.select.b <- totxv.mir.mtt.select.b[as.character(totxv.mir.select.b),]


totxvi.mir.mtt.select.b <- mtt.prep(drug)
source("drugcatadjust.R")
totxvi.mir.mtt.select.b <- subset (totxvi.mir.mtt.select.b, totxvi.mir.mtt.select.b$LIN == "B")
totxvi.mir.mtt.select.b <- subset (totxvi.mir.mtt.select.b, totxvi.mir.mtt.select.b$LC50.GROUP %in% c(1,3))
rownames(totxvi.mir.mtt.select.b) <- totxvi.mir.mtt.select.b$MRN

totxvi.mir.select.b <- intersect (colnames (mir.data.totxvi.som), totxvi.mir.mtt.select.b$MRN)

mir.data.totxvi.som <- t(scale(t(mir.data.totxvi.som)))

totxvi.mirdata.b <- mir.data.totxvi.som[,as.character(totxvi.mir.select.b)]
totxvi.mir.mtt.select.b <- totxvi.mir.mtt.select.b[as.character(totxvi.mir.select.b),]


mir.meta.result <- read.delim ("mir_lc50.tsv", as.is=TRUE, stringsAsFactors=FALSE)

mir.anno <- read.csv ("miranno.csv", as.is=TRUE, stringsAsFactors=FALSE)


mir.meta.result <- merge (mir.meta.result, mir.anno)
mir.meta.result <- subset (mir.meta.result, mir.meta.result$Int == 1)

mir.meta.result <- mir.meta.result[order(mir.meta.result$p.b),]

mir.meta.result <- subset (mir.meta.result, mir.meta.result$p.b < 0.05)

#if (drug == "PRED"){
#mir.meta.result <- subset (mir.meta.result, mir.meta.result$p.b < 0.05)
#}


##########


xv.mir.clust.p <- rep (NA, times=nrow(mir.meta.result))
xvi.mir.clust.p <- rep (NA, times=nrow(mir.meta.result))




for (i in 2:nrow(mir.meta.result)){
  mirdata.heat.b <- totxv.mirdata.b
  mir.probes.int.b <- mir.meta.result[1:i,"Name"]
  mirdata.heat.b <- mirdata.heat.b[as.character(mir.probes.int.b),]

#mirdata.heat.b <- t(mirdata.heat.b)
#mirdata.heat.b <- scale (mirdata.heat.b)
#mirdata.heat.b <- t (mirdata.heat.b)


totxv.mir.hc1 <- hcluster(mirdata.heat.b, method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxv.mir.hc2 <- hcluster(t(mirdata.heat.b), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xv.mir.clust.p[i] <- fisher.test(table(cutree(totxv.mir.hc2, k=2), totxv.mir.mtt.select.b$LC50.GROUP))$p.value


mirdata.heat.b <- totxvi.mirdata.b


mirdata.heat.b <- mirdata.heat.b[as.character(mir.probes.int.b),]

#mirdata.heat.b <- t(mirdata.heat.b)
#mirdata.heat.b <- scale (mirdata.heat.b)
#mirdata.heat.b <- t (mirdata.heat.b)

totxvi.mir.hc1 <- hcluster(mirdata.heat.b, method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxvi.mir.hc2 <- hcluster(t(mirdata.heat.b), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xvi.mir.clust.p[i] <- fisher.test(table(cutree(totxvi.mir.hc2, k=2), totxvi.mir.mtt.select.b$LC50.GROUP))$p.value

cat (i)
cat ("\n")
}

##############




results <- cbind (mir.meta.result$Name, xv.mir.clust.p, xvi.mir.clust.p, mir.meta.result$p.b)

colnames (results) <- c("Name", "xv.clust", "xvi.clust", "p.b")

results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)

results <- results[-1,]

results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2
results$order <- (1:nrow(results))+1


write.csv (results, file="mircutoffresult.csv")

results <- subset(results, results$order >= 10)

cutoff <- results[min(which (results$meta == min (results$meta)))+1,"p.b"]
mirprobes <- mir.meta.result[mir.meta.result$p.b < cutoff,"Name"]
myorder <- results[min(which (results$meta == min (results$meta))),"order"]
myminp <- results[min(which (results$meta == min (results$meta))),"meta"]

write.table (mirprobes, file="mir_sig_probes.tsv", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table (cutoff, "mircutoff.txt", col.names=FALSE, row.names=FALSE)



columns.int <- c("xv.clust", "xvi.clust", "p.b", "meta")

plims <- -log10(c(0.99, min (unlist(results[,columns.int]))))
plims[2] <- plims[2]+2

pdf(file=paste(drug, "_mircutoff.pdf", sep=""))
plot (results$order, -log10(results$xv.clust), type="l", ylim=plims, col=2, xlab="Number of Features", ylab="-log10(p-value)")
points (results$order, -log10(results$xvi.clust), type="l", col=3)
points (results$order, -log10(results$meta), type="l", col=4)
points (results$order, -log10(results$p.b), type="l", col=5)

mytext <- c("St. Jude Protocol XV", "St. Jude Protocol XVI", "Meta Clustering", "Alpha value")
mycols <- 2:5

mylegendx <- ifelse((dim(results)[1])/2 < myorder, 10, (dim(results)[1])/2)
legend(mylegendx, plims[2], mytext, col=mycols, pch=16)
arrows(myorder, -log10(myminp)+1.2, myorder, -log10(myminp)+0.2, length=0.1)
text(myorder, -log10(myminp)+1.2, paste("Meta Clustering p=", signif(myminp, 3),"\nAlpha=", signif(cutoff, 3), "\nFeatures=", myorder, sep=""), pos=4, cex=0.7)
dev.off()


