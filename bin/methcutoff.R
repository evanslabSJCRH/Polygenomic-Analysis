rm(list=ls())

if(file.exists("methcutoff.txt")){
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


mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")

all.result <- read.delim("meth_lc50.tsv", as.is=TRUE, stringsAsFactors=FALSE)

mtt.select <- subset (mtt.select, mtt.select$LC50.GROUP %in% c(1,3))

mtt.select <- subset (mtt.select, mtt.select$LIN == "B")
rownames (mtt.select) <- mtt.select$MRN
 
all.result.small <- subset (all.result, all.result$p.b < 0.05)
all.result.small <- all.result.small[order(all.result.small$p.b),]
all.result.small <- all.result.small[1:500,]

xv.clust <- rep (NA, times=nrow(all.result.small))
xvi.clust <- rep (NA, times=nrow(all.result.small))
dim (all.result.small)

for (i in 2:500){

probes.int.b <- all.result.small[1:i,"IlmnID"]
                                       
pt.int.xv.b <- intersect(rownames(mtt.select), colnames(meth.data.totxv.som))

mtt.select.xv <- mtt.select[pt.int.xv.b, ]

totxv.meth.data.b <- meth.data.totxv.som[probes.int.b,as.character(pt.int.xv.b)]

xv.meth.pt.hc <- hcluster(t(totxv.meth.data.b), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xv.clust[i] <- fisher.test(table(cutree(xv.meth.pt.hc, k=2), mtt.select[pt.int.xv.b, "LC50.GROUP"]))$p.value




pt.int.xvi.b <- intersect(rownames(mtt.select), colnames(meth.data.totxvi.som))

mtt.select.xvi <- mtt.select[pt.int.xvi.b, ]

totxvi.meth.data.b <- meth.data.totxvi.som[probes.int.b,as.character(pt.int.xvi.b)]

xvi.meth.pt.hc <- hcluster(t(totxvi.meth.data.b), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xvi.clust[i] <- fisher.test(table(cutree(xvi.meth.pt.hc, k=2), mtt.select[pt.int.xvi.b, "LC50.GROUP"]))$p.value



cat (i)
cat ("\n")
}


results <- cbind (xv.clust, xvi.clust, all.result.small$p.b)


results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)


colnames (results) <- c("xv.clust","xvi.clust","p.b")



results <- results[-1,]

results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2

results$order <- (1:nrow(results))+1

write.csv (results, file="methcutoffresult.csv")

results <- subset (results, results$order >= 10)
results[min(which (results$meta == min (results$meta)))+1,"p.b"]
write.table (results[min(which (results$meta == min (results$meta)))+1,"p.b"], "methcutoff.txt", col.names=FALSE, row.names=FALSE)




