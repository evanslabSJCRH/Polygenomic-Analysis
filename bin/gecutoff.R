r.lib <- "drworkflow_Rlib"
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


all.result <- read.delim("ge_lc50.tsv", as.is=TRUE, stringsAsFactors=FALSE)





mtt.select <- subset (mtt.select, mtt.select$LC50.GROUP %in% c(1,3))

mtt.select <- subset (mtt.select, mtt.select$LIN == "B")
rownames (mtt.select) <- mtt.select$MRN
 
all.result.small <- subset (all.result, all.result$p.b < 0.05)
all.result.small <- all.result.small[order(all.result.small$p.b),]
all.result.small <- all.result.small[1:500,]

#if (drug %in% c("PRED", "LASP", "VCR", "6MP", "6TG")){
#  dutch.clust <- rep (NA, times=nrow(all.result.small))
#}

xv.clust <- rep (NA, times=nrow(all.result.small))
xvi.clust <- rep (NA, times=nrow(all.result.small))
dim (all.result.small)

for (i in 2:500){

probes.int.b <- all.result.small[1:i,"Probe.Set.ID"]



all.ge.data.b <- cbind(t(scale(t(stjude.dxbm.hm.mas5.probe.log2[probes.int.b,]))), t(scale(t(stjude.dxbm.xvi.mas5.probe.log2[probes.int.b,]))))



                                        
pt.int.b <- intersect(rownames(mtt.select), colnames(all.ge.data.b))

mtt.select <- mtt.select[pt.int.b, ]

all.ge.data.b <- all.ge.data.b[, pt.int.b]



xv.ge.data <- all.ge.data.b[,as.character(mtt.select[mtt.select$PROT == "TOTXV", "MRN"])]
xv.ge.pt.hc <- hcluster(t(xv.ge.data), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xv.clust[i] <- fisher.test(table(cutree(xv.ge.pt.hc, k=2), mtt.select[mtt.select$PROT == "TOTXV", "LC50.GROUP"]))$p.value



xvi.ge.data <- all.ge.data.b[,as.character(mtt.select[mtt.select$PROT == "TOTXVI", "MRN"])]
xvi.ge.pt.hc <- hcluster(t(xvi.ge.data), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xvi.clust[i] <- fisher.test(table(cutree(xvi.ge.pt.hc, k=2), mtt.select[mtt.select$PROT == "TOTXVI", "LC50.GROUP"]))$p.value

cat (i)
cat ("\n")
}


results <- cbind (xv.clust, xvi.clust, all.result.small$p.b)

results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)




colnames (results) <- c("xv.clust","xvi.clust","p.b")


results <- results[-1,]




results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2

results$order <- (1:nrow(results))+1

write.csv (results, file="gecutoffresult.csv")

results <- subset (results, results$order > 10)




results$max <- pmax (results$xv.clust, results$xvi.clust)


results.orig <- results
results <- results[order(results$meta),]

if (results[1,"max"] < 0.01){
  cutoff <- results.orig[min(which (results.orig$meta == min (results.orig$meta)))+1,"p.b"]
  myorder <- results.orig[min(which (results.orig$meta == min (results.orig$meta))),"order"]
  myminp <- results.orig[min(which (results.orig$meta == min (results.orig$meta))),"meta"]
  write.table (cutoff, "gecutoff.txt", col.names=FALSE, row.names=FALSE)

}

save.image("test.RData")

if (results[1,"max"] >= 0.01 && drug != "VCR"){
 results.orig <- results
 results <- subset(results, results$max < 0.01)

}



    results$exponent.sum <- floor(-log10(results$xv.clust))+floor(-log10(results$xvi.clust))


cutoff <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]+1==results.orig$order),"p.b"]  
myorder <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]==results.orig$order),"order"]  
myminp <- results.orig[which(results[min(which (results$exponent.sum == max (results$exponent.sum))),"order"]==results.orig$order),"meta"]
### for 0.0001 cutoff
cutoff <- 0.0001


  write.table (cutoff, "gecutoff.txt", col.names=FALSE, row.names=FALSE)




columns.int <- c("xv.clust", "xvi.clust", "p.b", "meta")

plims <- -log10(c(0.99, min (unlist(results.orig[,columns.int]))))
plims[2] <- plims[2]+5

pdf(file=paste(drug, "_gecutoff.pdf", sep=""), width=10, height=10)

x.pos <- c(0,1)/1
y.pos <- c(0,1)/1

plot (0,0,type="n", axes=FALSE, xlab="", ylab="")

par(fig=c(x.pos[1], x.pos[2], y.pos[1],y.pos[2]), new=T)

plot (results.orig$order, -log10(results.orig$xv.clust), type="l", ylim=plims, col=2, xlab="Number of Features", ylab="-log10(p-value)")
points (results.orig$order, -log10(results.orig$xvi.clust), type="l", col=3)
points (results.orig$order, -log10(results.orig$meta), type="l", col=4)
points (results.orig$order, -log10(results.orig$p.b), type="l", col=5)

mytext <- c("St. Jude Protocol XV", "St. Jude Protocol XVI", "Meta Clustering", "Alpha value")
mycols <- 2:5



arrows(myorder, -log10(myminp)+1.2, myorder, -log10(myminp)+0.2, length=0.1)
text(myorder, -log10(myminp)+1.2, paste("Meta Clustering p=", signif(myminp, 3),"\nAlpha=", signif(cutoff, 3), "\nFeatures=", myorder, sep=""), pos=4)


dev.off()


