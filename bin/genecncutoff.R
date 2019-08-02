rm(list=ls())

if(file.exists("cncutoff.txt")){
  q(save="no")
}

r.lib <- "drworkflow_Rlib"


require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)

require (amap, lib.loc=r.lib)

wd <- getwd()
wdsplit <- unlist (strsplit (wd, "/"))
drug <- wdsplit[length(wdsplit)]
drug
##############

cn.result.stat <- read.delim ("cn_lc50.tsv", as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),] 
cn.result.stat <- cn.result.stat[1:1000,]

load ("cn.data.som.seg.biostat.201607a.RData")
colnames(cn.data.som) <- gsub(".somatic","",colnames(cn.data.som))
head(cn.data.som)
head(cn.result.stat)

totxv.cn.data <- cn.data.som[rownames(cn.data.som) %in% cn.result.stat$Probe.Set.ID,]
head(totxv.cn.data)


head(totxv.cn.data)

mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")
drug
rownames(mtt.select) <- mtt.select$MRN
mtt.select

mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXV")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))



cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxv.cn.data))
cn.pt.int
totxv.cn.data <- totxv.cn.data[,as.character(cn.pt.int)]

head(totxv.cn.data)

mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
rownames (mtt.select.b) <- mtt.select.b$MRN
head(mtt.select.b)
mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


totxv.cn.data <- as.matrix (totxv.cn.data)
#totxv.cn.smooth.data <- as.matrix(totxv.cn.smooth.data)
dim(totxv.cn.data)

totxv.drug.cat <- mtt.select.b$LC50.GROUP
head(totxv.drug.cat)
head(totxv.cn.data)

load ("cn.data.som.seg.biostat.201607a.RData")

colnames(cn.data.som) <- gsub(".somatic","",colnames(cn.data.som))


 totxvi.cn.data <- cn.data.som[cn.result.stat[,1],]
 #totxvi.cn.smooth.data <- cn.data.som[cn.result.stat[,1],]
head(totxvi.cn.data)



mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")
 rownames(mtt.select) <- mtt.select$MRN
 mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
 mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXVI")
 mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))


 cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxvi.cn.data))
cn.pt.int
 totxvi.cn.data <- totxvi.cn.data[,as.character(cn.pt.int)]
 #totxvi.cn.smooth.data <- totxvi.cn.smooth.data[,as.character(cn.pt.int)]

 mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
 rownames (mtt.select.b) <- mtt.select.b$MRN

 mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


 totxvi.cn.data <- as.matrix (totxvi.cn.data)

 dim (totxvi.cn.data)
head(totxvi.cn.data)

totxvi.drug.cat <- mtt.select.b$LC50.GROUP

###################################

xv.cn.clust.p <- rep (NA, times=nrow(cn.result.stat))
xvi.cn.clust.p <- rep (NA, times=nrow(cn.result.stat))

xv.cn.clust.cat.p <- rep (NA, times=nrow(cn.result.stat))
xvi.cn.clust.cat.p <- rep (NA, times=nrow(cn.result.stat))

clut.methods <- "correlation"

link.methods <- "ward"
                                        #


for (k in 1:length(link.methods)){
for (j in 1:length(clut.methods)){

for (i in 25:nrow(cn.result.stat)){
  cn.probes.int <- cn.result.stat[1:i,"Probe.Set.ID"]

 totxv.cn <- totxv.cn.data[cn.probes.int,]
 
  totxv.cn.hc2 <- hcluster (t(totxv.cn), method =clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

  totxv.cn.clust.p <- fisher.test(table(cutree(totxv.cn.hc2, k=2), totxv.drug.cat))$p.value
  totxv.cn.clust.p

  xv.cn.clust.p[i] <- totxv.cn.clust.p


 

  
 # totxvi.cn <- totxvi.cn.smooth.data[cn.probes.int,]
 totxvi.cn <- totxvi.cn.data[cn.probes.int,]

  totxvi.cn.hc2 <- hcluster (t(totxvi.cn), method =clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

  totxvi.cn.clust.p <- fisher.test(table(cutree(totxvi.cn.hc2, k=2), totxvi.drug.cat))$p.value
  totxvi.cn.clust.p

  xvi.cn.clust.p[i] <- totxvi.cn.clust.p

}
cat(clut.methods[j])
cat(link.methods[k])
cat(min(xv.cn.clust.p, na.rm=TRUE))
cat(min(xvi.cn.clust.p, na.rm=TRUE))
cat ("\n")
}

}

##############

#q()


 results <- cbind (xv.cn.clust.p, xvi.cn.clust.p, cn.result.stat$p.b)
results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)

colnames (results) <- c("xv.clust", "xvi.clust", "p.b")

results$better <- rep(NA, times=nrow(results))
for(i in 25:nrow(results)){

  results[i,"better"] <- dim (results[results[25:nrow(results),"xv.clust"] <= results[i,"xv.clust"] & results[25:nrow(results),"xvi.clust"] <= results[i,"xvi.clust"],])[1]

}



 results <- results[-1,]

 results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2
 results$order <- (1:nrow(results))+1


 write.csv (results, file="cncutoffresult.csv")

 results <- subset(results, results$order >= 25)

for (i in 2:nrow(results)){
results[i,"same"] <- ifelse (results[i,"meta"] == results[i+1,"meta"], TRUE, FALSE)
results[i,"same.p.b"] <- ifelse (results[i,"p.b"] == results[i+1,"p.b"], TRUE, FALSE)
}

threshold <- results[min(which (results$meta == min (results$meta))),"order"]

changepoint <- results[which(results$same == FALSE & results$same.p.b == FALSE),"order"]

changepoint <- subset (changepoint, changepoint > threshold)

changepoint <- min (changepoint)
 write.csv (results, file="cncutoffresult_test.csv")

write.table (results[which (results$better == min(results$better, na.rm=TRUE)),"p.b"], "cncutoff.txt", col.names=FALSE, row.names=FALSE)

write.table (rownames(results)[which (results$better == min(results$better, na.rm=TRUE))], "cncutoffindex.txt", col.names=FALSE, row.names=FALSE)



