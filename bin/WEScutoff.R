#rm(list=ls())

#if(file.exists("cncutoff.txt")){
#  q(save="no")
#}

#r.lib <- "drworkflow_Rlib"

require (EvansData)
require (EvansAnalysis)

require (amap)
require(dplyr)
require("stringr")
wd <- getwd()
wdsplit <- unlist (strsplit (wd, "/"))
drug <- wdsplit[length(wdsplit)]
drug <- "PREDANDDEX"
setwd("/home/rautry/workflow_results_final/PRED/")
date <- Sys.Date()
##############
cn.result.stat <- read.csv("all.lm.b.mutationPREDANDDEXwithprotocol_2018-07-31.csv", as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)
colnames(cn.result.stat) <- c("Probe.Set.ID","stat.b","p.b")
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),] 
cn.result.stat <- cn.result.stat[1:500,]
head(cn.result.stat)

cn.data.som.pre <- read.csv("/home/rautry/data/mutation_WES/mutation.call_gene-ALLnoncoding2018-07-30-.csv", header = TRUE, check.names = FALSE, row.names =1)
head(cn.data.som.pre)

cn.data.som <- cn.data.som.pre %>% 
  mutate_all(funs(str_replace(., "YES", "1"))) %>% mutate_all(funs(str_replace(., "NO", "2")))
rownames(cn.data.som) <- rownames(cn.data.som.pre)
write.csv(cn.data.som, file = paste("WES.data.som_replacenumeric",date,".csv", sep =""), quote = FALSE , row.names = TRUE)
head(cn.data.som)
head(cn.result.stat)

totxv.cn.data <- cn.data.som[rownames(cn.data.som) %in% cn.result.stat$Probe.Set.ID,]
head(totxv.cn.data)


drug <- "PREDANDDEX"
#head(totxv.cn.smooth.data)
mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")
if(drug == "PREDANDDEX"){
  mtt.select.PRED <- mtt.prep("PRED")
  mtt.select.PRED$LC50.GROUP <- ifelse (mtt.select.PRED$LC50 < 0.1, 1, ifelse(mtt.select.PRED$LC50 >= 0.1 & mtt.select.PRED$LC50 < 65, 2, ifelse(mtt.select.PRED$LC50 >= 64, 3, NA)))
  mtt.select.PRED <- subset(mtt.select.PRED, mtt.select.PRED$LIN == "B")
  
  mtt.select.DEX <- mtt.prep("DEX")
  mtt.select.DEX <- subset(mtt.select.DEX, mtt.select.DEX$LIN == "B")
  
  mtt.select.DEX$LC50 <- (mtt.select.DEX$LC50) *8
  mtt.select.DEX.select <- mtt.select.DEX[!(mtt.select.DEX$MRN %in% mtt.select.PRED$MRN),]
  mtt.select <- rbind(mtt.select.PRED, mtt.select.DEX.select)
}

rownames(mtt.select) <- mtt.select$MRN
mtt.select

mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXV")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))



cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxv.cn.data))
cn.pt.int
totxv.cn.data <- totxv.cn.data[,as.character(cn.pt.int)]
#totxv.cn.smooth.data <- totxv.cn.smooth.data[,as.character(cn.pt.int)]
head(totxv.cn.data)
mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
rownames (mtt.select.b) <- mtt.select.b$MRN
head(mtt.select.b)
mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


totxv.cn.data <- as.matrix (totxv.cn.data)
dim(totxv.cn.data)

totxv.drug.cat <- mtt.select.b$LC50.GROUP
head(totxv.drug.cat)
head(totxv.cn.data)


 totxvi.cn.data <- cn.data.som[cn.result.stat[,1],]
 #totxvi.cn.smooth.data <- cn.data.som[cn.result.stat[,1],]
head(totxvi.cn.data)



mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")
if(drug == "PREDANDDEX"){
  mtt.select.PRED <- mtt.prep("PRED")
  mtt.select.PRED$LC50.GROUP <- ifelse (mtt.select.PRED$LC50 < 0.1, 1, ifelse(mtt.select.PRED$LC50 >= 0.1 & mtt.select.PRED$LC50 < 65, 2, ifelse(mtt.select.PRED$LC50 >= 64, 3, NA)))
  mtt.select.PRED <- subset(mtt.select.PRED, mtt.select.PRED$LIN == "B")
  
  mtt.select.DEX <- mtt.prep("DEX")
  mtt.select.DEX <- subset(mtt.select.DEX, mtt.select.DEX$LIN == "B")
  
  mtt.select.DEX$LC50 <- (mtt.select.DEX$LC50) *8
  mtt.select.DEX.select <- mtt.select.DEX[!(mtt.select.DEX$MRN %in% mtt.select.PRED$MRN),]
  mtt.select <- rbind(mtt.select.PRED, mtt.select.DEX.select)
} 
rownames(mtt.select) <- mtt.select$MRN
 mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
 mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == c("TOTXVI"))
 mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))


 cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxvi.cn.data))
cn.pt.int
 totxvi.cn.data <- totxvi.cn.data[,as.character(cn.pt.int)]
 #totxvi.cn.smooth.data <- totxvi.cn.smooth.data[,as.character(cn.pt.int)]

 mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
 rownames (mtt.select.b) <- mtt.select.b$MRN

 mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


 totxvi.cn.data <- as.matrix (totxvi.cn.data)
 #totxvi.cn.smooth.data <- as.matrix(totxvi.cn.smooth.data)

 dim (totxvi.cn.data)
head(totxvi.cn.data)

totxvi.drug.cat <- mtt.select.b$LC50.GROUP
#################################################
# #############################

both.cn.data <- cn.data.som[cn.result.stat[,1],]
#totxvi.cn.smooth.data <- cn.data.som[cn.result.stat[,1],]
head(both.cn.data)
dim(both.cn.data)


mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")
if(drug == "PREDANDDEX"){
  mtt.select.PRED <- mtt.prep("PRED")
  mtt.select.PRED$LC50.GROUP <- ifelse (mtt.select.PRED$LC50 < 0.1, 1, ifelse(mtt.select.PRED$LC50 >= 0.1 & mtt.select.PRED$LC50 < 65, 2, ifelse(mtt.select.PRED$LC50 >= 64, 3, NA)))
  mtt.select.PRED <- subset(mtt.select.PRED, mtt.select.PRED$LIN == "B")
  
  mtt.select.DEX <- mtt.prep("DEX")
  mtt.select.DEX <- subset(mtt.select.DEX, mtt.select.DEX$LIN == "B")
  
  mtt.select.DEX$LC50 <- (mtt.select.DEX$LC50) *8
  mtt.select.DEX.select <- mtt.select.DEX[!(mtt.select.DEX$MRN %in% mtt.select.PRED$MRN),]
  mtt.select <- rbind(mtt.select.PRED, mtt.select.DEX.select)
} 
rownames(mtt.select) <- mtt.select$MRN
mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")
#mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == c("TOTXV","TOTXVI")
                        mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))
                        
                        
                        cn.pt.int <- intersect (mtt.select.b$MRN, colnames (both.cn.data))
                        cn.pt.int
                        both.cn.data <- both.cn.data[,as.character(cn.pt.int)]
                        #both.cn.smooth.data <- both.cn.smooth.data[,as.character(cn.pt.int)]
                        
                        mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
                        rownames (mtt.select.b) <- mtt.select.b$MRN
                        
                        mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]
                        
                        
                        both.cn.data <- as.matrix (both.cn.data)
                        #both.cn.smooth.data <- as.matrix(both.cn.smooth.data)
                        
                        dim (both.cn.data)
                        head(both.cn.data)
                        
                        both.drug.cat <- mtt.select.b$LC50.GROUP
###################################

xv.cn.clust.p <- rep (NA, times=nrow(cn.result.stat))
xvi.cn.clust.p <- rep (NA, times=nrow(cn.result.stat))
both.cn.clust.p <- rep (NA, times=nrow(cn.result.stat))
xv.cn.clust.cat.p <- rep (NA, times=nrow(cn.result.stat))
xvi.cn.clust.cat.p <- rep (NA, times=nrow(cn.result.stat))
both.cn.clust.cat.p <- rep (NA, times=nrow(cn.result.stat))
#clut.methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "abspearson", "correlation", "abscorrelation")#, "kendall", "spearman")
clut.methods <- "euclidean"

#link.methods <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")#"centroid2"
link.methods <- "ward"
                                        #


#save.image("cntrouble.RData")
for (k in 1:length(link.methods)){
for (j in 1:length(clut.methods)){

for (i in 10:nrow(cn.result.stat)){
  cn.probes.int <- cn.result.stat[1:i,"Probe.Set.ID"]

 # totxv.cn <- totxv.cn.smooth.data[cn.probes.int,]
 totxv.cn <- totxv.cn.data[cn.probes.int,]
 
 totxv.cn.hc2 <- hcluster (t(totxv.cn), method =clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

 
  totxv.cn.clust.p <- fisher.test(table(cutree(totxv.cn.hc2, k=2), totxv.drug.cat))$p.value
  totxv.cn.clust.p

  xv.cn.clust.p[i] <- totxv.cn.clust.p

###
#  totxv.cn <- totxv.cn.data[cn.probes.int,]

#  totxv.cn.hc2 <- hcluster (t(totxv.cn), method = clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

#  totxv.cn.clust.p <- fisher.test(table(cutree(totxv.cn.hc2, k=2), totxv.drug.cat))$p.value
#  totxv.cn.clust.p

#  xv.cn.clust.cat.p[i] <- totxv.cn.clust.p

###
 

  
 totxvi.cn <- totxvi.cn.data[cn.probes.int,]

  totxvi.cn.hc2 <- hcluster (t(totxvi.cn), method =clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

  totxvi.cn.clust.p <- fisher.test(table(cutree(totxvi.cn.hc2, k=2), totxvi.drug.cat))$p.value
  totxvi.cn.clust.p

  xvi.cn.clust.p[i] <- totxvi.cn.clust.p
###BOTH###
###

#  both.cn <- both.cn.data[cn.probes.int,]

#  both.cn.hc2 <- hcluster (t(both.cn), method = clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

#  both.cn.clust.p <- fisher.test(table(cutree(both.cn.hc2, k=2), both.drug.cat))$p.value
#  both.cn.clust.p

#  xvi.cn.clust.p[i] <- both.cn.clust.p

both.cn <- both.cn.data[cn.probes.int,]
dim(both.cn)
both.cn.hc2 <- hcluster (t(both.cn), method =clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

all.cn.clust.p <- fisher.test(table(cutree(both.cn.hc2, k=2), both.drug.cat))$p.value
all.cn.clust.p

both.cn.clust.p[i] <- all.cn.clust.p
###

#  both.cn <- both.cn.data[cn.probes.int,]

#  both.cn.hc2 <- hcluster (t(both.cn), method = clut.methods[j], diag = FALSE, upper = FALSE, link = link.methods[k], members = NULL, doubleprecision = FALSE)

#  both.cn.clust.p <- fisher.test(table(cutree(both.cn.hc2, k=2), both.drug.cat))$p.value
#  both.cn.clust.p

#  xvi.cn.clust.p[i] <- both.cn.clust.p

  ###
#cat (i)

}
cat(clut.methods[j])
cat(link.methods[k])
cat(min(xv.cn.clust.p, na.rm=TRUE))
cat(min(xvi.cn.clust.p, na.rm=TRUE))
cat(min(both.cn.clust.p, na.rm=TRUE))
cat ("\n")
}

}

##############

#q()


 results <- cbind (xv.cn.clust.p, xvi.cn.clust.p,both.cn.clust.p, cn.result.stat$p.b)
results <- as.data.frame (results, as.is=TRUE, stringsAsFactors=FALSE)

colnames (results) <- c("xv.clust", "xvi.clust","both.clust", "p.b")

results$better <- rep(NA, times=nrow(results))
for(i in 10:nrow(results)){

  results[i,"better"] <- dim (results[results[10:nrow(results),"xv.clust"] <= results[i,"xv.clust"]& results[10:nrow(results),"xv.clust"]<= 0.05 & results[10:nrow(results),"xvi.clust"] <= results[i,"xvi.clust"]&results[10:nrow(results),"xvi.clust"]<= 0.05,])[1]

}
results$better[results$better == 0 ]<- NA


 results <- results[-1,]

 results$meta <- pnorm((qnorm(results$xv.clust/2) + qnorm(results$xvi.clust/2))/sqrt(2))*2
 results$order <- (1:nrow(results))+1

#save.image("cncutoffdata.RData")
setwd("/home/rautry/workflow_results_final/PRED/")
 write.csv (results, file="WEScutoffresult.csv")

 results <- subset(results, results$order >= 10)

for (i in 2:nrow(results)){
results[i,"same"] <- ifelse (results[i,"meta"] == results[i+1,"meta"], TRUE, FALSE)
results[i,"same.p.b"] <- ifelse (results[i,"p.b"] == results[i+1,"p.b"], TRUE, FALSE)
}

threshold <- results[min(which (results$meta == min (results$meta))),"order"]

changepoint <- results[which(results$same == FALSE & results$same.p.b == FALSE),"order"]

changepoint <- subset (changepoint, changepoint > threshold)

changepoint <- min (changepoint)
 write.csv (results, file="WEScutoffresult_test.csv")

write.table ((results[which (results$better == min(results$better , na.rm=TRUE)),"p.b"]), "WEScutoff.txt", col.names=FALSE, row.names=FALSE)

write.table (rownames(results)[which (results$better == min(results$better, na.rm=TRUE))], "WEScutoffindex.txt", col.names=FALSE, row.names=FALSE)

write.table (results[which (results$both.clust == min(results$both.clust, na.rm=TRUE)),"p.b"], "WEScutoff_both.txt", col.names=FALSE, row.names=FALSE)

write.table (rownames(results)[which (results$both.clust == min(results$both.clust, na.rm=TRUE))], "WEScutoffindex_both.txt", col.names=FALSE, row.names=FALSE)



