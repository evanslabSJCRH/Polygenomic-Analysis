require(EvansAnalysis)
require(EvansData)
require(SJHMGEData)
require(dplyr)
require(gplots)
require(gdata)
require("amap")
require(coin)
require("stringr")
require(lattice)
require(grDevices)
require(ape)

date <- Sys.Date()
drug <- "PREDANDDEX"
<<fileprep>>
  cn.data.som.pre <- read.csv("mutation.call_gene-ALLnoncoding2018-07-30-.csv", header = TRUE, check.names = FALSE, row.names =1)
#cn.data.som.pre <- read.csv("/home/rautry/workflow_results_final/PRED/WES.data.som_replacenumeric2018-08-06.csv", header = TRUE, check.names = FALSE, row.names =1)

head(cn.data.som.pre)
cn.data.som <- cn.data.som.pre %>% 
  mutate_all(funs(str_replace(., "YES", "1"))) %>% mutate_all(funs(str_replace(., "NO", "2")))
rownames(cn.data.som) <- rownames(cn.data.som.pre)
cn.data.som <- apply(cn.data.som, 2, function(x){as.numeric(x)})
rownames(cn.data.som) <- rownames(cn.data.som.pre)
cn.data.som <- apply(cn.data.som, 2, function(x){as.numeric(x)})
rownames(cn.data.som) <- rownames(cn.data.som.pre)
dim(cn.data.som)
  @
<<bothcnsnpprep>>=
  cn.result.stat <- read.csv("all.lm.b.mutationPREDANDDEXwithprotocol_2018-07-31.csv", as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)

colnames(cn.result.stat) <- c("Probe.Set.ID","stat.b","p.b")
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),] 
cn.result.stat <- cn.result.stat[1:500,]
head(cn.result.stat)
dim(cn.result.stat)
both.cn.data <- cn.data.som[rownames(cn.data.som) %in%  cn.result.stat$Probe.Set.ID,]
dim(both.cn.data)


#bothi.cn.data <- cn.data.som[cn.result.stat[,1],]
#bothi.cn.smooth.data <- cn.data.som[cn.result.stat[,1],]
head(both.cn.data)

#rownames (both.cn.data) <- both.cn.data$Probe.Set.ID
#both.cn.data$Probe.Set.ID <- NULL

#both.cn.data <- both.cn.data[,grepl("_Somatic", colnames(both.cn.data))] 
#colnames (both.cn.data) <- gsub ("_Somatic", "", colnames (both.cn.data))

#both.cn.smooth.data <- cn.result.smooth[cn.result.stat[,1],]
#both.cn.smooth.data <- both.cn.smooth.data[,grepl("_Somatic", colnames(both.cn.smooth.data))] 
#colnames (both.cn.smooth.data) <- gsub ("_Somatic", "", colnames (both.cn.smooth.data))

#mtt.select <- mtt.prep(drug)
#source("drugcatadjust.R")
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



mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))
write.csv(mtt.select.b, file = paste("mtt.select.b",date,".csv", sep =""))


cn.pt.int <- intersect (mtt.select.b$MRN, colnames (both.cn.data))

both.cn.data <- both.cn.data[,as.character(cn.pt.int)]
#both.cn.smooth.data <- both.cn.smooth.data[,as.character(cn.pt.int)]

mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
rownames (mtt.select.b) <- mtt.select.b$MRN

mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


both.cn.data <- as.matrix (both.cn.data)
#both.cn.smooth.data <- as.matrix(both.cn.smooth.data)


both.drug.cat <- mtt.select.b$LC50.GROUP


cn.cutoff <- read.delim ("WEScutoff.txt", header=FALSE)
cn.cutoff <- unlist (cn.cutoff)
names (cn.cutoff) <- NULL

cn.result.stat <- subset (cn.result.stat, cn.result.stat$p.b <= max(cn.cutoff))
#cn.result.stat <- subset (cn.result.stat, cn.result.stat$p.b <= ifelse(max(cn.cutoff)<0.05,max(cn.cutoff),0.05))
dim(cn.result.stat)
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),]

cn.cutoff.index <- unlist(read.delim("WEScutoffindex.txt", as.is=TRUE, stringsAsFactors=FALSE, header=FALSE))

names (cn.cutoff.index) <- NULL

cn.result.stat <- cn.result.stat[1:max(cn.cutoff.index),]

dim(cn.result.stat)
cn.probes.int <- cn.result.stat[,1]

cn.probes.int <- cn.probes.int[!is.na(cn.probes.int)]
both.cn.data <- both.cn.data[cn.probes.int,]
dim(both.cn.data)
#both.cn.smooth.data <- both.cn.smooth.data[cn.probes.int,]

#R <- rev(row.names(both.cn.data))
#both.cn.data <- both.cn.data[R,]

dim (both.cn.data)

both.cn.hc1 <- hcluster (both.cn.data, method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

both.cn.hc2 <- hcluster (t(both.cn.data), method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

both.cn.clust.p <- fisher.test(table(cutree(both.cn.hc2, k=2), both.drug.cat))$p.value
both.cn.clust.p

both.cn.hc2 <- reorder(as.dendrogram(both.cn.hc2), mtt.select.b[,"LC50"], agglo.FUN=mean)


both.cn.hc1 <- as.dendrogram(both.cn.hc1)

both.cn.data <- both.cn.data[order.dendrogram(both.cn.hc1), order.dendrogram(both.cn.hc2)]
both.cn.phenotype.b <- mtt.select.b$LC50.GROUP
#rm (cn.result)
#rm (cn.result.smooth)
<<totxvmutprep>>=
drug <- "PREDANDDEX"
  #save.image("cntest.RData")


cn.result.stat <- read.csv("all.lm.b.mutationPREDANDDEXwithprotocol_2018-07-31.csv", as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)

colnames(cn.result.stat) <- c("Probe.Set.ID","stat.b","p.b")
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),] 
cn.result.stat <- cn.result.stat[1:500,]
head(cn.result.stat)
dim(cn.result.stat)
cn.data.som.pre <- read.csv("/home/rautry/data/mutation_WES/mutation.call_gene-ALLnoncoding2018-07-30-.csv", header = TRUE, check.names = FALSE, row.names =1)

head(cn.data.som.pre)
cn.data.som <- cn.data.som.pre %>% 
  mutate_all(funs(str_replace(., "YES", "1"))) %>% mutate_all(funs(str_replace(., "NO", "2")))
rownames(cn.data.som) <- rownames(cn.data.som.pre)
cn.data.som <- apply(cn.data.som, 2, function(x){as.numeric(x)})
rownames(cn.data.som) <- rownames(cn.data.som.pre)
cn.data.som <- apply(cn.data.som, 2, function(x){as.numeric(x)})
rownames(cn.data.som) <- rownames(cn.data.som.pre)
dim(cn.data.som)


totxv.cn.data <- cn.data.som[rownames(cn.data.som) %in%  cn.result.stat$Probe.Set.ID,]
dim(totxv.cn.data)



head(totxv.cn.data)


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

mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXV")

mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))



cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxv.cn.data))

totxv.cn.data <- totxv.cn.data[,as.character(cn.pt.int)]
#totxv.cn.smooth.data <- totxv.cn.smooth.data[,as.character(cn.pt.int)]

mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
rownames (mtt.select.b) <- mtt.select.b$MRN

mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


totxv.cn.data <- as.matrix (totxv.cn.data)


totxv.drug.cat <- mtt.select.b$LC50.GROUP


cn.cutoff <- read.delim ("WEScutoff.txt", header=FALSE)
cn.cutoff <- unlist (cn.cutoff)
names (cn.cutoff) <- NULL

cn.result.stat <- subset (cn.result.stat, cn.result.stat$p.b <= max(cn.cutoff))
#cn.result.stat <- subset (cn.result.stat, cn.result.stat$p.b <= ifelse(max(cn.cutoff)<0.05,max(cn.cutoff),0.05))

cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),]

cn.cutoff.index <- unlist(read.delim("WEScutoffindex.txt", as.is=TRUE, stringsAsFactors=FALSE, header=FALSE))

names (cn.cutoff.index) <- NULL

cn.result.stat <- cn.result.stat[1:max(cn.cutoff.index),]


cn.probes.int <- cn.result.stat[,1]

cn.probes.int <- cn.probes.int[!is.na(cn.probes.int)]
totxv.cn.data <- totxv.cn.data[cn.probes.int,]
dim(totxv.cn.data)
#totxv.cn.smooth.data <- totxv.cn.smooth.data[cn.probes.int,]

#R <- rev(row.names(totxv.cn.data))
#totxv.cn.data <- totxv.cn.data[R,]

dim (totxv.cn.data)

totxv.cn.hc1 <- hcluster (totxv.cn.data, method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxv.cn.hc2 <- hcluster (t(totxv.cn.data), method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxv.cn.clust.p <- fisher.test(table(cutree(totxv.cn.hc2, k=2), totxv.drug.cat))$p.value
totxv.cn.clust.p

totxv.cn.hc2 <- reorder(as.dendrogram(totxv.cn.hc2), mtt.select.b[,"LC50"], agglo.FUN=mean)


totxv.cn.hc1 <- as.dendrogram(totxv.cn.hc1)

#totxv.cn.data <- totxv.cn.data[order.dendrogram(totxv.cn.hc1), order.dendrogram(totxv.cn.hc2)]
totxv.cn.data <- totxv.cn.data[order.dendrogram(both.cn.hc1), order.dendrogram(totxv.cn.hc2)]

totxv.cn.phenotype.b <- mtt.select.b$LC50.GROUP
#rm (cn.result)
#rm (cn.result.smooth)

@ 

<<totxviWESMUTprep>>=
  


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
mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXVI")
mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))


cn.pt.int <- intersect (mtt.select.b$MRN, colnames (totxvi.cn.data))

totxvi.cn.data <- totxvi.cn.data[,as.character(cn.pt.int)]
dim(totxvi.cn.data)
mtt.select.b <- mtt.select.b[mtt.select.b$MRN %in% cn.pt.int,]
rownames (mtt.select.b) <- mtt.select.b$MRN
dim(mtt.select.b)
mtt.select.b <- mtt.select.b[as.character(cn.pt.int),]


totxvi.cn.data <- as.matrix (totxvi.cn.data)

dim (totxvi.cn.data)


totxvi.cn.hc1 <- hcluster (totxvi.cn.data, method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxvi.cn.hc2 <- hcluster (t(totxvi.cn.data), method = "euclidean", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)

totxvi.cn.clust.p <- fisher.test(table(cutree(totxvi.cn.hc2, k=2), mtt.select.b$LC50.GROUP))$p.value
totxvi.cn.clust.p

totxvi.cn.hc2 <- reorder(as.dendrogram(totxvi.cn.hc2), mtt.select.b[,"LC50"], agglo.FUN=mean)


totxvi.cn.hc1 <- as.dendrogram(totxvi.cn.hc1)

totxvi.cn.data <- totxvi.cn.data[order.dendrogram(both.cn.hc1), order.dendrogram(totxvi.cn.hc2)]
totxvi.cn.data
totxvi.cn.phenotype.b <- mtt.select.b$LC50.GROUP


all.cn.data <- cbind(totxv.cn.data, totxvi.cn.data)

write.csv(all.cn.data, file = "all.cn.data.heatmap.csv", quote =FALSE, row.names = TRUE)
write.csv(totxv.cn.data, file = "t15.cn.data.heatmap.csv", quote =FALSE, row.names = TRUE)
write.csv(totxvi.cn.data, file = "t16.cn.data.heatmap.csv", quote =FALSE, row.names = TRUE)

totxv.cn.select <- apply(totxv.cn.data, 1, function(x){ifelse(sum(x)== 2*ncol(totxv.cn.data),"YES","NO")})
#totxv.cn.select <- totxv.cn.select[totxv.cn.select== "YES",]
apply(totxvi.cn.data, 1, function(x){ifelse(sum(x)== 2*ncol(totxvi.cn.data),"YES","NO")})
head(all.cn.data)
@ 
<<sciNotationdef>>=
  
  sciNotation <- function(x=0, y=0, p="", mytext1="", mytext2="", digits=1, cex=1.2, type="twoline") {
    exponent <- floor(log10(p))
    base <- round(p / 10^exponent, digits)
    hght <- strheight("Here", cex=cex)
    p <- "p = "
    base <- as.character(base)
    exponent <- as.character(exponent)
    spacing <- 1.5
    
    if (type=="twoline"){
      text(x, y - (hght * spacing * 1), labels=bquote(bold(.(mytext1))), cex=cex)
      text(x, y - (hght * spacing * 2), labels=bquote(bold(.(mytext2))), cex=cex)
      text(x, y - (hght * spacing * 3), labels=bquote(bold(.(p)*.(base)*" × 10"^.(exponent))), cex=cex)
      
     
    }
    if (type=="oneline"){
      text (x, y, labels=bquote(bold(.(p)*.(base)*" × 10"^.(exponent))), cex=2)
    }
  }

@
#########################################
#############Figure Setup################
#########################################


bluered <- c("#0000FF", "#0D0DFF", "#1B1BFF", "#2828FF", "#3636FF", "#4343FF", 
             "#5151FF", "#5E5EFF", "#6B6BFF", "#7979FF", "#8686FF", "#9494FF", 
             "#A1A1FF", "#AEAEFF", "#BCBCFF", "#C9C9FF", "#D7D7FF", "#E4E4FF", 
             "#F2F2FF", "#FFFFFF", "#FFFFFF", "#FFF2F2", "#FFE4E4", "#FFD7D7", 
             "#FFC9C9", "#FFBCBC", "#FFAEAE", "#FFA1A1", "#FF9494", "#FF8686", 
             "#FF7979", "#FF6B6B", "#FF5E5E", "#FF5151", "#FF4343", "#FF3636", 
             "#FF2828", "#FF1B1B", "#FF0D0D", "#FF0000")



if (drug %in% c("PRED", "PREDANDDEX" ,"VCR", "LASP")){
  x.pos <- c(0, 1, 2, 2.1, 3.1, 3.35, 4.65, 5.95, 6.05, 7.35, 7.45, 8.75, 10.05)/10.05
  x.pos <- c(0, 1, 2, 2.1, 3.1, 3.35, 4.65, 5.45, 5.75, 6.75, 6.85, 7.85, 8.85)/8.85
  #10.3, 11.3, 11.4, 12.4, 13.4)/13.4
  #            1  2  3   4    5     6     7     8     9    10    11    12    13     14    15    16   17    18
  #            M  M  M   M    M
  #                                       G     G     G    G     G     G     G    
  
  y.pos <- c(0, 1,   2, 2.05, 3.05, 3.9425, 5.05, 5.1, 6.1, 7.6)/7.6
  
  #, 8.6,  9.2075, 9.2575, 10.5575, 12.0575, 12.65, 12.7, 13.7, 15.2, 15.3)/15.3
  #            1  2    3  4     5      6       7     8   9     10  11     12      13       14      15       16       17  18    19    20   
  
  
}


#if (!(drug %in% c("PRED", "VCR", "LASP"))){
#  x.pos <- c(0, 1, 2, 2.1, 3.1, 3.35, 4.65, 5.95, 4.65, 5.95, 6.05, 7.35,  8.65,  8.9, 9.9,    10, 11, 12)/12  
# y.pos <- c(0, 1,   2, 2.05, 3.05, 3.9425, 5.05, 5.1, 6.1, 7.6, 8.6,  9.2075, 9.2575, 10.5575, 12.0575, 12.65, 12.7, 13.7, 15.2, 15.3)/15.3
#}

par (mar=c(0,0,0,0), las=1)
plot (0,0,type="n", axes=FALSE, xlab="", ylab="")


#########################################
#########END Figure Setup################
#########################################
#########################################
#############Lower left panel############
#########################################

#par(fig=c(x.pos[2], x.pos[5], y.pos[9],y.pos[10]), new=T)
#par(fig=c(x.pos[2], x.pos[8], y.pos[9],y.pos[10]), new=T)
#plot (0,0, type="n", xlab="", ylab="", axes=FALSE, ylim=c(-1,2))
#text (0,1.25, "Mutation", cex=2, xpd=TRUE)
#par(fig=c(x.pos[2], x.pos[3], y.pos[9],y.pos[10]), new=T)
par(fig=c(x.pos[2], x.pos[5], y.pos[8]+0.07,y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=totxv.cn.clust.p, mytext1="SJCRH", mytext2="TOTXV")
#par(fig=c(x.pos[4], x.pos[5], y.pos[9],y.pos[10]), new=T)
par(fig=c(x.pos[6], x.pos[8], y.pos[8]+0.07,y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=totxvi.cn.clust.p, mytext1="SJCRH", mytext2="TOTXVI")
#par(fig=c(x.pos[2], x.pos[3], y.pos[2],y.pos[7]), new=T)
#par(fig=c(x.pos[2], x.pos[5], y.pos[2],y.pos[7]), new=T)
#image(1:ncol(totxv.cn.data), 1:nrow(totxv.cn.data), t(totxv.cn.data), xlim = 0.5 + c(0, ncol(totxv.cn.data)), ylim = 0.5 + c(0, nrow(totxv.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(0,0,1,0.5), rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
par(fig=c(x.pos[2], x.pos[5], y.pos[2],y.pos[7]), new=T)
image(1:ncol(totxv.cn.data), 1:nrow(totxv.cn.data), t(totxv.cn.data), xlim = 0.5 + c(0, ncol(totxv.cn.data)), ylim = 0.5 + c(0, nrow(totxv.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(128,0,128,maxColorValue = 255), rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
#par(fig=c(x.pos[4], x.pos[5], y.pos[2],y.pos[7]), new=T)
#par(fig=c(x.pos[6], x.pos[8], y.pos[2],y.pos[7]), new=T)
#image(1:ncol(totxvi.cn.data), 1:nrow(totxvi.cn.data), t(totxvi.cn.data), xlim = 0.5 + c(0, ncol(totxvi.cn.data)), ylim = 0.5 + c(0, nrow(totxvi.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(0,0,1,0.5),  rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
par(fig=c(x.pos[6], x.pos[8], y.pos[2],y.pos[7]), new=T)
image(1:ncol(totxvi.cn.data), 1:nrow(totxvi.cn.data), t(totxvi.cn.data), xlim = 0.5 + c(0, ncol(totxvi.cn.data)), ylim = 0.5 + c(0, nrow(totxvi.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(128,0,128,maxColorValue = 255),  rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
#par(fig=c(x.pos[1], x.pos[2], y.pos[2],y.pos[7]), new=T)
par(fig=c(x.pos[1], x.pos[2], y.pos[2],y.pos[7]), new=T)
plot (totxv.cn.hc1, leaflab="perpendicular", axes=FALSE, yaxs="i", horiz=TRUE)#need to make combined
#par(fig=c(x.pos[8], x.pos[10], y.pos[2],y.pos[7]), new=T)
#plot (totxvi.cn.hc1, leaflab="perpendicular", axes=FALSE, yaxs="i", horiz=TRUE, cex.lab =2)#need to make combinedordering

#par(fig=c(x.pos[2], x.pos[3], y.pos[7],y.pos[8]), new=T)
par(fig=c(x.pos[2], x.pos[5], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(totxv.cn.phenotype.b)), col = ifelse(totxv.cn.phenotype.b ==1, "green", ifelse(totxv.cn.phenotype.b == 3, "red", "gray"))[order.dendrogram(totxv.cn.hc2)], axes = FALSE)
#par(fig=c(x.pos[4], x.pos[5], y.pos[7],y.pos[8]), new=T)
par(fig=c(x.pos[6], x.pos[8], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(totxvi.cn.phenotype.b)), col = ifelse(totxvi.cn.phenotype.b ==1, "green", ifelse(totxvi.cn.phenotype.b == 3, "red", "gray"))[order.dendrogram(totxvi.cn.hc2)], axes = FALSE, frame.plot=TRUE)
#par(fig=c(x.pos[2], x.pos[3], y.pos[8],y.pos[9]), new=T)
par(fig=c(x.pos[2], x.pos[5], y.pos[8],y.pos[9]-0.06), new=T)
plot (totxv.cn.hc2, leaflab="none", axes=FALSE, xaxs="i")
#par(fig=c(x.pos[4], x.pos[5], y.pos[8],y.pos[9]), new=T)

par(fig=c(x.pos[6], x.pos[8], y.pos[8],y.pos[9]-0.06), new=T)
plot (totxvi.cn.hc2, leaflab="none", axes=FALSE, xaxs="i")

#par(fig=c(x.pos[2], x.pos[5], y.pos[2]-((x.pos[5]-x.pos[2])/10*18/13), y.pos[2]), new=T)
par(fig=c(x.pos[2], x.pos[8], y.pos[2]-((x.pos[5]-x.pos[2])/10*18/13), y.pos[2]), new=T)

plot (1,1, xlim=c(-10,10), ylim=c(-1,1), xaxs="i", yaxs="i", type="n", axes=FALSE)
curve(-sqrt(1-(x-9)^2)+1, 9, 10, add=TRUE)
curve(-sqrt(1-(x+9)^2)+1, -9, -10, add=TRUE)
curve(sqrt(1-(x-1)^2)-1, 0, 1, add=TRUE)
curve(sqrt(1-(x+1)^2)-1, -1, 0, add=TRUE)
segments(-9,0,-1,0)
segments(1,0,9,0)

cn.meta.p <- pnorm((qnorm((totxv.cn.clust.p)/2) + qnorm(totxvi.cn.clust.p/2))/sqrt(2))*2
#par(fig=c(x.pos[2], x.pos[5], y.pos[1],  y.pos[2]), new=T)
par(fig=c(x.pos[2], x.pos[8], y.pos[1],  y.pos[2]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=cn.meta.p, type="oneline")
#par(fig=c(x.pos[8], x.pos[10], y.pos[2],y.pos[7]), new=T)
#mycol <- c(rgb(0,0,1,0.5), rgb(255,191,0, maxColorValue=255))
#legend("right", fill = mycol,
       #legend = c("Mutation", "No Mutation"))

#########################################
#########END Lower left panel############
#########################################
