
r.lib <- "/home/rautry/drworkflow_Rlib"

require (EvansData)
require (EvansAnalysis)
require (SJHMGEData)
require (SJHMSNPData)

require (mvtnorm)
require (modeltools)
require (coin)
require(dplyr)
require("stringr")
<<geprep>>=
  drug <- "PRED"
mtt.select <- mtt.prep(drug)
source("drugcatadjust.R")





ge.cutoff <- read.delim ("gecutoff.txt", header=FALSE)
ge.cutoff <- unlist (ge.cutoff)
names (ge.cutoff) <- NULL

ge.lc50.result <- read.delim("ge_lc50.tsv", as.is=TRUE, stringsAsFactors=FALSE)
ge.lc50.result <- subset(ge.lc50.result, ge.lc50.result$p.b < ge.cutoff)
ge.lc50.result <- ge.lc50.result [order(ge.lc50.result$p.b),]

mtt.select <- subset (mtt.select, mtt.select$LC50.GROUP %in% c(1,3))
mtt.select <- subset (mtt.select, mtt.select$LIN == "B")
rownames (mtt.select) <- mtt.select$MRN



probes.int.b <- ge.lc50.result[,"Probe.Set.ID"]




all.ge.data.b <- cbind(t(scale(t(stjude.dxbm.hm.mas5.probe.log2[probes.int.b,]))), t(scale(t(stjude.dxbm.xvi.mas5.probe.log2[probes.int.b,]))))


pt.int.b <- intersect(rownames(mtt.select), colnames(all.ge.data.b))

mtt.select <- mtt.select[pt.int.b, ]

all.ge.data.b <- all.ge.data.b[, pt.int.b]

all.ge.hc2 <- hcluster(all.ge.data.b, method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
all.ge.hc2 <- as.dendrogram(all.ge.hc2)

if(drug == "PRED"){
  all.ge.hc2 <- rev(all.ge.hc2)
}

xv.ge.data <- all.ge.data.b[,as.character(mtt.select[mtt.select$PROT == "TOTXV", "MRN"])]
xv.ge.pt.hc <- hcluster(t(xv.ge.data), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xv.ge.clust.p <- fisher.test(table(cutree(xv.ge.pt.hc, k=2), mtt.select[mtt.select$PROT == "TOTXV", "LC50.GROUP"]))$p.value

xv.ge.pt.hc <- as.dendrogram(xv.ge.pt.hc)
xv.ge.pt.hc <- reorder (xv.ge.pt.hc, mtt.select[mtt.select$PROT == "TOTXV", "LC50"], agglo.FUN=sum)

xvi.ge.data <- all.ge.data.b[,as.character(mtt.select[mtt.select$PROT == "TOTXVI", "MRN"])]
xvi.ge.pt.hc <- hcluster(t(xvi.ge.data), method = "correlation", diag = FALSE, upper = FALSE, link = "ward", members = NULL, doubleprecision = FALSE)
xvi.ge.clust.p <- fisher.test(table(cutree(xvi.ge.pt.hc, k=2), mtt.select[mtt.select$PROT == "TOTXVI", "LC50.GROUP"]))$p.value


xvi.ge.pt.hc <- as.dendrogram(xvi.ge.pt.hc)
xvi.ge.pt.hc <- reorder (xvi.ge.pt.hc, mtt.select[mtt.select$PROT == "TOTXVI", "LC50"], agglo.FUN=sum)
xvi.ge.pt.hc <- merge(rev(xvi.ge.pt.hc[[1]]),xvi.ge.pt.hc[[2]])

if(drug == "LASP"){
  xvi.ge.pt.hc <- rev(xvi.ge.pt.hc)}


xv.ge.breaks <- 41
xv.ge.extreme <- max(abs(xv.ge.data), na.rm = TRUE)
xv.ge.breaks <- seq(-xv.ge.extreme, xv.ge.extreme, length=xv.ge.breaks)
xv.min.breaks <- min (xv.ge.breaks)
xv.max.breaks <- max (xv.ge.breaks)
xv.ge.data[xv.ge.data < xv.min.breaks] <- xv.min.breaks
xv.ge.data[xv.ge.data > xv.max.breaks] <- xv.max.breaks
xv.ge.data <- xv.ge.data[order.dendrogram(all.ge.hc2), order.dendrogram(xv.ge.pt.hc)]
xv.ge.col <- mtt.select[mtt.select$PROT == "TOTXV", "LC50.GROUP"]


xvi.ge.breaks <- 41
xvi.ge.extreme <- max(abs(xvi.ge.data), na.rm = TRUE)
xvi.ge.breaks <- seq(-xvi.ge.extreme, xvi.ge.extreme, length=xvi.ge.breaks)
xvi.min.breaks <- min (xvi.ge.breaks)
xvi.max.breaks <- max (xvi.ge.breaks)
xvi.ge.data[xvi.ge.data < xvi.min.breaks] <- xvi.min.breaks
xvi.ge.data[xvi.ge.data > xvi.max.breaks] <- xvi.max.breaks
xvi.ge.data <- xvi.ge.data[order.dendrogram(all.ge.hc2), order.dendrogram(xvi.ge.pt.hc)]
xvi.ge.col <- mtt.select[mtt.select$PROT == "TOTXVI", "LC50.GROUP"]
write.csv(xv.ge.data, file = "xv.ge.data.csv", quote = FALSE)
head(xv.ge.data)
#TODO(spaugh): Add rm here to remove cruft

@ 
drug <- "PREDANDDEX"
#save.image("cntest.RData")
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



#mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))





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
cn.data.som.orig <- cn.data.som

cn.result.stat <- read.csv("/home/rautry/data/mutation_WES/all.lm.b.mutationPREDANDDEXwithprotocol_2018-07-31.csv", as.is=TRUE, header=TRUE, stringsAsFactors=FALSE)
colnames(cn.result.stat) <- c("Probe.Set.ID","stat.b","p.b")
cn.result.stat <- cn.result.stat[order(cn.result.stat$p.b),]
subset(cn.result.stat, cn.result.stat$Probe.Set.ID == "NRAS")
cn.results <- cn.result.stat
cn.cutoff <- read.delim ("WEScutoff.txt", header=FALSE)
cn.cutoff <- unlist (max(cn.cutoff))
names (cn.cutoff) <- NULL

cn.results <- subset (cn.results, cn.results$p.b <= cn.cutoff)
head(cn.results)
dim(cn.results)
cn.probe.int <- cn.results[,1]
ge.lc50 <- read.csv("/home/rautry/workflow_results_final/PRED/ge.lc50.result.b.csv", header = TRUE, stringsAsFactors = FALSE)
head(ge.lc50)
ge.probe.int <- ge.lc50[,3] 
ge.data.som.xv <-stjude.dxbm.hm.mas5.probe.log2
ge.data.som.xv.orig <-stjude.dxbm.hm.mas5.probe.log2
head(ge.data.som.xv.orig)
ge.data.som.xv.orig <- ge.data.som.xv.orig[ge.probe.int,]
dim(ge.data.som.xv.orig)
head(ge.data.som.xv)
pt.int.xv <- intersect (colnames(ge.data.som.xv.orig), colnames(cn.data.som.orig))
pt.int.xv <- intersect (pt.int.xv, tot.xv[tot.xv$LIN == "B","MRN"])
ge.data.som.xvi.orig <- stjude.dxbm.xvi.mas5.probe.log2
ge.data.som.xvi.orig <- ge.data.som.xvi.orig[ge.probe.int,]
head(ge.data.som.xv.orig)
pt.int.xvi <- intersect (colnames(ge.data.som.xvi.orig), colnames(cn.data.som.orig))
pt.int.xvi <- intersect (pt.int.xvi, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])
cn.data.som.xv <- cn.data.som.orig[cn.probe.int, as.character(pt.int.xv)]

cn.data.som.xvi <- cn.data.som.orig[cn.probe.int, as.character(pt.int.xvi)]
cn.data.som<- cbind(cn.data.som.xv,cn.data.som.xvi)
prot <- c(rep("TOTXV", times = length(pt.int.xv)), rep("TOTXVI", times = length(pt.int.xvi)))
out <- NULL
for(j in 1:nrow(ge.data.som.xv.orig)){

probe <- (unlist(rownames(ge.data.som.xv.orig)))[j]  

ge.data.som.xv <- unlist(ge.data.som.xv.orig[probe,as.character(pt.int.xv)])



ge.data.som.xvi <- unlist(ge.data.som.xvi.orig[probe,as.character(pt.int.xvi)])


ge.data.som <- unlist(c(ge.data.som.xv,ge.data.som.xvi))


head(ge.data.som)
#########################



str(ge.data.som)
str (cn.data.som)
dim (cn.data.som)
ge.data.som
#q(save="no")


quantile (ge.data.som)

phenotype <- ifelse (ge.data.som <= quantile(ge.data.som)[2], 1, ifelse(ge.data.som >= quantile(ge.data.som)[4], 3, 2)) 
##########

t.p.b <- rep (NA, times=dim(cn.data.som)[1])
t.stat.b <- rep (NA, times=dim(cn.data.som)[1])
w.p.b <- rep (NA, times=dim(cn.data.som)[1])
w.stat.b <- rep (NA, times=dim(cn.data.som)[1])
lm.stat.b <- rep(NA, times= dim(cn.data.som)[1])
lm.p.b <- rep(NA, times=dim(cn.data.som)[1])
phenotype

test <- unlist(cn.data.som[1,])

test[phenotype %in% c(1,3)]

for (i in 1:nrow(cn.data.som)){
  
  x.b <- unlist(cn.data.som[i,])
 
  if (length (unique(x.b[phenotype %in% c(1,3)])) > 1){
    t.p.b[i] <-try(t.test (x.b[phenotype==3], x.b[phenotype==1])$p.value,silent = TRUE)
    t.result <- try(t.test (x.b[phenotype==3], x.b[phenotype==1])$statistic,silent =TRUE )
    names (t.result) <- NULL
    t.stat.b[i] <- t.result
    result <- wilcox_test(unlist(x.b[phenotype %in% c(1,3)])~factor(phenotype[phenotype %in% c(1,3)], levels=c(3,1)), distribution="exact")
    w.p.b[i]  <- pvalue (result)
    w.stat.b[i] <- statistic (result)
   
      lm.stat.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,1]
      lm.p.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,4]

    
  } 
  cat (i)  


cn.result<- cbind (rownames (cn.data.som), rep(probe, times=dim(cn.data.som)[1]), t.p.b, t.stat.b, w.p.b, w.stat.b,lm.stat.b,lm.p.b)



colnames (cn.result) <- c("SNPProbeSetID", "GEProbeSetID", "t.p.b", "t.stat.b", "w.p.b", "w.stat.b","lm.stat.b","lm.p.b")

cn.result <- as.data.frame (cn.result, stringsAsFactors=FALSE, as.is=TRUE)
}
out <- rbind(out,cn.result)
}
#output CSV file of CN to GE connections


totxv.cn.data <- read.csv("t15.cn.data.heatmap.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
dim(totxv.cn.data)
rownames(totxv.cn.data)<- (totxv.cn.data)[,1]
totxv.cn.data[,1] <- NULL
head(totxv.cn.data)
totxvi.cn.data <- read.csv("t16.cn.data.heatmap.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, )
rownames(totxvi.cn.data)<- (totxvi.cn.data)[,1]
totxvi.cn.data[,1] <- NULL
head(totxvi.cn.data)

if (drug %in% c("LASP","6MP","6TG")){
  cnge.connect.cutoff <- 0.001
}

cnge.connect.cutoff <- 0.001
cnge.connect <- read.csv ("MutationvsGE.csv", as.is=TRUE, stringsAsFactors=FALSE, header=TRUE)

cnge.connect.cutoff
head(cnge.connect)
dim(cnge.connect)
#Matching rownames of data to connection file
cnge.connect2 <- subset (cnge.connect, cnge.connect$SNPProbeSetID %in% rownames(totxv.cn.data))
head(cnge.connect2)


cnge.connect.sig <- cnge.connect2[cnge.connect2$lm.p.b < cnge.connect.cutoff & !is.na(cnge.connect2$lm.p.b) ,]
head(cnge.connect.sig)
#applying columns to dataframe
ge.anno <- as.data.frame(HG.U133_Plus_2.na31.annot.long, header = TRUE, stringsAsFactors = FALSE)
ge.info <- read.delim("/home/rautry/drworkflow_NODUTCH/bin/U133_gene_pos.txt", as.is = TRUE , stringsAsFactors = FALSE,na.strings=c(""," ","NA"))
head(ge.info)
head(ge.anno)
head(ge.info)
ge.info <- ge.info[,c("Probe.Set.ID","Start","End", "Chr","Gene.Symbol")]
ge.info$Start <- as.numeric(ge.info$Start)
ge.info$End <- as.numeric(ge.info$End)
ge.anno <- ge.info

ge.anno.match <- ge.anno[!duplicated(ge.anno$Gene.Symbol),]
cn.ge.indx <- match(cnge.connect.sig$SNPProbeSetID, ge.anno.match$Gene.Symbol)

cn.ge.indx <- cn.ge.indx[!is.na(cn.ge.indx)]

ge.anno.cn <- ge.anno[cn.ge.indx,]

ge.anno.cn <- left_join(cnge.connect.sig,ge.anno.match, , by = c("SNPProbeSetID" ="Gene.Symbol"))
head(ge.anno.cn)
is.na(ge.anno.cn$Chr)
dim(ge.anno.cn)
ge.cn.indx <- match(cnge.connect.sig$GEProbeSetID, ge.anno$Probe.Set.ID)

cn.anno.ge <- ge.anno[ge.cn.indx,]
head(cn.anno.ge)
dim(cn.anno.ge)

head(cnge.connect.sig)
dim(cnge.connect.sig)


cnge.connect.sig$cnorder <- apply (cnge.connect.sig, 1, function (x){which(x["SNPProbeSetID"] == rownames(totxv.cn.data))})
cnge.connect.sig$georder <- apply (cnge.connect.sig, 1, function (x){which(x["GEProbeSetID"] == rownames(xv.ge.data))})
colnames(cnge.connect.sig) <- c("SNPProbeSetID","GEProbeSetID","t.p.b","t.stat.b","w.p.b","w.stat.b","lm.stat.b","lm.p.b","Gene.match","cnorder","georder")
head(cnge.connect.sig)
cnge.connect.sig$loc.start <- unlist(as.numeric(ge.anno.cn$Start))
cnge.connect.sig$loc.end <- unlist(as.numeric(ge.anno.cn$End))
cnge.connect.sig$ge.start <- unlist(as.numeric(cn.anno.ge$Start))
cnge.connect.sig$ge.end <- unlist(as.numeric(cn.anno.ge$End))
cnge.connect.sig$cnchr <- unlist(ge.anno.cn$Chr)
cnge.connect.sig$gechr <- unlist(cn.anno.ge$Chr)
cnge.connect.sig$samechr <- apply (cnge.connect.sig, 1, function (x){ifelse((x["cnchr"] == x["gechr"])== "TRUE", "TRUE","FALSE")})
cnge.connect.sig$distance <-  pmin(abs(cnge.connect.sig$loc.start - cnge.connect.sig$ge.start),abs(cnge.connect.sig$loc.end - cnge.connect.sig$ge.start), abs(cnge.connect.sig$loc.start - cnge.connect.sig$ge.end),abs(cnge.connect.sig$loc.end - cnge.connect.sig$ge.end))

#output CSV file of CN to GE connections

write.csv(cnge.connect.sig, file = "TableS10.csv", row.names = FALSE, quote = FALSE)

@
<<plots, include = FALSE>>
###plots
ge.data.som.select <- unlist(c(ge.data.som.xv.orig["209970_x_at",as.character(pt.int.xv)], ge.data.som.xvi.orig["209970_x_at",as.character(pt.int.xvi)]))
cn.data.som["NRAS",]
plot.meth.p <- summary(lm(cn.data.som["NRAS",]~ge.data.som.select+factor(prot)))$coef[2,4]


boxplot(ge.data.som.select~factor(cn.data.som["NRAS",], levels =c(2,1)), names = c("Wildtype", "Mutant"),ylab = "log2(CASP1 expression)", xlab = "NRAS")
stripchart(ge.data.som.select~factor(cn.data.som["NRAS",], levels =c(2,1)),vertical = TRUE, method = "jitter", pch = 21,cex =1, col= "black",bg ="gray",add = TRUE,overplot = TRUE)
plot.exp <- floor(log10(plot.meth.p))
plot.base <- round(plot.meth.p / 10^plot.exp, 1)
plot.base <- as.character(plot.base)
plot.exp <- as.character(plot.exp)
mtext(bquote(bold(.("p=")*.(plot.base)*" × 10"^.(plot.exp))), cex=1.6, line=-0.1)
###########
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
<<mutvsgefig, fig = TRUE, width=18, height=15.5, include=FALSE>>
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
sciNotation(p=totxv.cn.clust.p, mytext1="SJCRH", mytext2="TOTXV", type = twoline)
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
plot (totxv.cn.hc1, leaflab="none", axes=FALSE, yaxs="i", horiz=TRUE)#need to make combined
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

par(mfrow(c(1,8)), new = T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
mycol <- c(rgb(128,0,128,maxColorValue = 255), rgb(255,191,0, maxColorValue=255))
legend("right", fill = mycol,
       legend = c("Mutation", "No Mutation"))




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



if (drug %in% c("PRED", "VCR", "LASP")){
  x.pos <- c(0, 1, 2, 2.1, 3.1, 3.35, 4.65, 5.95, 6.05, 7.35, 7.45, 8.75, 10.05)/10.05
  x.pos <- c(0, 1, 2, 2.1, 3.1, 3.35, 4.65, 5.65, 5.75, 6.75, 6.85, 7.85, 8.85)/8.85
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

#########################################
#########END Figure Setup################
#########################################
par(fig=c(x.pos[2], x.pos[5], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE, ylim=c(-1,2))
text (0,1.25, "Mutations", cex=2, xpd=TRUE)

par(fig=c(x.pos[2], x.pos[3], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=totxv.cn.clust.p, mytext1="SJCRH", mytext2="TOTXV")

par(fig=c(x.pos[4], x.pos[5], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=totxvi.cn.clust.p, mytext1="SJCRH", mytext2="TOTXVI")

par(fig=c(x.pos[2], x.pos[3], y.pos[2],y.pos[7]), new=T)
image(1:ncol(totxv.cn.data), 1:nrow(totxv.cn.data), t(totxv.cn.data), xlim = 0.5 + c(0, ncol(totxv.cn.data)), ylim = 0.5 + c(0, nrow(totxv.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(128,0,128,maxColorValue = 255),  rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
par(fig=c(x.pos[4], x.pos[5], y.pos[2],y.pos[7]), new=T)
image(1:ncol(totxvi.cn.data), 1:nrow(totxvi.cn.data), t(totxvi.cn.data), xlim = 0.5 + c(0, ncol(totxvi.cn.data)), ylim = 0.5 + c(0, nrow(totxvi.cn.data)), xlab = "", ylab = "", axes = FALSE, col=c(rgb(0,0,1,1), rgb(128,0,128,maxColorValue = 255),  rgb(255,191,0, maxColorValue=255), rgb(1,0,0,0.5), rgb(1,0,0,1)), breaks=c(-1,0.5,1.5,2.5,3.5,4.5), frame.plot=TRUE)
table(totxv.cn.data)
par(fig=c(x.pos[1], x.pos[2], y.pos[2],y.pos[7]), new=T)
plot (totxv.cn.hc1, leaflab="none", axes=FALSE, yaxs="i", horiz=TRUE)#need to make combined ordering

par(fig=c(x.pos[2], x.pos[3], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(totxv.cn.phenotype.b)), col = ifelse(totxv.cn.phenotype.b ==1, "green", ifelse(totxv.cn.phenotype.b == 3, "red", "gray"))[order.dendrogram(totxv.cn.hc2)], axes = FALSE)

par(fig=c(x.pos[4], x.pos[5], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(totxv.cn.phenotype.b)), col = ifelse(totxv.cn.phenotype.b ==1, "green", ifelse(totxv.cn.phenotype.b == 3, "red", "gray"))[order.dendrogram(totxv.cn.hc2)], axes = FALSE, frame.plot=TRUE)

par(fig=c(x.pos[2], x.pos[3], y.pos[8],y.pos[9]), new=T)
plot (totxv.cn.hc2, leaflab="none", axes=FALSE, xaxs="i")

par(fig=c(x.pos[4], x.pos[5], y.pos[8],y.pos[9]), new=T)
plot (totxvi.cn.hc2, leaflab="none", axes=FALSE, xaxs="i")

par(fig=c(x.pos[2], x.pos[5], y.pos[2]-((x.pos[5]-x.pos[2])/10*18/13), y.pos[2]), new=T)
plot (1,1, xlim=c(-10,10), ylim=c(-1,1), xaxs="i", yaxs="i", type="n", axes=FALSE)
curve(-sqrt(1-(x-9)^2)+1, 9, 10, add=TRUE)
curve(-sqrt(1-(x+9)^2)+1, -9, -10, add=TRUE)
curve(sqrt(1-(x-1)^2)-1, 0, 1, add=TRUE)
curve(sqrt(1-(x+1)^2)-1, -1, 0, add=TRUE)
segments(-9,0,-1,0)
segments(1,0,9,0)

cn.meta.p <- pnorm((qnorm((totxv.cn.clust.p)/2) + qnorm(totxvi.cn.clust.p/2))/sqrt(2))*2

par(fig=c(x.pos[2], x.pos[5], y.pos[1],  y.pos[2]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=cn.meta.p, type="oneline")

#########################################
#########END Lower left panel############
#########################################



#########################################
################Middle panel#############
#########################################

#par(fig=c(x.pos[7], x.pos[12], y.pos[17],y.pos[18]), new=T)
#plot(0,0, type="n", xlab="", ylab="", axes=FALSE, ylim=c(-1,2))
#text (0,1, ifelse (drug == "PRED", "Prednisolone", ifelse (drug == "VCR", "Vincristine", ifelse (drug == "LASP", "Asparaginase", ifelse (drug == "CYTU", "Cytarabine", ifelse (drug == "DEX", "Dexamethasone", ifelse (drug == "6MP","6-Mercaptopurine", ifelse (drug == "6TG", "6-Thioguanine", ""))))))), cex=3)

par(fig=c(x.pos[7], x.pos[10], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE, ylim=c(-1,2))
text (0,1.25, "mRNA Expression", cex=2, xpd=TRUE)

#if (drug %in% c("PRED", "VCR", "LASP")){
# par(fig=c(x.pos[7], x.pos[8], y.pos[9],y.pos[10]), new=T)
# plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
# sciNotation(p=dutch.ge.clust.p, mytext2="Dutch")
#}

#if (drug %in% c("PRED", "VCR", "LASP")){
#  par(fig=c(x.pos[7], x.pos[8], y.pos[8],y.pos[9]), new=T)
# plot (dutch.ge.pt.hc, leaflab="none", axes=FALSE, xaxs="i")
#}

par(fig=c(x.pos[7], x.pos[8], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=xv.ge.clust.p, mytext1="SJCRH", mytext2="TOTXV")

par(fig=c(x.pos[7], x.pos[8], y.pos[8],y.pos[9]), new=T)
plot (xv.ge.pt.hc, leaflab="none", axes=FALSE, xaxs="i")

par(fig=c(x.pos[9], x.pos[10], y.pos[9],y.pos[10]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)
sciNotation(p=xvi.ge.clust.p, mytext1="SJCRH", mytext2="TOTXVI")

par(fig=c(x.pos[9], x.pos[10], y.pos[8],y.pos[9]), new=T)
plot (xvi.ge.pt.hc, leaflab="none", axes=FALSE, xaxs="i")

#if (drug %in% c("PRED", "VCR", "LASP")){
# par(fig=c(x.pos[7], x.pos[8], y.pos[7],y.pos[8]), new=T)
# image(cbind(1:length(dutch.ge.col)), col = ifelse(dutch.ge.col ==1, "green", ifelse(dutch.ge.col == 3, "red", "gray"))[order.dendrogram(dutch.ge.pt.hc)], axes = FALSE)
#}

par(fig=c(x.pos[7], x.pos[8], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(xv.ge.col)), col = ifelse(xv.ge.col ==1, "green", ifelse(xv.ge.col == 3, "red", "gray"))[order.dendrogram(xv.ge.pt.hc)], axes = FALSE)

par(fig=c(x.pos[9], x.pos[10], y.pos[7],y.pos[8]), new=T)
image(cbind(1:length(xvi.ge.col)), col = ifelse(xvi.ge.col ==1, "green", ifelse(xvi.ge.col == 3, "red", "gray"))[order.dendrogram(xvi.ge.pt.hc)], axes = FALSE)

#if (drug %in% c("PRED", "VCR", "LASP")){
#par(fig=c(x.pos[7], x.pos[8], y.pos[2],y.pos[7]), new=T)
#image(1:ncol(dutch.ge.data), 1:nrow(dutch.ge.data), t(dutch.ge.data), xlim = 0.5 + c(0, ncol(dutch.ge.data)), ylim = 0.5 + c(0, nrow(dutch.ge.data)), xlab = "", ylab = "", axes = FALSE, col = bluered, breaks = dutch.ge.breaks, frame.plot=TRUE)
#}

par(fig=c(x.pos[7], x.pos[8], y.pos[2],y.pos[7]), new=T)
image(1:ncol(xv.ge.data), 1:nrow(xv.ge.data), t(xv.ge.data), xlim = 0.5 + c(0, ncol(xv.ge.data)), ylim = 0.5 + c(0, nrow(xv.ge.data)), xlab = "", ylab = "", axes = FALSE, col = bluered, breaks = xv.ge.breaks, frame.plot=TRUE)

par(fig=c(x.pos[9], x.pos[10], y.pos[2],y.pos[7]), new=T)
image(1:ncol(xvi.ge.data), 1:nrow(xvi.ge.data), t(xvi.ge.data), xlim = 0.5 + c(0, ncol(xvi.ge.data)), ylim = 0.5 + c(0, nrow(xvi.ge.data)), xlab = "", ylab = "", axes = FALSE, col = bluered, breaks = xvi.ge.breaks, frame.plot=TRUE)

par(fig=c(x.pos[10], x.pos[12], y.pos[2],y.pos[7]), new=T)
plot (all.ge.hc2,  axes=FALSE, yaxs="i", leaflab ="none",horiz=TRUE, xlim=c(0,attributes(all.ge.hc2)$height))

#TODO(spaugh): Watch ratio change here for final output dimensions
par(fig=c(x.pos[7], x.pos[10], y.pos[2]-((x.pos[5]-x.pos[2])/10*18/13), y.pos[2]), new=T)
plot (1,1, xlim=c(-10,10), ylim=c(-1,1), xaxs="i", yaxs="i", type="n", axes=FALSE)
curve(-sqrt(1-(x-9)^2)+1, 9, 10, add=TRUE)
curve(-sqrt(1-(x+9)^2)+1, -9, -10, add=TRUE)
curve(sqrt(1-(x-1)^2)-1, 0, 1, add=TRUE)
curve(sqrt(1-(x+1)^2)-1, -1, 0, add=TRUE)
segments(-9,0,-1,0)
segments(1,0,9,0)

#if (drug %in% c("PRED", "VCR", "LASP")){
# ge.meta.p <- pnorm((qnorm(xv.ge.clust.p/2) + qnorm(xvi.ge.clust.p/2)+qnorm(dutch.ge.clust.p/2))/sqrt(3))*2
#}

#if (!(drug %in% c("PRED", "VCR", "LASP"))){
ge.meta.p <- pnorm((qnorm(xv.ge.clust.p/2) + qnorm(xvi.ge.clust.p/2))/sqrt(2))*2
#}

par(fig=c(x.pos[7], x.pos[10], y.pos[1], y.pos[2]), new=T)
plot (0,0, type="n", xlab="", ylab="", axes=FALSE)

#if (drug %in% c("PRED", "VCR", "LASP")){
# sciNotation(p=ge.meta.p, type="oneline")
#}

#if (!(drug %in% c("PRED", "VCR", "LASP"))){
sciNotation(p=ge.meta.p, type="oneline")
#}

#########################################
#############END Middle panel############
#########################################


#########################################
####START lower left connections ########
#########################################
#par(fig=c(x.pos[5],x.pos[7],y.pos[6],y.pos[16]), new=T)
#par(fig=c(x.pos[12],x.pos[14],y.pos[2],y.pos[12]), new=T)
par(fig=c(x.pos[5],x.pos[7],y.pos[2],y.pos[7]), new=T)


unitsize <- nrow(totxvi.cn.data)/(y.pos[7]-y.pos[2])
unitsize

boxsize <- y.pos[7]-y.pos[2]
boxsize

#lowerlimit <- 0.5+floor(unitsize*boxsize)
#lowerlimit
lowerlimit <- 0.5
lowerlimit

#upperlimit <- 0.5-unitsize*boxsize+floor(unitsize*boxsize)
#upperlimit
upperlimit <- unitsize*boxsize+0.5
upperlimit

geunitsize <- (y.pos[7]-y.pos[2])/nrow(xv.ge.data)

offset <- floor(unitsize*boxsize)-nrow(totxvi.cn.data)
#offset
#offset <- unitsize*((y.pos[12]-y.pos[2])-(y.pos[12]-y.pos[6]))
offset

multi.unit <- (upperlimit-lowerlimit)/(boxsize/geunitsize)
multi.unit


plot(0,0, xlim=c(0,1), ylim=c(lowerlimit, upperlimit), axes=FALSE, xaxs="i", yaxs="i", type="n")

myleftcon <- sample (1:nrow(xv.ge.data), nrow(xv.ge.data)-5, replace=TRUE)
myrightcon <- sample(1:nrow(totxvi.cn.data), nrow(xv.ge.data)-5, replace=TRUE)
numcon <- (myleftcon) + (myrightcon)

#for (i in 1:nrow(snp.ge.result)){
#for (i in 1:length(myleftcon)){

#num.con <- sample((nrow(xv.ge.data)-10):nrow(xv.ge.data), 1)
#   lines (c(0,1), c(snp.ge.result[i,"snporder"],(snp.ge.result[i,"georder"]-0.5)*multi.unit+lowerlimit+offset), col=rgb(0,1,0,0.75))
#lines (c(1,0), c(myrightcon[i],(myleftcon[i]-0.5)*multi.unit+lowerlimit+offset), col=rgb(1,0,1,0.75))

#}

#lines (c(0,1), c(1, (1-0.5)*multi.unit+lowerlimit+offset), col="black")
#for (i in 1:nrow(snp.ge.result)){
# lines (c(0,1), c(cnge.connect.sig[i,"cnorder"],(cnge.connect.sig[i,"georder"]-0.5)*multi.unit+lowerlimit+offset), col=rgb(1,0,1,0.75))
#}

##Drawing cis lines after trans lines will help them show over top of the many trans connections
cnge.connect.sig.cis <- subset(cnge.connect.sig, cnge.connect.sig$cistrans == "CIS")
cnge.connect.sig.trans <- subset(cnge.connect.sig, cnge.connect.sig$cistrans == "TRANS")
#for (i in 1:nrow(snp.ge.result)){
# lines (c(0,1), c(cnge.connect.sig[i,"cnorder"],(cnge.connect.sig[i,"georder"]-0.5)*multi.unit+lowerlimit+offset), col=rgb(1,0,1,0.75))
#}

##Drawing cis lines after trans lines will help them show over top of the many trans connections
cnge.connect.sig.cis <- subset(cnge.connect.sig, cnge.connect.sig$cistrans == "CIS")
cnge.connect.sig.trans <- subset(cnge.connect.sig, cnge.connect.sig$cistrans == "TRANS")
if(nrow(cnge.connect.sig.trans) >0 ){
  for (i in 1:nrow(cnge.connect.sig.trans)){
    if (drug == "PRED"){
      lines (c(1,0), c(cnge.connect.sig.trans[i,"cnorder"],(cnge.connect.sig.trans[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col= rgb(0,1,1,0.5))
    }else{
      lines (c(1,0), c(cnge.connect.sig.trans[i,"cnorder"],(cnge.connect.sig.trans[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col= rgb(0,1,1,0.5))
    }
  }
}
if(nrow(cnge.connect.sig.cis) >0 ){
  for (i in 1:nrow(cnge.connect.sig.cis)){
    if (drug == "PRED"){
      lines (c(1,0), c(cnge.connect.sig.cis[i,"cnorder"],(cnge.connect.sig.cis[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col= rgb(1,0,1,1))
    }else{
      lines (c(1,0), c(cnge.connect.sig.cis[i,"cnorder"],(cnge.connect.sig.cis[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col= rgb(1,0,1,0.5))
    }
  }
}
#for (i in 1:nrow(cnge.connect.sig)){
#if (drug == "PRED"){
# lines (c(1,0), c(cnge.connect.sig[i,"cnorder"],(cnge.connect.sig[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col=ifelse(cnge.connect.sig$cistrans[i] =="CIS" ,rgb(1,0,1,1), rgb(0.73,1,1,0.05)))
#}else{
#lines (c(1,0), c(cnge.connect.sig[i,"cnorder"],(cnge.connect.sig[i,"georder"]-0.5)*multi.unit+lowerlimit+offset),col=ifelse(cnge.connect.sig$cistrans[i] =="CIS" ,rgb(1,0,1,1), rgb(0.73,1,1,0.5)))
#}
#}
#lty= ifelse(cnge.connect.sig$cistrans[i] =="CIS",1,3
#lines(c(1,0),c(0.5,762.4802), col = "red")
#lines(c(1,0),c(300,(56-0.5)*multi.unit+lowerlimit+offset), col = "green", lty = "dotted"))
#plot(0,0, xlim=c(0,1), ylim=c(lowerlimit, upperlimit), axes=FALSE, xaxs="i", yaxs="i", type="n")
#for (i in 1:nrow(cnge.connect.sig)){
 # lines (c(1,0), c(cnge.connect.sig[i,"cnorder"]+offset,(cnge.connect.sig[i,"georder"]-0.5)*multi.unit+lowerlimit), col=rgb(1,0,0,0.75))
#}
#########################################
####END lower left connections ##########
#########################################

@
