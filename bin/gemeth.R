invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig

protocol <- tolower(as.character(gsub("--", "", args.orig[10])))
idx.orig <- gsub("--", "", args.orig[9])
idx.orig <- as.numeric(idx.orig)
split.num <- 200

protocol
idx.orig

args.orig

#q(save="no")

if(file.exists(paste(protocol, "_ge_meth/output.", sprintf("%03.0f", idx.orig), sep=""))){
  q(save="no")
}

r.lib <- "drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (modeltools, lib.loc=r.lib)
require (coin, lib.loc=r.lib)

ge.results <- read.table("ge_lc50.tsv", header=TRUE, stringsAsFactors=FALSE)
meta.alpha <- read.delim ("gecutoff.txt", header=FALSE)
meta.alpha <- unlist (meta.alpha)
names (meta.alpha) <- NULL
ge.results <- subset (ge.results, ge.results$p.b < meta.alpha)

ge.probe.int <- ge.results$Probe.Set.ID

ge.probe.int

meth.results <- read.table("meth_sig_probes.tsv", header=FALSE, stringsAsFactors=FALSE)

ntests <- nrow(meth.results)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1  
end.row <-  min (block.size*idx.orig, ntests)
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row

meth.probe.int <- meth.results[begin.row:end.row,]

meth.probe.int

if (protocol == "totxv"){
  ge.data.som <- stjude.dxbm.hm.mas5.probe.log2[ge.probe.int,]
  meth.data.som <- meth.data.totxv.som
  colnames(meth.data.som)[colnames(meth.data.som)== "meth.data.totxv.som[, \"25699\"]"] <- "25699"

  pt.int <- intersect (colnames (ge.data.som), colnames(meth.data.som))
  pt.int <- intersect (pt.int, tot.xv[tot.xv$LIN == "B","MRN"])
  
  ge.data.som <- ge.data.som[,pt.int]
  
  meth.data.som <- meth.data.som[unlist(meth.probe.int), as.character(pt.int)]
  class(pt.int)
}


if (protocol == "totxvi"){
  ge.data.som <- stjude.dxbm.xvi.mas5.probe.log2[ge.probe.int,]
  meth.data.som <- meth.data.totxvi.som
  
  pt.int <- intersect (colnames(ge.data.som), colnames(meth.data.som))
  pt.int <- intersect (pt.int, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])
  
  ge.data.som <- ge.data.som[,as.character(pt.int)]

  meth.data.som <- meth.data.som[unlist(meth.probe.int),as.character(pt.int)]
}

if (protocol == "all"){
 ge.data.som.xv <- stjude.dxbm.hm.mas5.probe.log2[ge.probe.int,]
  meth.data.som.xv <- meth.data.totxv.som
  colnames(meth.data.som.xv)[colnames(meth.data.som.xv)== "meth.data.totxv.som[, \"25699\"]"] <- "25699"

  pt.int.xv <- intersect (colnames(ge.data.som.xv), colnames(meth.data.som.xv))
  pt.int.xv <- intersect (pt.int.xv, tot.xv[tot.xv$LIN == "B","MRN"])
  
  ge.data.som.xv <- ge.data.som.xv[,pt.int.xv]
  meth.data.som.xv <- meth.data.som.xv[unlist(meth.probe.int),as.character(pt.int.xv)]

 
 ge.data.som.xvi <- stjude.dxbm.xvi.mas5.probe.log2[ge.probe.int,]
  meth.data.som.xvi <- meth.data.totxvi.som
  
  pt.int.xvi <- intersect (colnames (ge.data.som.xvi), colnames(meth.data.som.xvi))
  pt.int.xvi <- intersect (pt.int.xvi, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])
  
  ge.data.som.xvi <- ge.data.som.xvi[,as.character(pt.int.xvi)]

 meth.data.som.xvi <- meth.data.som.xvi[unlist(meth.probe.int),as.character(pt.int.xvi)]
 meth.data.som <- cbind(meth.data.som.xv, meth.data.som.xvi)
 ge.data.som <- cbind(ge.data.som.xv,ge.data.som.xvi)
 prot <- c(rep ("TOTXV", times=ncol(ge.data.som.xv)), rep ("TOTXVI", times=ncol(ge.data.som.xvi)))
}


#########################
for (i in 1){
  meth.low.pts <- colnames (meth.data.som)[which(meth.data.som[i,] < 0.25)]
  meth.high.pts <-  colnames (meth.data.som)[which(meth.data.som[i,] >= 0.25)]

  meth.low.pts
  meth.high.pts

 if (min (c(length(meth.low.pts), length(meth.high.pts)))>1){

    methge.t <- do.call(rbind, apply(ge.data.som, 1, function(x){t.test(x[meth.high.pts],x[meth.low.pts], na.action=na.rm)}[c("p.value", "statistic")]))

    methge.t <- as.data.frame (methge.t, as.is=TRUE, stringsAsFactors=FALSE)
    methge.t$p.value <- unlist (methge.t$p.value)
    methge.t$statistic <- unlist (methge.t$statistic)
    methge.t$Probe.Set.ID <- rownames (methge.t)
    methge.t <- methge.t[,c("Probe.Set.ID", "p.value", "statistic")]
    head (methge.t)


    methge.w.p <- rep (NA, times=nrow(ge.data.som))
    methge.w.stat <- rep (NA, times=nrow(ge.data.som))

    for (j in 1:nrow(ge.data.som)){
      result <- wilcox_test(unlist(ge.data.som[j,c(meth.high.pts,meth.low.pts)])~factor(c(rep(3, times=length(meth.high.pts)), rep(1, times=length(meth.low.pts))), levels=c(3,1)), distribution="exact")
      methge.w.p[j]  <- pvalue (result)
      methge.w.stat[j] <- statistic (result)
    }

    methge.w <- cbind (Probe.Set.ID=rownames (ge.data.som), methge.w.p, methge.w.stat)
    methge.w <- as.data.frame (methge.w, as.is=TRUE, stringsAsFactors=FALSE)
    methge.w$methge.w.p <- as.numeric(methge.w$methge.w.p)
    methge.w$methge.w.stat <- as.numeric(methge.w$methge.w.stat)

    head (methge.w)
  
if (protocol == "all"){

   methge.lm.p <- rep (NA, times=nrow(ge.data.som))
    methge.lm.stat <- rep (NA, times=nrow(ge.data.som))
           methge.lm <- as.data.frame(cbind(methge.lm.p, methge.lm.stat),header = TRUE, stringsAsFactors = FALSE)

    for(j in 1:nrow(ge.data.som)){
      x.b <- unlist(meth.data.som[i,])
       methge.lm$lm.stat[j] <- summary(lm(unlist(ge.data.som[j,])~x.b+factor(prot)))$coef[2,1]
       methge.lm$lm.p[j] <- summary(lm(unlist(ge.data.som[j,])~x.b+factor(prot)))$coef[2,4]
     }
}
if (protocol != "all"){
   methge.lm.p <- rep (NA, times=nrow(ge.data.som))
    methge.lm.stat <- rep (NA, times=nrow(ge.data.som))
      methge.lm <- as.data.frame(cbind(methge.lm.p, methge.lm.stat),header = TRUE, stringsAsFactors = FALSE)
    for(j in 1:nrow(ge.data.som)){
      x.b <- unlist(meth.data.som[i,])
       methge.lm$lm.stat[j] <- summary(lm(unlist(ge.data.som[j,])~x.b))$coef[2,1]
       methge.lm$lm.p[j] <- summary(lm(unlist(ge.data.som[j,])~x.b))$coef[2,4]
     }
 }
  

 


methge.lm$Probe.Set.ID <- rownames (ge.data.som)
methge.lm <- methge.lm[,3:5]



    methge <- merge (methge.t, methge.w)
    methge <- merge(methge,methge.lm)

    methge <- cbind (methge, rep (rownames(meth.data.som)[i], times=nrow(methge)))

    colnames (methge) <- c("Probe.Set.ID", "t.p", "t.stat", "w.p", "w.stat","lm.stat", "lm.p", "methprobe")

    methge <- methge[,c("Probe.Set.ID", "methprobe", "t.p", "t.stat", "w.p", "w.stat", "lm.stat","lm.p")]

    colnames (methge)[3:8] <- paste(protocol, ".", colnames (methge)[3:8], sep="")

    head (methge)
    methge$methprobe <- as.character(methge$methprobe)

  }
  if (min (c(length(meth.low.pts), length(meth.high.pts)))<2){

  methge <- cbind (rownames(ge.data.som), rep(rownames(meth.data.som)[i], times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)),rep(NA, times=nrow(ge.data.som)),rep(NA, times=nrow(ge.data.som)))   

  colnames (methge) <- c("Probe.Set.ID", "methprobe", "t.p", "t.stat", "w.p", "w.stat","lm.stat","lm.p")
  methge <- as.data.frame (methge, as.is=TRUE, stringsAsFactors=FALSE)
 methge$t.p <- as.numeric(methge$t.p)
  methge$t.stat <- as.numeric (methge$t.stat)
  methge$w.p <- as.numeric (methge$w.p)
 methge$w.stat <- as.numeric (methge$w.stat)
  methge$lm.p <- as.numeric(methge$lm.p)
  methge$lm.stat <- as.numeric(methge$lm.stat)
 colnames (methge)[3:8] <- paste(protocol, ".", colnames (methge)[3:8], sep="")
    
}

}

methgeresult <- methge
#########################
#########################
for (i in 2:nrow(meth.data.som)){
#for (i in 2){  
  meth.low.pts <- colnames (meth.data.som)[which(meth.data.som[i,] < 0.25)]
  meth.high.pts <-  colnames (meth.data.som)[which(meth.data.som[i,] >= 0.25)]

  meth.low.pts
  meth.high.pts

 if (min (c(length(meth.low.pts), length(meth.high.pts)))>1){

    methge.t <- do.call(rbind, apply(ge.data.som, 1, function(x){t.test(x[meth.high.pts],x[meth.low.pts], na.action=na.rm)}[c("p.value", "statistic")]))

    methge.t <- as.data.frame (methge.t, as.is=TRUE, stringsAsFactors=FALSE)
    methge.t$p.value <- unlist (methge.t$p.value)
    methge.t$statistic <- unlist (methge.t$statistic)
    methge.t$Probe.Set.ID <- rownames (methge.t)
    methge.t <- methge.t[,c("Probe.Set.ID", "p.value", "statistic")]
    head (methge.t)


    methge.w.p <- rep (NA, times=nrow(ge.data.som))
    methge.w.stat <- rep (NA, times=nrow(ge.data.som))

    for (j in 1:nrow(ge.data.som)){
      result <- wilcox_test(unlist(ge.data.som[j,c(meth.high.pts,meth.low.pts)])~factor(c(rep(3, times=length(meth.high.pts)), rep(1, times=length(meth.low.pts))), levels=c(3,1)), distribution="exact")
      methge.w.p[j]  <- pvalue (result)
      methge.w.stat[j] <- statistic (result)
    }

    methge.w <- cbind (Probe.Set.ID=rownames (ge.data.som), methge.w.p, methge.w.stat)
    methge.w <- as.data.frame (methge.w, as.is=TRUE, stringsAsFactors=FALSE)
    methge.w$methge.w.p <- as.numeric(methge.w$methge.w.p)
    methge.w$methge.w.stat <- as.numeric(methge.w$methge.w.stat)

    head (methge.w)
 
    if (protocol == "all"){

   methge.lm.p <- rep (NA, times=nrow(ge.data.som))
    methge.lm.stat <- rep (NA, times=nrow(ge.data.som))
       methge.lm <- as.data.frame(cbind(methge.lm.p, methge.lm.stat),header = TRUE, stringsAsFactors = FALSE)

   
    for(j in 1:nrow(ge.data.som)){
      x.b <- unlist(meth.data.som[i,])
       methge.lm$lm.stat[j] <- summary(lm(unlist(ge.data.som[j,])~x.b+factor(prot)))$coef[2,1]
       methge.lm$lm.p[j] <- summary(lm(unlist(ge.data.som[j,])~x.b+factor(prot)))$coef[2,4]
     }
 }
if (protocol != "all"){
  
   methge.lm.p <- rep (NA, times=nrow(ge.data.som))
    methge.lm.stat <- rep (NA, times=nrow(ge.data.som))
      methge.lm <- as.data.frame(cbind(methge.lm.p, methge.lm.stat),header = TRUE, stringsAsFactors = FALSE)

    for(j in 1:nrow(ge.data.som)){
      x.b <- unlist(meth.data.som[i,])
       methge.lm$lm.stat[j] <- summary(lm(unlist(ge.data.som[j,])~x.b))$coef[2,1]
       methge.lm$lm.p[j] <- summary(lm(unlist(ge.data.som[j,])~x.b))$coef[2,4]
     }
 }

 


methge.lm$Probe.Set.ID <- rownames (ge.data.som)
methge.lm <- methge.lm[,3:5]


    methge <- merge (methge.t, methge.w)
    methge <- merge(methge,methge.lm)

    methge <- cbind (methge, rep (rownames(meth.data.som)[i], times=nrow(methge)))

    colnames (methge) <- c("Probe.Set.ID", "t.p", "t.stat", "w.p", "w.stat","lm.stat", "lm.p", "methprobe")

    methge <- methge[,c("Probe.Set.ID", "methprobe", "t.p", "t.stat", "w.p", "w.stat", "lm.stat","lm.p")]

    colnames (methge)[3:8] <- paste(protocol, ".", colnames (methge)[3:8], sep="")

    head (methge)
    methge$methprobe <- as.character(methge$methprobe)
  }

  if (min (c(length(meth.low.pts), length(meth.high.pts)))<2){

  methge <- cbind (rownames(ge.data.som), rep(rownames(meth.data.som)[i], times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)), rep(NA, times=nrow(ge.data.som)),rep(NA, times=nrow(ge.data.som)),rep(NA, times=nrow(ge.data.som)))   

  colnames (methge) <- c("Probe.Set.ID", "methprobe", "t.p", "t.stat", "w.p", "w.stat","lm.stat","lm.p")
  methge <- as.data.frame (methge, as.is=TRUE, stringsAsFactors=FALSE)
 methge$t.p <- as.numeric(methge$t.p)
  methge$t.stat <- as.numeric (methge$t.stat)
  methge$w.p <- as.numeric (methge$w.p)
 methge$w.stat <- as.numeric (methge$w.stat)
  methge$lm.p <- as.numeric(methge$lm.p)
  methge$lm.stat <- as.numeric(methge$lm.stat)
 colnames (methge)[3:8] <- paste(protocol, ".", colnames (methge)[3:8], sep="")
    
}


methgeresult <- rbind(methgeresult,methge)
}



######################



write.table (methgeresult, file=paste(protocol, "_ge_meth/output.", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")
