invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig

protocol <- tolower(as.character(gsub("--", "", args.orig[11])))
mirprobe <- gsub("--", "", args.orig[9])

protocol
mirprobe

if(file.exists(paste(protocol, "_ge_mir/", mirprobe, ".tsv", sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (modeltools, lib.loc=r.lib)
require (coin, lib.loc=r.lib)

ge.results <- read.table("ge_lc50.tsv", header=TRUE, stringsAsFactors=FALSE)
ge.results <- subset (ge.results, ge.results$p.b < 0.05)

ge.probe.int <- ge.results$Probe.Set.ID

ge.probe.int

colnames (ge.results)

if (protocol == "totxv"){
ge.data.som <- stjude.dxbm.hm.mas5.probe.log2[ge.probe.int,]
mir.data.som <- mir.data.totxv.som

pt.int <- intersect (colnames (ge.data.som), colnames(mir.data.som))
pt.int <- intersect (pt.int, tot.xv[tot.xv$LIN == "B","MRN"])

ge.data.som <- ge.data.som[,pt.int]

mir.data.som <- mir.data.som[mirprobe, pt.int]
}


if (protocol == "totxvi"){
  ge.data.som <- stjude.dxbm.xvi.mas5.probe.log2[ge.probe.int,]
  mir.data.som <- mir.data.totxvi.som
  
  pt.int <- intersect (colnames (ge.data.som), colnames(mir.data.som))
  pt.int <- intersect (pt.int, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])
  
  ge.data.som <- ge.data.som[,as.character(pt.int)]

  mir.data.som <- mir.data.som[mirprobe,as.character(pt.int)]
}



if (protocol == "all"){

  ge.data.som.xv <- stjude.dxbm.hm.mas5.probe.log2[ge.probe.int,]
  mir.data.som.xv <- mir.data.totxv.som

  pt.int.xv <- intersect (colnames (ge.data.som.xv), colnames(mir.data.som.xv))
  pt.int.xv <- intersect (pt.int.xv, tot.xv[tot.xv$LIN == "B","MRN"])

  ge.data.som.xv <- ge.data.som.xv[,pt.int.xv]

  mir.data.som.xv <- mir.data.som.xv[mirprobe, pt.int.xv]


  
  ge.data.som.xvi <- stjude.dxbm.xvi.mas5.probe.log2[ge.probe.int,]
  mir.data.som.xvi <- mir.data.totxvi.som
  
  pt.int.xvi <- intersect (colnames (ge.data.som.xvi), colnames(mir.data.som.xvi))
  pt.int.xvi <- intersect (pt.int.xvi, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])
  
  ge.data.som.xvi <- ge.data.som.xvi[,as.character(pt.int.xvi)]
  
  mir.data.som.xvi <- mir.data.som.xvi[mirprobe,as.character(pt.int.xvi)]

  ge.data.som <- cbind(ge.data.som.xv, ge.data.som.xvi)
  mir.data.som <- c(mir.data.som.xv, mir.data.som.xvi)
  prot <- c(rep ("TOTXV", times=ncol(ge.data.som.xv)), rep ("TOTXVI", times=ncol(ge.data.som.xvi)))
}

if (protocol == "all"){
  mirge.lm <- t(apply(ge.data.som, 1, function(x){summary(lm(x~mir.data.som+prot))$coef}[2,c(1,4)]))
}

if (protocol != "all"){
  mirge.lm <- t(apply(ge.data.som, 1, function(x){summary(lm(x~mir.data.som))$coef}[2,c(1,4)]))
}

mirge.lm <- as.data.frame(mirge.lm, as.is=TRUE, stringsAsFactors=FALSE)
mirge.lm$Probe.Set.ID <- rownames (mirge.lm)

mir.data.som

quantile(mir.data.som)

which (mir.data.som <= quantile(mir.data.som)[2])

names (mir.data.som)

mir.low.pts <- names (mir.data.som)[which(mir.data.som <= quantile(mir.data.som)[2])]
mir.high.pts <-  names (mir.data.som)[which(mir.data.som >= quantile(mir.data.som)[4])]


mir.low.pts
mir.high.pts


mirge.t <- do.call(rbind, apply(ge.data.som, 1, function(x){t.test(x[mir.high.pts],x[mir.low.pts], na.action=na.rm)}[c("p.value", "statistic")]))

mirge.t <- as.data.frame (mirge.t, as.is=TRUE, stringsAsFactors=FALSE)
mirge.t$p.value <- unlist (mirge.t$p.value)
mirge.t$statistic <- unlist (mirge.t$statistic)
mirge.t$Probe.Set.ID <- rownames (mirge.t)
mirge.t <- mirge.t[,c("Probe.Set.ID", "p.value", "statistic")]
head (mirge.t)


mirge.w.p <- rep (NA, times=nrow(ge.data.som))
mirge.w.stat <- rep (NA, times=nrow(ge.data.som))

for (i in 1:nrow(ge.data.som)){
  result <- wilcox_test(unlist(ge.data.som[i,c(mir.high.pts,mir.low.pts)])~factor(c(rep(3, times=length(mir.high.pts)), rep(1, times=length(mir.low.pts))), levels=c(3,1)), distribution="exact")
  mirge.w.p[i]  <- pvalue (result)
  mirge.w.stat[i] <- statistic (result)
}

mirge.w <- cbind (Probe.Set.ID=rownames (ge.data.som), mirge.w.p, mirge.w.stat)
mirge.w <- as.data.frame (mirge.w, as.is=TRUE, stringsAsFactors=FALSE)
mirge.w$mirge.w.p <- as.numeric(mirge.w$mirge.w.p)
mirge.w$mirge.w.stat <- as.numeric(mirge.w$mirge.w.stat)

head (mirge.w)

mirge <- merge (mirge.t, mirge.w)

mirge <- merge (mirge, mirge.lm)

mirge <- cbind (mirge, rep (mirprobe, times=nrow(mirge)))

colnames (mirge) <- c("Probe.Set.ID", "t.p", "t.stat", "w.p", "w.stat", "lm.stat", "lm.p", "mirprobe")

mirge <- mirge[,c("Probe.Set.ID", "mirprobe", "t.p", "t.stat", "w.p", "w.stat", "lm.p", "lm.stat")]

colnames (mirge)[3:ncol(mirge)] <- paste(protocol, ".", colnames (mirge)[3:ncol(mirge)], sep="")

head (mirge)

write.table (mirge, file=paste(protocol, "_ge_mir/", mirprobe, ".tsv", sep=""), row.names=FALSE, quote=FALSE, sep="\t")

