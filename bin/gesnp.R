invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()

args.orig

protocol <- tolower(as.character(gsub("--", "", args.orig[12])))
idx.orig <- gsub("--", "", args.orig[9])
idx.orig <- as.numeric(idx.orig)
split.num <- 4

probe <- gsub("--", "", args.orig[11])

protocol
idx.orig
split.num
probe
cprobe <- gsub("/", "_", probe)

if(file.exists(paste(protocol, "_ge_snp/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"
#"/nfs_exports/apps/pharmsci/evans_lab/Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)
require (SJHMSNPData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (modeltools, lib.loc=r.lib)
require (coin, lib.loc=r.lib)

getwd()
snp.results <- read.table("snp_lc50.tsv", header=TRUE, stringsAsFactors=FALSE)

snp.alpha <- read.delim ("snpcutoff.txt", header=FALSE)
snp.alpha <- unlist (snp.alpha)
names (snp.alpha) <- NULL


snp.results <- subset (snp.results, snp.results[,"p.b"] < snp.alpha)

ntests <- nrow(snp.results)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1  
end.row <-  min (block.size*idx.orig, ntests)
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row




snp.probe.int <- snp.results[begin.row:end.row,1]

#snp.probe.int


if (protocol == "totxv"){
ge.data.som <- stjude.dxbm.hm.mas5.probe.log2
data ("20111216totxvsnpdatasom")
snp.data.som <- snp.data.totxv.som
rm(snp.data.totxv.som)

pt.int <- intersect (colnames (ge.data.som), colnames(snp.data.som))
pt.int <- intersect (pt.int, tot.xv[tot.xv$LIN == "B","MRN"])

ge.data.som <- unlist(ge.data.som[probe,as.character(pt.int)])

snp.data.som <- snp.data.som[snp.probe.int, as.character(pt.int)]

}


if (protocol == "totxvi"){
ge.data.som <- stjude.dxbm.xvi.mas5.probe.log2

data ("20111216totxvisnpdatasom")
snp.data.som <- snp.data.totxvi.som
rm(snp.data.totxvi.som)

pt.int <- intersect (colnames (ge.data.som), colnames(snp.data.som))
pt.int <- intersect (pt.int, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])

ge.data.som <- unlist(ge.data.som[probe,as.character(pt.int)])

snp.data.som <- snp.data.som[snp.probe.int, as.character(pt.int)]
}

if (protocol == "all"){

 ge.data.som.xv <- stjude.dxbm.hm.mas5.probe.log2
data ("20111216totxvsnpdatasom")
snp.data.som.xv <- snp.data.totxv.som
rm(snp.data.totxv.som)

pt.int.xv <- intersect (colnames (ge.data.som.xv), colnames(snp.data.som.xv))
pt.int.xv <- intersect (pt.int.xv, tot.xv[tot.xv$LIN == "B","MRN"])

ge.data.som.xv <- unlist(ge.data.som.xv[probe,as.character(pt.int.xv)])

snp.data.som.xv <- snp.data.som.xv[snp.probe.int, as.character(pt.int.xv)]

ge.data.som.xvi <- stjude.dxbm.xvi.mas5.probe.log2

data ("20111216totxvisnpdatasom")
snp.data.som.xvi <- snp.data.totxvi.som
rm(snp.data.totxvi.som)

pt.int.xvi <- intersect (colnames (ge.data.som.xvi), colnames(snp.data.som.xvi))
pt.int.xvi <- intersect (pt.int.xvi, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])

ge.data.som.xvi <- unlist(ge.data.som.xvi[probe,as.character(pt.int.xvi)])

snp.data.som.xvi <- snp.data.som.xvi[snp.probe.int, as.character(pt.int.xvi)]
  ge.data.som <- unlist(c(ge.data.som.xv,ge.data.som.xvi))
   snp.data.som <- cbind(snp.data.som.xv,snp.data.som.xvi)
 prot <- c(rep("TOTXV", times = length(pt.int.xv)), rep("TOTXVI", times = length(pt.int.xvi)))
}

#########################



str(ge.data.som)
str (snp.data.som)
dim (snp.data.som)
ge.data.som
#q(save="no")


quantile (ge.data.som)

phenotype <- ifelse (ge.data.som <= quantile(ge.data.som)[2], 1, ifelse(ge.data.som >= quantile(ge.data.som)[4], 3, 2)) 

fisher.p.b <- rep (NA, times=dim(snp.data.som)[1])
#fisher.or.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.01.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.03.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.11.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.13.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.21.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.23.b <- rep (NA, times=dim(snp.data.som)[1])
lm.stat.b <- rep (NA, times=dim(snp.data.som)[1])
lm.p.b <- rep (NA, times=dim(snp.data.som)[1])
for (i in 1:nrow(snp.data.som)){

  x.b <- unlist(snp.data.som[i,])

  fisher.freq.b <- as.data.frame(table (factor (x.b, levels=c("AA", "AB", "BB")), factor(phenotype, levels=c(1,3))), as.is=TRUE, stringsAsFactors=FALSE)$Freq

  fisher.01.b[i] <- fisher.freq.b[1]
  fisher.11.b[i] <- fisher.freq.b[2]
  fisher.21.b[i] <- fisher.freq.b[3]

  fisher.03.b[i] <- fisher.freq.b[4]
  fisher.13.b[i] <- fisher.freq.b[5]
  fisher.23.b[i] <- fisher.freq.b[6]

  x.b <- gsub("NN", NA, x.b)
  x.b <- gsub("NC", NA, x.b)


  if(max(fisher.freq.b)>0){
    if (min(dim(table(x.b))) >= 2){
      snp.result.temp <-  fisher.test (table(factor(x.b, levels=c("AA", "AB", "BB")), factor(phenotype, levels=c(1,3))))
      fisher.p.b[i] <- snp.result.temp$p.value

                                        #if (!is.null(snp.result.temp$estimate)){
                                        #fisher.or.b[i] <- snp.result.temp$estimate
                                        #}

    }
  }
}

for (i in 1:nrow(snp.data.som)){

  x.b <- unlist(snp.data.som[i,])
    x.b <- gsub("NN", NA, x.b)
  x.b <- gsub("NC", NA, x.b)
  x.b <- gsub("AA", 0, x.b)
  x.b <- gsub("AB", 1, x.b)
  x.b <- gsub("BB", 2, x.b)
  x.b <- as.numeric(x.b)
   if (min(dim(table(x.b))) >= 2){
   if (protocol != "all" ){
    lm.stat.b[i] <- summary(lm(x.b~ge.data.som))$coef[2,1]
    lm.p.b[i] <- summary(lm(x.b~ge.data.som))$coef[2,4]
  }
    if (protocol == "all" ){
      if(as.numeric(unlist(table(is.na(x.b))[1])) > length(which(prot == "TOTXV"))){
    lm.stat.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,1]
    lm.p.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,4]
  }
     }
 }
}
#snp.result <- cbind (rownames (snp.data.som), fisher.p.b, fisher.or.b, fisher.01.b, fisher.03.b, fisher.11.b, fisher.13.b, fisher.21.b, fisher.23.b)
snp.result <- cbind (rownames (snp.data.som), rep(probe, times=dim(snp.data.som)[1]), fisher.p.b, fisher.01.b, fisher.03.b, fisher.11.b, fisher.13.b, fisher.21.b, fisher.23.b,lm.stat.b,lm.p.b)

colnames (snp.result) <- c("SNPProbeSetID", "GEProbeSetID", "fisher.p.b", "fisher.01.b", "fisher.03.b", "fisher.11.b", "fisher.13.b", "fisher.21.b", "fisher.23.b","lm.stat.b","lm.p.b")
snp.result <- as.data.frame (snp.result, stringsAsFactors=FALSE, as.is=TRUE)

colnames (snp.result)[3:11] <- paste (protocol, colnames (snp.result)[3:11], sep=".")
##################


write.table (snp.result, file=paste(protocol, "_ge_snp/output.", cprobe,".", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
