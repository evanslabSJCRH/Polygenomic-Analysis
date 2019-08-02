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

if(file.exists(paste(protocol, "_ge_cn/output.", cprobe, ".", sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "drworkflow_Rlib"
require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMGEData, lib.loc=r.lib)
require (SJHMSNPData, lib.loc=r.lib)

require (mvtnorm, lib.loc=r.lib)
require (modeltools, lib.loc=r.lib)
require (coin, lib.loc=r.lib)

cn.results <- read.table("cn_lc50.tsv", header=TRUE, stringsAsFactors=FALSE)

cn.cutoff <- read.delim ("cncutoff.txt", header=FALSE)
cn.cutoff <- unlist (cn.cutoff)
names (cn.cutoff) <- NULL

cn.results <- subset (cn.results, cn.results$p.b <= cn.cutoff)


ntests <- nrow(cn.results)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1  
end.row <-  min (block.size*idx.orig, ntests)
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row




cn.probe.int <- cn.results[begin.row:end.row,1]

#snp.probe.int

if (protocol == "totxv"){
ge.data.som <- stjude.dxbm.hm.mas5.probe.log2

load ("/home/rautry/drworkflow_data/cn.data.som.seg.biostat.201607a.RData")
colnames(cn.data.som) <- gsub(".somatic","",colnames(cn.data.som))
pt.int <- intersect (colnames (ge.data.som), colnames(cn.data.som))
pt.int <- intersect (pt.int, tot.xv[tot.xv$LIN == "B","MRN"])

ge.data.som <- unlist(ge.data.som[probe,as.character(pt.int)])

cn.data.som <- cn.data.som[cn.probe.int, as.character(pt.int)]

}


if (protocol == "totxvi"){
ge.data.som <- stjude.dxbm.xvi.mas5.probe.log2


load ("cn.data.som.seg.biostat.201607a.RData")
colnames(cn.data.som) <- gsub(".somatic","",colnames(cn.data.som))
pt.int <- intersect (colnames (ge.data.som), colnames(cn.data.som))
pt.int <- intersect (pt.int, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])

ge.data.som <- unlist(ge.data.som[probe,as.character(pt.int)])

cn.data.som <- cn.data.som[cn.probe.int, as.character(pt.int)]
}

if (protocol == "all"){
  
load ("cn.data.som.seg.biostat.201607a.RData")
colnames(cn.data.som) <- gsub(".somatic","",colnames(cn.data.som))
ge.data.som.xv <-stjude.dxbm.hm.mas5.probe.log2
   pt.int.xv <- intersect (colnames(ge.data.som.xv), colnames(cn.data.som))
pt.int.xv <- intersect (pt.int.xv, tot.xv[tot.xv$LIN == "B","MRN"])

ge.data.som.xv <- unlist(ge.data.som.xv[probe,as.character(pt.int.xv)])

cn.data.som.xv <- cn.data.som[cn.probe.int, as.character(pt.int.xv)]

  ge.data.som.xvi <- stjude.dxbm.xvi.mas5.probe.log2
   pt.int.xvi <- intersect (colnames (ge.data.som.xvi), colnames(cn.data.som))
pt.int.xvi <- intersect (pt.int.xvi, tot.xvi[tot.xvi$lineage_from_immunophenotype == "B","mrn"])

ge.data.som.xvi <- unlist(ge.data.som.xvi[probe,as.character(pt.int.xvi)])

cn.data.som.xvi <- cn.data.som[cn.probe.int, as.character(pt.int.xvi)]
   ge.data.som <- unlist(c(ge.data.som.xv,ge.data.som.xvi))
   cn.data.som <- cbind(cn.data.som.xv,cn.data.som.xvi)
   prot <- c(rep("TOTXV", times = length(pt.int.xv)), rep("TOTXVI", times = length(pt.int.xvi)))
 }

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
    if (protocol != "all" ){
    lm.stat.b[i] <- summary(lm(x.b~ge.data.som))$coef[2,1]
    lm.p.b[i] <- summary(lm(x.b~ge.data.som))$coef[2,4]
  }
     if (protocol == "all" ){
    lm.stat.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,1]
    lm.p.b[i] <- summary(lm(x.b~ge.data.som+factor(prot)))$coef[2,4]
  }
  
  
} 
  cat (i)  

}


cn.result <- cbind (rownames (cn.data.som), rep(probe, times=dim(cn.data.som)[1]), t.p.b, t.stat.b, w.p.b, w.stat.b,lm.stat.b,lm.p.b)

colnames (cn.result) <- c("SNPProbeSetID", "GEProbeSetID", "t.p.b", "t.stat.b", "w.p.b", "w.stat.b","lm.stat.b","lm.p.b")

cn.result <- as.data.frame (cn.result, stringsAsFactors=FALSE, as.is=TRUE)


############################


write.table (cn.result, file=paste(protocol, "_ge_cn/output.", cprobe,".", sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
