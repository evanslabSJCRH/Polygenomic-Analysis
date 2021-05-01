invisible(options(echo = TRUE)) #Turn on echoing

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))
drug <- as.character(gsub("--", "", args.orig[10])) 
protocol <- tolower(as.character(gsub("--", "", args.orig[11])))

split.num <- 200


if(file.exists(paste(protocol, "_cn_lc50/output.",sprintf("%03.0f", idx.orig), sep=""))){
q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"

require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (mvtnorm, lib.loc=r.lib)
require (modeltools, lib.loc=r.lib)
require (coin, lib.loc=r.lib)




##########################
if (protocol == "totxv"){
  load ("/research/rgs01/project_space/relligrp/pgx/common/data/patient/snpchip/cnv/cn.data.som.seg.biostat.20160707.RData")
}
#########################

#########################
if (protocol == "totxvi"){
    load ("/research/rgs01/project_space/relligrp/pgx/common/data/patient/snpchip/cnv/cn.data.som.seg.biostat.20160707.RData") 
}
#########################

#########################
if (protocol == "all"){
  load ("/research/rgs01/project_space/relligrp/pgx/common/data/patient/snpchip/cnv/cn.data.som.seg.biostat.20160707.RData")
}
#########################

#########################


 #cn.data.som <- read.delim("/home/pgx/data/patient/snpchip/cnv/cn.data.som.gene.RData", all =TRUE)


#########################


ntests <- nrow(cn.data.som)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1  
end.row <-  min (block.size*idx.orig, ntests)
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row


mtt.select <- mtt.prep (drug)
rownames (mtt.select) <- mtt.select$MRN
if (drug == "PRED"){
mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 < 0.1, 1, ifelse(mtt.select$LC50 >= 0.1 & mtt.select$LC50 < 65, 2, ifelse(mtt.select$LC50 >= 65, 3, NA)))
}

mtt.select.b <- subset (mtt.select, mtt.select$LIN == "B")

if (protocol == "totxv"){
mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXV")
}
if (protocol == "totxvi"){
mtt.select.b <- subset (mtt.select.b, mtt.select.b$PROT == "TOTXVI")
}
if (protocol == "all"){
mtt.select.b <- mtt.select.b
}

head(mtt.select.b)
#totxvi.mtt <- subset(mtt.select.b, mtt.select.b$PROT == "TOTXV")
#dim(totxv.mtt)


cn.mtt.select <- mtt.select
cn.mtt.select.b <- mtt.select.b
dim(cn.mtt.select)
dim(cn.mtt.select.b)
mtt.select <- subset (mtt.select, mtt.select$LC50.GROUP %in% c(1,3))
mtt.select.b <- subset (mtt.select.b, mtt.select.b$LC50.GROUP %in% c(1,3))

pt.int.lm.b <- intersect (rownames(cn.mtt.select.b), colnames(cn.data.som))
pt.int.lm <- intersect (rownames(cn.mtt.select), colnames(cn.data.som))
pt.int.lm.b
pt.int <- intersect (rownames(mtt.select), colnames(cn.data.som))
cn.data.all <- cn.data.som[begin.row:end.row,as.character(pt.int.lm.b)]
cn.data.som <- cn.data.som[begin.row:end.row,as.character(pt.int)]
#cn.data.all <- cn.data.som[9311:9405,as.character(pt.int.lm.b)]
dim(cn.data.all)
#cn.data.som <- cn.data.som[9311:9405,as.character(pt.int)]
dim(cn.data.som)
mtt.select <- mtt.select[as.character(pt.int),]
cn.mtt.select.b <- cn.mtt.select.b[as.character(pt.int.lm.b),]
all.prot.b <- as.factor(cn.mtt.select.b$PROT)
all.prot.b

##################

phenotype <- mtt.select$LC50.GROUP

t.p <- rep(NA, times=nrow(cn.data.som))
t.stat <- rep(NA, times=nrow(cn.data.som))


w.p <- rep(NA, times=nrow(cn.data.som))
w.stat <- rep(NA, times=nrow(cn.data.som))


for (i in 1:nrow(cn.data.som)){

if (length(unique (unlist(cn.data.som[i,]))) > 1){
  t.p[i] <- t.test(cn.data.som[i, phenotype==3], cn.data.som[i,phenotype==1])$p.value
  stat <-  t.test(cn.data.som[i, phenotype==3], cn.data.som[i,phenotype==1])$statistic
  names (stat) <- NULL
  t.stat[i] <- stat
}
}

for (i in 1:nrow(cn.data.som)){
if (length(unique (unlist(cn.data.som[i,]))) > 1){

  result <- wilcox_test(unlist(cn.data.som[i,phenotype %in% c(1,3)])~factor(phenotype[phenotype %in% c(1,3)], levels=c(3,1)), distribution="exact")
  w.p[i]  <- pvalue (result)
  w.stat[i] <- statistic (result)
}
}

#########
#save.image("CNtro.RData")
###################
cn.lm.p.b <- rep (NA, times=dim(cn.data.all)[1])
cn.lm.stat.b <- rep (NA, times=dim(cn.data.all)[1])

for (i in 1:nrow(cn.data.all)){

  x.b <- unlist(cn.data.all[i,])
 
  if (protocol != "all"){  
    if(length (table(x.b)) >1){
      curr.result <- summary(lm(log2(cn.mtt.select.b$LC50)~x.b))$coef[2,3:4]
      cn.lm.p.b[i] <- curr.result[2]
      cn.lm.stat.b[i] <- curr.result[1]
    }
  }


  if (protocol == "all"){  
    curr.prot.b <- all.prot.b[!is.na(x.b)]
    lc50.data <- cn.mtt.select.b$LC50[!is.na(x.b)]
    x.b <- x.b[!is.na(x.b)]
    if(length (table(x.b)) >1){
	if(length(unique(curr.prot.b))>1){
        curr.result <- summary(lm(log2(lc50.data)~curr.prot.b+x.b))$coef["x.b",3:4]
        cn.lm.p.b[i] <- curr.result[2]
        cn.lm.stat.b[i] <- curr.result[1]
      }
      }
      if(length(unique(curr.prot.b))==1){
        curr.result <- summary(lm(log2(lc50.data)~x.b))$coef["x.b",3:4]
        cn.lm.p.b[i] <- curr.result[2]
        cn.lm.stat.b[i] <- curr.result[1]
      }
    }
  }





cn.result <- cbind (rownames (cn.data.som), t.p, t.stat, w.p, w.stat, cn.lm.p.b, cn.lm.stat.b)

colnames (cn.result) <- c("Probe.Set.ID", "t.p.b", "t.stat.b", "w.p.b", "w.stat.b", "cn.lm.p.b", "cn.lm.stat.b")

##################

write.table (cn.result, file=paste(protocol, "_cn_lc50/output.",sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")


#TODO(spaugh): Make these tests for BandT Tonly and Bonly
