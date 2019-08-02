invisible(options(echo = TRUE))

args.orig <- commandArgs()
idx.orig <- as.numeric(gsub("--", "", args.orig[9]))
drug <- as.character(gsub("--", "", args.orig[10])) 
protocol <- tolower(as.character(gsub("--", "", args.orig[11])))

args.orig 

# TODO(spaugh): Incorporate permutation analysis?

split.num <- 200


if(file.exists(paste(protocol, "_snp_lc50/output.",sprintf("%03.0f", idx.orig), sep=""))){
  q(save="no")
}

r.lib <- "/home/rautry/drworkflow_Rlib"

require (EvansData, lib.loc=r.lib)
require (EvansAnalysis, lib.loc=r.lib)
require (SJHMSNPData, lib.loc=r.lib)

if (protocol == "totxv"){
  data ("20111216totxvsnpdatasom")
  snp.data.som <- snp.data.totxv.som
  rm(snp.data.totxv.som)
}

if (protocol == "totxvi"){
  data ("20111216totxvisnpdatasom")
  snp.data.som <- snp.data.totxvi.som
  rm(snp.data.totxvi.som)
}

if (protocol == "all"){
  data ("20111216totxvsnpdatasom")
  data ("20111216totxvisnpdatasom")
  
  xv.only.snps <- setdiff(rownames(snp.data.totxv.som), rownames(snp.data.totxvi.som))
  xvi.add.snps <- as.data.frame(matrix(rep (NA, times=length(xv.only.snps)*ncol(snp.data.totxvi.som)), ncol=ncol(snp.data.totxvi.som)), as.is=TRUE, stringsAsFactors=FALSE)
  colnames (xvi.add.snps) <- colnames (snp.data.totxvi.som)
  rownames (xvi.add.snps) <- xv.only.snps
  snp.data.totxvi.som <- rbind (snp.data.totxvi.som, xvi.add.snps)
  snp.data.totxvi.som <- snp.data.totxvi.som[rownames(snp.data.totxv.som),]
 

  snp.data.som <- cbind (snp.data.totxv.som, snp.data.totxvi.som, stringsAsFactors=FALSE)
  
}

if (protocol == "imputedxv"){
  load (paste("/home/rautry/genomic_results/imputedXV/region_", idx.orig, ".RData", sep=""))
             
  snp.data.som <- genotype.data
  rm(genotype.data)
  rm(info.data)

}


if (protocol == "imputedxvi"){
  my.files <- list.files ("/home/rautry/RepoTest/test/impute_pipeline/results/temp/")
  
  curr.file <- my.files[idx.orig]
  
  curr.file.info <- gsub ("geno", "info", curr.file)
  genotype.data <- read.delim (paste("/home/rautry/RepoTest/test/impute_pipeline/results/temp/", curr.file, sep=""), as.is=TRUE, stringsAsFactor=FALSE, sep=" ", header=FALSE)
  
  genotype.data <- t (genotype.data)
  
  genotype.data.names <- genotype.data[1,]
  genotype.data.names <- strsplit (genotype.data.names, "->")
  genotype.data.names <- as.data.frame(genotype.data.names, as.is=TRUE, stringsAsFactors=FALSE)
  
  genotype.data.names <- genotype.data.names[1,]
  genotype.data.names <- unlist(genotype.data.names)
  
  names(genotype.data.names) <- NULL
  
  genotype.data.names
  
  colnames (genotype.data) <- genotype.data.names
  
  genotype.data <- genotype.data[-c(1:3),]
  
  info.data <- read.delim (paste("/home/rautry/RepoTest/test/impute_pipeline/results/info/", curr.file.info, sep=""), as.is=TRUE, stringsAsFactors=FALSE)
  
  rownames (genotype.data) <- info.data$SNP
  
 # chip.key <- read.csv("/home/rautry/genomic_results/RESOURCES/affysnp/KEYS/20111103_SOM_GERM_SNP6_TOTXVI.csv", as.is=TRUE, stringsAsFactors=FALSE)
  chip.key <- read.csv("/home/rautry/drworkflow/bin/20111103_SOM_GERM_SNP6_TOTXVI.csv", as.is=TRUE, stringsAsFactors=FALSE)
  chip.key$geno.file <- gsub(".CEL", ".birdseed-v2.chp.txt", chip.key$Filename)
  rownames (chip.key) <- chip.key$geno.file
  

  chip.int <- intersect (colnames (genotype.data), rownames (chip.key))
  
  chip.key <- chip.key[chip.int,]
  
  chip.key <- chip.key[colnames(genotype.data),]
  
  colnames (genotype.data) <- chip.key$mrn
  

  snp.data.som <- genotype.data
  rm(genotype.data)
  rm(info.data)

}



ntests <- nrow(snp.data.som)
ntests

block.size <- ntests/split.num
block.size <- floor (block.size)
block.size

begin.row <- ((idx.orig-1)*block.size)+1
#begin.row <- 8000  
end.row <-  min (block.size*idx.orig, ntests)
#end.row <-  9000
if (idx.orig == split.num){ end.row <- ntests}

begin.row
end.row

mtt.select <- mtt.prep (drug)
mtt.select <- subset (mtt.select, mtt.select$LIN == "B")

source("drugcatadjust.R")


rownames (mtt.select) <- mtt.select$MRN

snp.mtt.select.b <- mtt.select

mtt.select <- subset (mtt.select, mtt.select$LC50.GROUP %in% c(1,3))

pt.int.lm <- intersect (rownames(snp.mtt.select.b), colnames(snp.data.som))

pt.int <- intersect (rownames(mtt.select), colnames(snp.data.som))

if (protocol %in% c("totxv", "totxvi", "all")){
snp.data.all <- snp.data.som[begin.row:end.row,as.character(pt.int.lm)]
snp.data.som <- snp.data.som[begin.row:end.row,as.character(pt.int)]
}

if (protocol %in% c("imputedxv", "imputedxvi")){
snp.data.all <- snp.data.som[,as.character(pt.int.lm)]
snp.data.som <- snp.data.som[,as.character(pt.int)]
}



mtt.select <- mtt.select[as.character(pt.int),]
snp.mtt.select.b <- snp.mtt.select.b[as.character(pt.int.lm),]

all.prot.b <- as.factor(snp.mtt.select.b$PROT)

##################

phenotype <- mtt.select$LC50.GROUP

fisher.p.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.01.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.03.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.11.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.13.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.21.b <- rep (NA, times=dim(snp.data.som)[1])
fisher.23.b <- rep (NA, times=dim(snp.data.som)[1])

################
#save.image("trouble.RData")
#q(save="no")
for (i in 1:nrow(snp.data.som)){
  
  x.b <- unlist(snp.data.som[i,])
  if (protocol %in% c("imputedxv", "imputedxvi")){
    x.b <- factor(x.b, levels=c("A/A", "A/C", "C/A", "A/G", "G/A", "A/T", "T/A", "C/C", "C/G", "G/C", "C/T", "T/C", "G/G", "G/T", "T/G", "T/T"))
    x.b <- factor(x.b)
  }
  
  if (protocol %in% c("totxv", "totxvi", "all")){
    fisher.freq.b <- as.data.frame(table (factor (x.b, levels=c("AA", "AB", "BB")), factor(phenotype, levels=c(1,3))), as.is=TRUE, stringsAsFactors=FALSE)$Freq
    

                                        #if (protocol %in% c("imputedxv", "imputedxvi")){
  #  fisher.freq.b <- as.data.frame(table (x.b, factor(phenotype, levels=c(1,3))), as.is=TRUE, stringsAsFactors=FALSE)$Freq
  #}


    fisher.01.b[i] <- fisher.freq.b[1] 
    fisher.11.b[i] <- fisher.freq.b[2]
    fisher.21.b[i] <- fisher.freq.b[3] 
    
    fisher.03.b[i] <- fisher.freq.b[4] 
    fisher.13.b[i] <- fisher.freq.b[5] 
    fisher.23.b[i] <- fisher.freq.b[6]
  
    x.b <- gsub("NN", NA, x.b)
    x.b <- gsub("NC", NA, x.b)
    
    
    if (min(dim(table(x.b))) >= 2){
      snp.result.temp <-  fisher.test (table(factor(x.b, levels=c("AA", "AB", "BB")), factor(phenotype, levels=c(1,3))))
      fisher.p.b[i] <- snp.result.temp$p.value
    }
  }

}
###################
snp.lm.p.b <- rep (NA, times=dim(snp.data.all)[1])
snp.lm.stat.b <- rep (NA, times=dim(snp.data.all)[1])

for (i in 1:nrow(snp.data.all)){
  x.b <- unlist(snp.data.all[i,])
  if (protocol %in% c("totxv", "totxvi", "all")){
    x.b <- gsub("NN", NA, x.b)
    x.b <- gsub("NC", NA, x.b)
    x.b <- gsub("AA", 0, x.b)
    x.b <- gsub("AB", 1, x.b)
    x.b <- gsub("BB", 2, x.b)
    x.b <- as.numeric(x.b)
     af = mean(x.b, na.rm=TRUE)/2
  maf = ifelse(af > 0.5, 1-af, af)
 
    if (protocol != "all"){  
      if(length (table(x.b)) >1 & max(table(x.b)>1)){
         
        curr.result <- summary(lm(log2(snp.mtt.select.b$LC50)~x.b))$coef[2,3:4]
        snp.lm.p.b[i] <- curr.result[2]
        snp.lm.stat.b[i] <- curr.result[1]
      }
    }
    
    if (protocol == "all"){  
      curr.prot.b <- all.prot.b[!is.na(x.b)]
      lc50.data <- snp.mtt.select.b$LC50[!is.na(x.b)]
      x.b <- x.b[!is.na(x.b)]
      if(length (table(x.b)) >1){
        if(length(unique(curr.prot.b))>1 & min(apply(table (curr.prot.b, x.b),1, min))>0 & length(subset(x.b,!is.na(x.b)))>(0.95*length(x.b))&maf>0.01){
         # if(length(unique(curr.prot.b))>1 & max(table(x.b)) >1 & length(subset(x.b,!is.na(x.b)))>(0.95*length(x.b))&maf>0.01){
          curr.result <- summary(lm(log2(lc50.data)~curr.prot.b+x.b))$coef["x.b",3:4]
          snp.lm.p.b[i] <- curr.result[2]
          snp.lm.stat.b[i] <- curr.result[1]
        }
        if(length(unique(curr.prot.b))==1){
          curr.result <- summary(lm(log2(lc50.data)~x.b))$coef["x.b",3:4]
          snp.lm.p.b[i] <- curr.result[2]
          snp.lm.stat.b[i] <- curr.result[1]
        }
      }
    }
  }
  if (protocol %in% c("imputedxv", "imputedxvi")){
    x.b <- factor(x.b, levels=c("A/A", "A/C", "C/A", "A/G", "G/A", "A/T", "T/A", "C/C", "C/G", "G/C", "C/T", "T/C", "G/G", "G/T", "T/G", "T/T"))
    x.b <- factor(x.b)
    curr.prot.b <- all.prot.b[!is.na(x.b)]
    lc50.data <- snp.mtt.select.b$LC50[!is.na(x.b)]
    x.b <- x.b[!is.na(x.b)]
    if(length (table(x.b)) >1){
        curr.result <-
          summary(lm(log2(lc50.data)~x.b))$coef[2,3:4]
        snp.lm.p.b[i] <- curr.result[2]
        snp.lm.stat.b[i] <- curr.result[1]

    }
  }
}

if (!(protocol %in% c("imputedxv", "imputedxvi"))){
  snp.result <- cbind (rownames (snp.data.som), fisher.p.b, fisher.01.b, fisher.03.b, fisher.11.b, fisher.13.b, fisher.21.b, fisher.23.b, snp.lm.p.b, snp.lm.stat.b)
  
  colnames (snp.result) <- c("ProbeSetID", "fisher.p.b", "fisher.01.b", "fisher.03.b", "fisher.11.b", "fisher.13.b", "fisher.21.b", "fisher.23.b", "snp.lm.p.b", "snp.lm.stat.b")
}

if (protocol %in% c("imputedxv", "imputedxvi")){
  snp.result <- cbind (rownames (snp.data.som), snp.lm.p.b, snp.lm.stat.b)
  
  colnames (snp.result) <- c("ProbeSetID", "snp.lm.p.b", "snp.lm.stat.b")
}


snp.result <- as.data.frame (snp.result, stringsAsFactors=FALSE, as.is=TRUE)

##################

write.table (snp.result, file=paste(protocol, "_snp_lc50/output.",sprintf("%03.0f", idx.orig), sep=""), row.names=FALSE, quote=FALSE, sep="\t")


