require(EvansAnalysis)
require(EvansData)
require(SJHMGEData)
require(dplyr)
require(gplots)
require(gdata)
require("amap")
require(coin)

require("stringr")
date <- Sys.Date()

###PREPROCESSING OF WESDATA
###generation of mutation.call##
mutation1 <- 
mutation2 <- 
head(mutation2)
class(mutation2$GeneName)
PCGPID.file.2 <- 
class(PCGPID.file.2$SJID)
class(mutation2$SJID)
head(mutation2)

colnames(mutation2)[1] <- "SJID"

colnames(mutation1)[3] <- "SJID"
mutation1.order <- mutation1[,c("MRN","GeneName", "SJID","Chr","WU_HG19_Pos","Class","AAChange","ReferenceAllele","MutantAllele")]
mutation2.order <- mutation2[,c("MRN","GeneName", "SJID","Chr","WU_HG19_Pos","Class","AAChange","ReferenceAllele","MutantAllele")]
mutation1.order.sub <- subset(mutation1.order, !is.na(mutation1.order$MRN)& mutation1.order$Class != "silent")
mutation2.order.sub <- subset(mutation2.order, !is.na(mutation2.order$MRN)& mutation2.order$Class != "silent")
#mutation1.order.sub <- subset(mutation1.order.sub,mutation1.order.sub$Class != "UTR_3")
#mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$Class != "UTR_3")
#mutation1.order.sub <- subset(mutation1.order.sub, mutation1.order.sub$Class != "UTR_5")
#mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$Class != "UTR_5")
mutation1.order.sub <- subset(mutation1.order.sub, mutation1.order.sub$Class != "intron")
mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$Class != "intron")
mutation1.order.sub <- subset(mutation1.order.sub, mutation1.order.sub$Class != "")
mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$Class != "")
mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$GeneName != "NO_Somatic_INDELs_FOUND")
mutation2.order.sub <- subset(mutation2.order.sub, mutation2.order.sub$GeneName != "NO_Somatic_SNVs_FOUND")

table(mutation2.order.sub$Class)
table(mutation1.order.sub$Class)

dim(mutation1.order.sub)[1] +dim(mutation2.order.sub)[1]
length(unique(mutation1.order.sub$GeneName))+length(unique(mutation2.order.sub$GeneName))
head(mutation1.order.sub)
mutation.merge.notremoved <- rbind.data.frame(mutation1.order,mutation2.order)
table(mutation.merge.notremoved$Class)
#write.csv(mutation.merge.notremoved, file = "mutation.merge.notremoved.csv")
mutation.merge.notremoved.annovar <- mutation.merge.notremoved[,c("Chr","WU_HG19_Pos","ReferenceAllele","MutantAllele")]
colnames(mutation.merge.notremoved.annovar) <- c("CHROM", "POS", "REF", "ALT")
mutation.merge.notremoved.annovar$QUAL <- rep(NA, times = nrow(mutation.merge.notremoved.annovar))
mutation.merge.notremoved.annovar$ID <- rep(NA, times = nrow(mutation.merge.notremoved.annovar))
mutation.merge.notremoved.annovar$INFO <- rep(NA, times = nrow(mutation.merge.notremoved.annovar))
mutation.merge.notremoved.annovar$FILTER <- rep("PASS", times = nrow(mutation.merge.notremoved.annovar))
mutation.merge.notremoved.annovar <- mutation.merge.notremoved.annovar[,c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")]
#write.csv(mutation.merge.notremoved.annovar,file ="mutation.merge.notremoved.annovar_72018.csv", quote = FALSE, row.names = FALSE)
mutation.merge <- rbind.data.frame(mutation1.order.sub,mutation2.order.sub )

mutation.merge$GeneName <- as.character(mutation.merge$GeneName)
#write.csv(mutation.merge, file = "mutation.file.TAP.allsilentremoved.csv", quote = FALSE, row.names = TRUE)
mutation1.order.sub.trim <-  mutation1.order.sub[,c("MRN","GeneName")]
mutation2.order.sub.trim <- mutation2.order.sub[,c("MRN","GeneName")]
mutation.merge.full <- rbind(mutation1.order.sub.trim,mutation2.order.sub.trim)
dim(mutation.merge)
head(mutation.merge.full)
unique(mutation.merge.full$MRN)
mutation.call <- matrix(NA,length(unique(mutation.merge.full$GeneName)),length(unique(mutation.merge.full$MRN)))
mutation.call <- as.data.frame(mutation.call)
colnames(mutation.call)<- unique(mutation.merge.full$MRN)
row.names(mutation.call)<- sort(unique(as.character(mutation.merge.full$GeneName)))
head(row.names(mutation.call))
mutation.call["ABCD2",]

class(rownames(mutation.call))
colnames(mutation.call)[2]




#for(i in 1:length(unique(mutation.merge.full$GeneName))){
for(j in 1:length(unique(mutation.merge.full$MRN))){
 #for(j in 1:2){

mutation.call[,paste(unique(mutation.merge.full$MRN))[j]] <-do.call(rbind,lapply(split(mutation.merge.full,mutation.merge.full$GeneName),function(chunk){ifelse(grepl(unique(mutation.merge.full$MRN)[j],as.character(chunk))[1],"YES","NO")}))

}
#}

####
#j <-1
test.df <- split(mutation.merge.full,mutation.merge.full$GeneName)
test.df[1]
grepl(unique(mutation.merge.full$MRN)[j],as.character(test.df[i])[1])
test <- do.call(rbind,lapply(split(mutation.merge.full,mutation.merge.full$GeneName),function(chunk){ifelse(grepl(unique(mutation.merge.full$MRN)[j],as.character(chunk))[1],"YES","NO")}))
head(test)

head(mutation.call)

date <- Sys.Date()
#write.csv(mutation.call, file = paste("mutation.call-",date,"-.csv", sep = ""), quote = FALSE, row.names = TRUE)
#mutation.call$Gene <- unlist(row.names(mutation.call))
#write.csv(mutation.call, file = paste("mutation.call_gene-ALLnoncoding",date,"-.csv", sep = ""), quote = FALSE, row.names = TRUE)

#mutation.call <- read.csv("/rgs01/project_space/evansgrp/dr_genomics/common/data/mutation_WES/mutation.call_gene-ALLnoncoding2018-07-30-.csv", header = TRUE, check.names = FALSE, row.names =1)
#dim(mutation.call)
cn.data.som.pre <- read.csv("/rgs01/project_space/evansgrp/dr_genomics/common/data/mutation_WES/mutation.call_gene-ALLnoncoding2018-07-30-.csv", header = TRUE, check.names = FALSE, row.names =1)
head(cn.data.som.pre)

cn.data.som <- cn.data.som.pre %>% 
  mutate_all(funs(str_replace(., "YES", "1"))) %>% mutate_all(funs(str_replace(., "NO", "2")))
rownames(cn.data.som) <- rownames(cn.data.som.pre)
write.csv(cn.data.som , file = "mutation.call_gene-ALLnoncoding2018-07-30-numeric.csv", quote = FALSE,row.names = TRUE)
######end file creation
cn.data.som  <- read.csv("mutation.call_gene-ALLnoncoding2018-07-30-numeric.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
head(cn.data.som)
mutation.call <- cn.data.som
head(mutation.call)
dim(mutation.call)
colnames(mutation.call)
rownames(mutation.call)

##mtt.prep
analysis <- "NOprotocol_allNC"
#BONLY
drug <- "PREDANDDEX"
mtt.select <- mtt.prep(drug)
mtt.select <- subset(mtt.select, mtt.select$LIN == "B")
if (drug == "6MP"){
  mtt.select <- read.csv("6MP_LC50_AllCohorts_022414.csv", as.is=TRUE, stringsAsFactors=FALSE)
  mtt.select <- subset(mtt.select, mtt.select$LIN == "B")
  mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 <= 305, 1, ifelse(mtt.select$LC50 > 305 & mtt.select$LC50 < 1170, 2, ifelse(mtt.select$LC50 >= 1170, 3, NA)))
  mtt.select <- mtt.select[!grepl("M", mtt.select$MRN),]
  mtt.select$LIN <- gsub("C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("Down C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-B-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-T-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("T-ALL", "T", mtt.select$LIN)
  mtt.select <- subset(mtt.select, mtt.select$LIN %in% c("B", "T"))
  mtt.select <- mtt.select[!grepl("M", mtt.select$MRN),]
}
if (drug == "6TG"){
  mtt.select <- read.csv("6TG_LC50_AllCohorts_022414.csv", as.is=TRUE, stringsAsFactors=FALSE)
  mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 <= 16.7, 1, ifelse(mtt.select$LC50 > 16.7 & mtt.select$LC50 < 51.4, 2, ifelse(mtt.select$LC50 >= 51.4, 3, NA)))
  mtt.select <- subset(mtt.select, mtt.select$LIN == "B")
  mtt.select <- mtt.select[!grepl("M", mtt.select$MRN),]
  mtt.select$LIN <- gsub("C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("Down C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-B-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-T-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("T-ALL", "T", mtt.select$LIN)
  mtt.select <- subset(mtt.select, mtt.select$LIN %in% c("B", "T"))
  
}
if(drug == "PRED"){
  mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 < 0.1, 1, ifelse(mtt.select$LC50 >= 0.1 & mtt.select$LC50 < 65, 2, ifelse(mtt.select$LC50 >= 64, 3, NA)))
}
if(drug == "PREDANDDEX"){
  mtt.select.PRED <- mtt.prep("PRED")
  mtt.select.PRED <- subset(mtt.select.PRED, mtt.select.PRED$LIN == "B")
  mtt.select.PRED$LC50.GROUP <- ifelse (mtt.select.PRED$LC50 < 0.1, 1, ifelse(mtt.select.PRED$LC50 >= 0.1 & mtt.select.PRED$LC50 < 65, 2, ifelse(mtt.select.PRED$LC50 >= 64, 3, NA)))
 # mtt.select.PRED$LC50 <- scale(mtt.select.PRED$LC50)

  mtt.select.DEX <- mtt.prep("DEX")
  mtt.select.DEX <- subset(mtt.select.DEX, mtt.select.DEX$LIN == "B")
  #mtt.select.DEX$LC50 <- scale(mtt.select.DEX$LC50)
 mtt.select.DEX$LC50 <- (mtt.select.DEX$LC50)*8
  mtt.select.DEX.select <- mtt.select.DEX[!(mtt.select.DEX$MRN %in% mtt.select.PRED$MRN),]
  mtt.select <- rbind(mtt.select.PRED,mtt.select.DEX.select)
  #mtt.select$LC50 <- scale(mtt.select$LC50)

}
dim(mtt.select)
mutation.call.mtt <-  mutation.call[,colnames(mutation.call) %in% mtt.select$MRN]
head(mutation.call.mtt)
mtt.select.mutation <- mtt.select[mtt.select$MRN %in% colnames(mutation.call),]
mtt.select.mutation$MRN
head(mtt.select.mutation)
#mtt.select.mutation <- mtt.select[mtt.select$MRN %in% mutation.call.mtt,]
dim(mtt.select.mutation)
dim(mutation.call.mtt)
mutation.call.mtt <- mutation.call.mtt[,mtt.select.mutation$MRN]
mutation.call.mtt <- as.data.frame(t(mutation.call.mtt))
rownames((mutation.call))
head(mutation.call.mtt)
mutation.call.mtt$MRN <- unlist(rownames(mutation.call.mtt))
mutation.call.mtt.MRN <- unlist(rownames(mutation.call.mtt))
mutation.call.LC50 <- merge(mtt.select,mutation.call.mtt, by = "MRN")
mutation.call.LC50[,1]
mutation.call.mtt$MRN <- NULL
#mutation.call.LC50 <- merge(mtt.select,mutation.call.mtt, by = "MRN")
mutation.call.LC50$MRN

mutation.call.mtt.names <- unlist(rownames(mutation.call.mtt))
mutation.call.mtt.genes <- unlist(colnames(mutation.call.mtt))
colnames(mutation.call.mtt)
#mutation.call.mtt <- mutation.call.mtt[,1:6892]
head(mutation.call.LC50)
rownames(mutation.call.mtt) <- NULL
lm.stat <- rep(NA, times = (ncol(mutation.call.mtt)))
lm.p <- rep(NA, times = (ncol(mutation.call.mtt)))
table(mutation.call.mtt == "YES")
for(i in 1:ncol(mutation.call.mtt)){
if((length(levels(factor(mutation.call.mtt[,i]))) >1)){
lm.stat[i]<-  unlist(summary(lm(unlist(log10(mtt.select.mutation$LC50) ~ unlist(factor(mutation.call.mtt[,i]))+mtt.select.mutation$PROT)))$coef[2,1])
lm.p[i]<- unlist(summary(lm(unlist(log10(mtt.select.mutation$LC50) ~ unlist(factor(mutation.call.mtt[,i]))+mtt.select.mutation$PROT)))$coef[2,4])
#lm.stat[i]<-  unlist(summary(lm(unlist(log10(mtt.select.mutation$LC50) ~ unlist(factor(mutation.call.mtt[,i])))))$coef[2,1])
#lm.p[i]<- unlist(summary(lm(unlist(log10(mtt.select.mutation$LC50) ~ unlist(factor(mutation.call.mtt[,i])))))$coef[2,4])


}
if((length(levels(factor(mutation.call.mtt[,i]))) <2)){
lm.stat[i]<- "NA"
lm.p[i]<- "NA"

}
}
t.test(mutation.call.mtt[mtt.select.mutation$LC50.GROUP==3],mutation.call.mtt[mtt.select.mutation$LC50.GROUP==1], na.action=na.rm)
snp.result.temp <-  fisher.test (table(factor(mutation.call.mtt[,3751], levels=c( 1, 2)), factor(mtt.select.mutation$LC50.GROUP, levels=c(1,3))))
lm.stat[5]
length(lm.stat)
unlist(summary(lm(unlist(log10(mutation.call.LC50$LC50) ~ unlist(factor(mutation.call.mtt[,"KRAS"])))))$coef[2,4])
plot(unlist(unlist(log10(mutation.call.LC50$LC50) ~ unlist(factor(mutation.call.LC50[,"ZSCAN20"])))))

all.lm.b.mutation <- cbind.data.frame(unlist(colnames(mutation.call.mtt)),lm.stat,lm.p)
head(all.lm.b.mutation)
colnames(all.lm.b.mutation)[1] <-  "Gene"
all.lm.b.mutation <- all.lm.b.mutation[order(lm.p),]
colnames(mutation.call.LC50)[1:10]
write.csv(mutation.call.LC50, file = paste("mutation.call.LC50.b",drug,analysis,"_",date,".csv", sep =""), quote = FALSE, row.names = FALSE)
write.csv(all.lm.b.mutation, file = paste("all.lm.b.mutation",drug,analysis,"_",date,".csv", sep =""), quote = FALSE, row.names = FALSE)

