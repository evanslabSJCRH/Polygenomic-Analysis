if (drug == "PRED"){
  mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 < 0.1, 1, ifelse(mtt.select$LC50 >= 0.1 & mtt.select$LC50 < 65, 2, ifelse(mtt.select$LC50 >= 65, 3, NA)))
}


if (drug == "6TG"){
  mtt.select <- read.csv("6TG_LC50_AllCohorts_022414.csv", as.is=TRUE, stringsAsFactors=FALSE)
  mtt.select$LC50.GROUP <- ifelse (mtt.select$LC50 <= 16.7, 1, ifelse(mtt.select$LC50 > 16.7 & mtt.select$LC50 < 51.4, 2, ifelse(mtt.select$LC50 >= 51.4, 3, NA)))
  mtt.select <- mtt.select[!grepl("M", mtt.select$MRN),]
  mtt.select$LIN <- gsub("C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("Down C-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-B-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("PRE-T-ALL", "B", mtt.select$LIN)
  mtt.select$LIN <- gsub("T-ALL", "T", mtt.select$LIN)
  mtt.select <- subset(mtt.select, mtt.select$LIN %in% c("B", "T"))

}


if (drug == "6MP"){
  mtt.select <- read.csv("6MP_LC50_AllCohorts_022414.csv", as.is=TRUE, stringsAsFactors=FALSE)
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
