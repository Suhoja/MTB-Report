#######################################################
##   Unit Testing toolkit for MTB helper functions   ##
## According to the process denoted in localmaster.r ##
#######################################################
library(testthat)
source("../helpers/get_druggable.r")
source("../helpers/get_levels.r")

########## Load testing data for get_druggable.r ##########
db_GDKD = read.delim("../data/GDKD.csv",sep="\t")
db_CIVIC = read.delim("../data/CIViC.csv",sep="\t")
db_TARGETMERIC = read.delim("../data/TARGET_MERIC.csv",sep="\t")
db_OncoKB = read.csv2("../data/OncoKB.csv", header=T, sep=";")
synonyms = read.csv('../data/cancer_types.csv', header = TRUE,sep="\t")

SNV = read.csv("test_data/test_SNV.csv", header=T, row.names=NULL, stringsAsFactors=F)
CNV = read.csv("test_data/test_CNV.csv", header=T, row.names=NULL, stringsAsFactors=F)
TX = as.array(as.matrix(read.csv("test_data/test_TX.csv", header=T, row.names=NULL, stringsAsFactors=F)))

########## Load testing data for get_levels.r ##########
cancer           = "lymphoma"
cancer_GDKD      = unique(as.character(synonyms[grep(cancer,synonyms$tcga_cancer,ignore.case = T),"knowledge"]))
cancer_CIVIC     = sapply(cancer,function(x) as.character(na.omit(synonyms[grep(x,synonyms$tcga_cancer,ignore.case = T),"civic"])[1]))
cancer_CIVIC     = paste(cancer_CIVIC,collapse = ",")
cancer_CIVIC     = unique(strsplit(cancer_CIVIC,",")[[1]])
cancer_ONCO      = unique(as.character(synonyms[grep(cancer,synonyms$oncokb,ignore.case = T),"oncokb"]))

druggableGDKD = read.csv("test_data/dr_GDKD.csv", header=T, row.names=NULL, stringsAsFactors=F)
druggableCIVIC = read.csv("test_data/dr_CIVIC.csv", header=T, row.names=NULL, stringsAsFactors=F)
druggableTARGET = read.csv("test_data/dr_TARGET.csv", header=T, row.names=NULL, stringsAsFactors=F)
druggableONCO = read.csv("test_data/dr_ONCO.csv", header=T, row.names=NULL, stringsAsFactors=F)

levelsGDKD = get_levels_GDKD(druggableGDKD, cancer_GDKD)
levelsCIVIC = get_levels_CIVIC(druggableCIVIC, cancer_CIVIC)
levelsONCO = get_levels_ONCO(druggableONCO, cancer_ONCO)
levels = merge_levels(levelsGDKD, levelsCIVIC)



#######################################################
##         get_druggable.r function unit tests       ##
#######################################################

#Tests if a dataframe created by a match_XXX_YYY function from mock input data has the correct bounds
testthat::test_that(
  "druggable DF GDKD output bounds test",
  {
  bounds_gdkd_snv = match_SNV_GDKD(SNV, db_GDKD)
  bounds_gdkd_cnv = match_CNV_GDKD(CNV, db_GDKD)
  bounds_gdkd_tx = match_TX_GDKD(TX, db_GDKD)
  bounds_gdkd_wt = match_WT_GDKD(SNV, CNV, cancer_GDKD, db_GDKD)
  
  expect_lte(nrow(bounds_gdkd_snv), nrow(SNV))
  expect_lte(nrow(bounds_gdkd_cnv), nrow(CNV))
  expect_lte(nrow(bounds_gdkd_tx), length(TX))
  expect_lte(nrow(bounds_gdkd_wt), (nrow(CNV) + nrow(SNV)))
  
})

testthat::test_that(
  "druggable DF CIViC output bounds test",
  {
    bounds_civic_snv = match_SNV_CIVIC(SNV, db_CIVIC)
    bounds_civic_cnv = match_CNV_CIVIC(CNV, db_CIVIC)
    bounds_civic_wt = match_WT_CIVIC(SNV, CNV, cancer_CIVIC, db_CIVIC)
    
    expect_lte(nrow(bounds_civic_snv), nrow(SNV))
    expect_lte(nrow(bounds_civic_cnv), nrow(CNV))
    expect_lte(nrow(bounds_civic_wt), (nrow(CNV) + nrow(SNV)))
    
  })

testthat::test_that(
  "druggable DF OncoKB output bounds test",
  {
    bounds_onco_snv = match_SNV_OncoKB(SNV, db_OncoKB)
    bounds_onco_cnv = match_CNV_OncoKB(CNV, db_OncoKB)
    bounds_onco_tx = match_TX_OncoKB(TX, db_OncoKB)
    
    expect_lte(nrow(bounds_onco_snv), nrow(SNV))
    expect_lte(nrow(bounds_onco_cnv), nrow(CNV))
    expect_lte(nrow(bounds_onco_tx), length(TX))
    
  })

testthat::test_that(
  "druggable DF TARGET-MERIC output bounds test",
  {
    bounds_tm = match_TARGET_MERIC(SNV, CNV, TX, db_TARGETMERIC)
    
    expect_lte(nrow(bounds_tm), (nrow(SNV) + nrow(CNV) + length(TX)))
    
  })

#Checks if empty input on a function produces warnings
testthat::test_that(
  "druggable DF empty input test",
  {
    expect_error(match_SNV_CIVIC(data.frame(), db_CIVIC))
    expect_error(match_CNV_CIVIC(data.frame(), db_CIVIC))
    expect_error(match_WT_CIVIC(data.frame(), data.frame(), cancer_CIVIC, db_CIVIC))
    
    expect_error(match_SNV_GDKD(data.frame(), db_GDKD))
    expect_error(match_CNV_GDKD(data.frame(), db_GDKD))
    expect_error(match_TX_GDKD(data.frame(), db_GDKD))
    expect_error(match_WT_GDKD(data.frame(), data.frame(), cancer_GDKD, db_GDKD))
    
    expect_error(match_SNV_OncoKB(data.frame(), db_OncoKB))
    expect_error(match_CNV_OncoKB(data.frame(), db_OncoKB))
    expect_error(match_TX_OncoKB(data.frame(), db_OncoKB))
    
    expect_warning(match_TARGET_MERIC(data.frame(), data.frame(), data.frame(), db_TARGETMERIC))
    
  })

#Tests some use cases with insufficient columns
testthat::test_that(
  "druggable DF wrong format test",
  {
    #For this, it needs to be cleared up how to treat wrongly formatted inputs etc.
    expect_error(match_SNV_GDKD(SNV[,c("Hugo_Symbol", "Variant_Classification")], db_GDKD))
    expect_error(match_CNV_GDKD(CNV[,c("Gene")], db_GDKD))
    expect_error(match_WT_GDKD(SNV[,c("Hugo_Symbol", "Variant_Classification")], CNV[,c("Gene")], cancer_GDKD, db_GDKD))
    
    expect_error(match_SNV_CIVIC(SNV[,c("Hugo_Symbol", "Variant_Classification")], db_CIVIC))
    expect_error(match_CNV_CIVIC(CNV[,c("Gene")], db_CIVIC))
    expect_error(match_WT_CIVIC(SNV[,c("Hugo_Symbol", "Variant_Classification")], CNV[,c("Gene")], cancer_CIVIC, db_CIVIC))
    
    expect_error(match_SNV_OncoKB(SNV[,c("Hugo_Symbol", "Variant_Classification")], db_OncoKB))
    expect_error(match_CNV_OncoKB(CNV[,c("Gene")], db_OncoKB))
    
    expect_error(match_TARGET_MERIC(SNV[,c("Hugo_Symbol", "Variant_Classification")], CNV[,c("Gene")], TX, db_TARGETMERIC))
  })

#######################################################
##         get_levels.r function unit tests          ##
#######################################################

table1 = clean_levels(levels, levelsONCO, synonyms, sort_by="levels")
table2 = clean_levels(levels, levelsONCO, synonyms, sort_by="drug_freq")
table3 = clean_levels(levels, levelsONCO, synonyms, sort_by="genes")

testthat::test_that(
  "no elements missing after sorting",
  {
    expect_equal(nrow(table1), nrow(table2))
    expect_equal(nrow(table1), nrow(table3))
    expect_equal(nrow(table2), nrow(table3))
  }
)

#Could add tests for sorting here.



