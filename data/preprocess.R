### Preposition from Charlotte Höltermann ###
### aiming to replace preprocess.r. and update_civic_db.r by summarizing both
### outputs a txt file with current versions, too 

library(timeSeries)
library(timeDate)
library(openxlsx)

# Checks if folder exists
folder_exists <- function(folder_path){
  return (file.exists(folder_path))
}


# Checks if a folder contains a databasefile; ignores meta.txt
file_exists <- function(folder_path){
  files <- list.files(path=folder_path)
  files <- files[files != "meta.txt"]
  return (length(files) >= 1)
}

download_civic_nightly <- function(path_civic, filename){
  
  #download the civic-nightly-file
  url <- "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
  download.file(url, paste(path_civic, filename, sep=""))
  
  #create version- and date-file
  meta_civic_path = paste(path_civic, "meta.txt", sep="")
  writer <- file(meta_civic_path)
  date_line <- paste("date", Sys.Date(), sep=":")
  version_line <- paste("version", "1", sep=":")
  writeLines(c(date_line, version_line), writer)
  close(writer)
}


database_path = "/opt/MTB/data/"
path_civic <- "CIVIC/"
path_gdkd <- "GDKD/"
path_oncokb <- "ONCOKB/"

read_flag_gdkd <- FALSE
read_flag_oncokb = FALSE

output_path <- "/opt/MTB/data/"


if(folder_exists(database_path)){
  print("Database-Path exists.")
  setwd(database_path)
} else {
  print("No Database-Path found. Download CIViC-Nightly File into the automatically created Folder.")
  setwd(database_path)
  download_civic_nightly(path_civic, "nightly-ClinicalEvidenceSummaries_automatic_download.tsv")
}

##########################
###  Load CIViC-Files
###  CIViC: downloaded from https://civic.genome.wustl.edu/releases
##########################

if(!file_exists(path_civic)){
  download_civic_nightly(path_civic, "nightly-ClinicalEvidenceSummaries_automatic_download.tsv")
}

filenames <- list.files(path=path_civic)
filename_civic <- filenames[filenames != "meta.txt"][1]
civic_file <- paste(path_civic, filename_civic, sep="")
civic <- read.delim(civic_file, header=T, stringsAsFactors=F, sep="\t", quote="")

for (c in grep("citation_id", colnames(civic))){
  colnames(civic)[c] = gsub("citation_id", "pubmed_id", colnames(civic)[c])
}

civic <- civic[which(civic$evidence_type=="Predictive"),]

civic_agg  = aggregate(civic, by = list(civic$gene,civic$disease,civic$drugs,civic$evidence_level,civic$clinical_significance),
                       FUN = function(X) paste(unique(X), collapse=", "))[,6:41]

civic = civic_agg

write.table(civic, file=paste(output_path, "CIViC.csv", sep=""),col.names=T,row.names=F,sep="\t")

reader = file(paste(path_civic, "meta.txt", sep=""))
civic_lines = readLines(reader)
close(reader)

##########################
###  Load GDKD-Files
###  GDKD can be dowloaded from https://www.synapse.org/#!Synapse:syn2370773/files/
##########################

if(file_exists(path_gdkd)){
  filename_gdkd <- list.files(path=path_gdkd)[1]
  gdkd_file <- paste(path_gdkd, list.files(path=path_gdkd)[1], sep="")
  gdkd <- read.xlsx(gdkd_file, sheet = 1, colNames = TRUE, cols = c(1:45))
  
  # Small changes
  gdkd$Gene  = gsub(" ", "", gdkd$Gene)
  for (c in colnames(gdkd)){
    gdkd[,c] <- gsub("^ | $", "", gdkd[,c])
  }
  
  for (c in grep("PMID",colnames(gdkd))){
    gdkd[,c] <- gsub("\\^|_","",gdkd[,c])
  }
  
  # Aggregate databases (join variants with same evidence of: disease, association, drug, evidence level)
  gdkd_agg <- aggregate(gdkd, by = list(gdkd$Disease,gdkd$Gene,gdkd$Description,gdkd$Association_1,gdkd$Therapeutic.context_1),
                         FUN = function(X) paste(unique(X), collapse=", "))[,6:50]
  
  gdkd <- gdkd_agg
  
  write.table(gdkd, file=paste(output_path, "GDKD.csv", sep=""),col.names=T,row.names=F,sep="\t")
  
  reader = file(paste(path_gdkd, "meta.txt", sep=""))
  gdkd_lines = readLines(reader)
  close(reader)
  read_flag_gdkd = TRUE
}

##########################
###  Load OncoKB-Files
##########################

if(file_exists(path_oncokb)){
  filename_oncokb <- list.files(path=path_oncokb)[1]
  oncokb_file <- paste(path_oncokb, "/", list.files(path=path_oncokb)[1], sep="")
  oncokb <- read.delim(oncokb_file, header=T, stringsAsFactors = F, sep="\t")
  
  # Rearrange OncoKB to resemble CIVIC
  oncokb = oncokb[c(4,3,5,7,9,8,10,11,1,2)]
  colnames(oncokb) = c("Gene", "Entrez_ID", "Variant", "Disease", "Drugs", "Evidence Level", "Pubmed ID (for drug)", "Abstract (for drug)", "Isoform", "RefSeq")
  
  oncokb_agg = aggregate(oncokb, by = list(oncokb$Gene,oncokb$Disease,oncokb$Drugs,oncokb$`Evidence Level`),
                         FUN = function(X) paste(unique(X), collapse=", "))[,5:14]
  
  oncokb = oncokb_agg
  
  write.csv2(oncokb, file=paste(output_path, "OncoKB.csv", sep=""),row.names=F)

  reader = file(paste(path_oncokb, "meta.txt", sep=""))
  oncokb_lines = readLines(reader)
  close(reader)
  read_flag_oncokb = TRUE
}

##########################
###  Create Data/Version-Info-File
##########################

split = strsplit(civic_lines, ":")
date = split[[1]][2]
version = split[[2]][2]

ver = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = FALSE)

civic_ver = c("CIViC.csv", civic_file, "https://civic.genome.wustl.edu/releases", date, version, "Griffith et al. Nat Genet (2017)", "10.1038/ng.3774")
ver = rbind(ver, civic_ver)

if(read_flag_gdkd){
  split = strsplit(gdkd_lines, ":")
  date = split[[1]][2]
  version = split[[2]][2]
  
  gdkd_ver = c("GDKD.csv", gdkd_file, "https://www.synapse.org/#!Synapse:syn2370773", date, version,"Dienstmann et al. Cancer discovery 5.2 (2015) v20.0", "10.1158/2159-8290.CD-14-1118")
  ver = rbind(ver, gdkd_ver)
}

if(read_flag_oncokb){
  split = strsplit(oncokb_lines, ":")
  date = split[[1]][2]
  version = split[[2]][2]
  
  oncokb_ver = c("OncoKB.csv", oncokb_file, "https://www.oncokb.org/", date, version,"Chakravarty et al. JCO Precis Oncol. (2017)", "10.1200/PO.17.00011")
  ver = rbind(ver, oncokb_ver)
}

ver = setNames(ver, c("Database", "Downloaded file", "URL", "Date", "Version", "Publication", "DOI"))


write.csv(ver, file=paste(output_path,"MTB_Databases_versions.csv",sep=""), row.names = FALSE)


# write table with version information here:
#gdkd_ver = c("GDKD.csv", gdkd_file, "https://www.synapse.org/#!Synapse:syn2370773", "Dienstmann et al. Cancer discovery 5.2 (2015) v20.0", "10.1158/2159-8290.CD-14-1118")
#civic_ver = c("CIViC.csv", civic_file, "https://civic.genome.wustl.edu/releases", "Griffith et al. Nat Genet (2017)", "10.1038/ng.3774")
#target_ver = c("TARGET_MERIC.csv", "", "http://archive.broadinstitute.org/cancer/cga/target", "Van Allen et al. Nat Med (2014) v3 & Meric-Bernstam et al. J Natl Cancer Inst.(2015)", "10.1038/nm.3559 10.1093/jnci/djv098") 
#oncokb_ver = c("OncoKB.csv", oncokb_file, "https://www.oncokb.org/", "Chakravarty et al. JCO Precis Oncol. (2017)", "10.1200/PO.17.00011")

#ver = data.frame(matrix(ncol = 5, nrow = 0), stringsAsFactors = FALSE)
#ver = rbind(ver, gdkd_ver, civic_ver, target_ver, oncokb_ver)
#ver = setNames(ver, c("Database", "Downloaded file", "URL", "Publication", "DOI"))
#write.csv(ver, file="MTB_Databases_versions.csv", row.names = FALSE)

#onco_test = read.csv2("OncoKB.csv",header=T, sep=";")
