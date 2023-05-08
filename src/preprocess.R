
packages <- c("timeSeries", "timeDate", "openxlsx")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages, repos="https://ftp.fau.de/cran/")
  BiocManager::install("openxlsx")
}

library(timeSeries)
library(timeDate)
library(openxlsx)

#setwd("/mnt/c/Users/kevin/OneDrive/Desktop/Git-Repositorys/mtb-report/data")
setwd("/opt/MTB/data/")
print(getwd())

civic_file <- "nightly-ClinicalEvidenceSummaries_automatic.tsv"
gdkd_file <- "Knowledge_database_v20.0.xlsx"
oncokb_file <- "allActionableVariants.txt"

civic_output_file <- "CIViC.csv"
gdkd_output_file <- "GDKD.csv"
oncokb_output_file <- "OncoKB.csv"

civic_meta_file <- "meta.txt"

##########################
###  Load CIViC-Files
###  CIViC: downloaded from https://civic.genome.wustl.edu/releases
##########################

# If there is no CIVIC-File download the nightly-version and add 
# a meta-file with date and version
read_flag_civic <- FALSE
if(!file.exists(civic_output_file)){
  print('Create CIViC')
  print("No database file for CIViC found. Downloading CIViC...")
  print(civic_output_file)
  url <- "https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv"
  download.file(url, civic_file)

  # Meta file generation
  meta_civic_path = paste(output_path, "meta.txt", sep="")
  writer <- file(meta_civic_path)
  date_line <- paste("date", Sys.Date(), sep=":")
  version_line <- paste("version", "1", sep=":")
  writeLines(c(date_line, version_line), writer)
  close(writer)
  read_flag_civic <- TRUE

  # Modify the CIVIC-Table
  civic <- read.delim(civic_file, header=T, stringsAsFactors=F, sep="\t", quote="")

  # Rename 'citation_id' to 'pubmed_id'
  for (c in grep("citation_id", colnames(civic))){
    colnames(civic)[c] = gsub("citation_id", "pubmed_id", colnames(civic)[c])
  }

  # Filter evidence column for rows with Predictive
  civic <- civic[which(civic$evidence_type=="Predictive"),] ## seems like only predicitve rows contain drugs
  civic_agg  = aggregate(civic, by = list(civic$gene,civic$disease,civic$drugs,civic$evidence_level,civic$clinical_significance),
                        FUN = function(X) paste(unique(X), collapse=", "))[,6:41]
  civic = civic_agg
  civic$variant = gsub("OVEREXPRESSION", "Amplification (Overexpr.)", civic$variant, ignore.case = TRUE)
  civic$variant = gsub("UNDEREXPRESSION", "Deletion (Underexpr.)", civic$variant, ignore.case = TRUE)
  civic$variant = gsub("EXPRESSION", "Amplification (Expr.)", civic$variant, ignore.case = TRUE)
  write.table(civic, file=civic_output_file,col.names=T,row.names=F,sep="\t")

  # Read Civic meta file
  reader = file(civic_meta_file)
  civic_lines = readLines(reader)
  close(reader)
}

##########################
###  Load GDKD-Files
###  GDKD can be downloaded from https://www.synapse.org/#!Synapse:syn2370773/files/
##########################

read_flag_gdkd <- FALSE

if(!file.exists(gdkd_output_file)){
  print('Create GDKD')
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
  write.table(gdkd, file=gdkd_output_file,col.names=T,row.names=F,sep="\t")
  
  #reader = file(paste(path_gdkd, "meta.txt", sep=""))
  #gdkd_lines = readLines(reader)
  #close(reader)
  read_flag_gdkd = TRUE
}

##########################
###  Load OncoKB-Files
##########################

read_flag_oncokb = FALSE

#if(!file.exists(oncokb_output_file)){
#  print('Create OncoKB')
#  oncokb <- read.delim(oncokb_file, header=T, stringsAsFactors = F, sep="\t")
#  
#  # Rearrange OncoKB to resemble CIVIC
#  oncokb = oncokb[c(4,3,5,7,9,8,10,11,1,2)]
#  colnames(oncokb) = c("Gene", "Entrez_ID", "Variant", "Disease", "Drugs", "Evidence Level", "Pubmed ID (for drug)", "Abstract (for drug)", "Isoform", "RefSeq")
#  
#  oncokb_agg = aggregate(oncokb, by = list(oncokb$Gene,oncokb$Disease,oncokb$Drugs,oncokb$`Evidence Level`),
#                         FUN = function(X) paste(unique(X), collapse=", "))[,5:14]
#  oncokb = oncokb_agg
#  write.csv2(oncokb, file=oncokb_output_file,row.names=F)
#  #reader = file(paste(path_oncokb, "meta.txt", sep=""))
#  #oncokb_lines = readLines(reader)
#  #close(reader)
#  read_flag_oncokb = TRUE
#}

##########################
###  Load TARGET_MERIC Files
##########################


##########################
###  Create Meta-Info-File
##########################

#foo = civic_lines[1]
#x = strsplit(foo, ":")

#ver = data.frame(matrix(ncol = 7, nrow = 0), stringsAsFactors = FALSE)

#civic_ver = c("CIViC.csv", civic_file, "https://civic.genome.wustl.edu/releases", "Griffith et al. Nat Genet (2017)", "10.1038/ng.3774")
#ver = rbind(ver, civic_ver)

#if(read_flag_gdkd){
#  gdkd_ver = c("GDKD.csv", gdkd_file, "https://www.synapse.org/#!Synapse:syn2370773", "Dienstmann et al. Cancer discovery 5.2 (2015) v20.0", "10.1158/2159-8290.CD-14-1118")
#  ver = rbind(ver, gdkd_ver)
#}

#if(read_flag_oncokb){
#  oncokb_ver = c("OncoKB.csv", oncokb_file, "https://www.oncokb.org/", "Chakravarty et al. JCO Precis Oncol. (2017)", "10.1200/PO.17.00011")
#  ver = rbind(ver, oncokb_ver)
#}

#target_ver = c("TARGET_MERIC.csv", "", "http://archive.broadinstitute.org/cancer/cga/target", "Van Allen et al. Nat Med (2014) v3 & Meric-Bernstam et al. J Natl Cancer Inst.(2015)", "10.1038/nm.3559 10.1093/jnci/djv098") 

#ver = rbind(ver, gdkd_ver, civic_ver, target_ver, oncokb_ver)
#ver = setNames(ver, c("Database", "Downloaded file", "URL", "Date", "Version", "Publication", "DOI"))
#ver = setNames(ver, c("Database", "Downloaded file", "URL", "Publication", "DOI"))
#write.csv(ver, file="MTB_Databases_versions.csv", row.names = FALSE)


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

