# Checks if the Genom Librarys are already installed or need to be installed

#check if BiocManager is already installed
if(!require("BiocManager")){
  print("Biocmanager not installed!");
  print("Install BiocManager:");
  install.packages("BiocManager");
} else{
  print("BiocManager already installed.")
}

#check if genom-data already installed
if(!require("BSgenome.Hsapiens.UCSC.hg19")){
  print("BSgenome.Hsapiens.UCSC.hg19 not installed!");
  print("Install BSgenome.Hsapiens.UCSC.hg19:");
  #set timeout to 600 seconds due long download time
  options(timeout=600); 
  BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19'));
} else{
  print("BSgenome.Hsapiens.UCSC.hg19 already installed.")
}

#check if genom-data already installed
if(!require("BSgenome.Hsapiens.UCSC.hg38")){
  print("BSgenome.Hsapiens.UCSC.hg38 not installed!");
  print("Install BSgenome.Hsapiens.UCSC.hg38:");
  #set timeout to 600 seconds due long download time
  options(timeout=600); 
  BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38'));
} else{
  print("BSgenome.Hsapiens.UCSC.hg38 already installed.")
}

