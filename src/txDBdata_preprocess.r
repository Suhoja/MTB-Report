library(GenomicFeatures)
options(timeout=1200)

download_GFF_fromURI<-function(GFFURI,save_path){
  #' Base database load function
  #'
  gff_File <- paste(save_path , basename(GFFURI), sep = "")
  if(!file.exists(gff_File)){
    download.file(GFFURI,gff_File)
  }
  return(gff_File)
}

TxDB_generate_from_URI <- function(GFFURI,save_path){
  #' Creates databases from a URI and saves it to a data path
  #'
  gff_File<-download_GFF_fromURI(GFFURI,save_path)
  txDBtempPath<-dirname(gff_File)
  txdbName <- tools::file_path_sans_ext(gff_File)
  txDBPath <-paste(txdbName,"txdb.sqlite",sep = "_")
  if(!file.exists(txDBPath)){
    txdb <- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
    seqlevelsStyle(txdb) <- "UCSC"
    saveDb(txdb, file=txDBPath)
  }
  return(txDBPath)
}

GFFURI_hg38<-"ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.109.gff3.gz" #changed from 100 to 109
GFFURI_hg19<- "ftp://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.gff3.gz"
save_path <- "data/"
txDBPath_hg19 <- TxDB_generate_from_URI(GFFURI_hg19,save_path)
txDBPath_hg38 <- TxDB_generate_from_URI(GFFURI_hg38,save_path)
