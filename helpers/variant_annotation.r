suppressPackageStartupMessages(library(maftools, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(VariantAnnotation, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(GenomicFeatures, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(org.Hs.eg.db, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tidyr, quietly=TRUE, warn.conflicts = FALSE))


#############################
# #Data load method
# download_GFF_fromURI<-function(GFFURI,save_path){
#   
#   # path = system.file(package="CoordinateTest", "extdata")
#   # if(path==""){
#   #   path="~/Downloads"
#   # }
#   gff_File <<- paste(save_path , "/", basename(GFFURI), sep = "")
#   #gfffilename <- basename(GFFURI)
#   if(!file.exists(gff_File)){
#     download.file(GFFURI,gff_File)
#   }
#   
#   return(gff_File)
# }
# GFFTXDBDataload <- function(gff_File){
# 
#   if(exists("txDBPath") & file.exists("txDBPath")  ){
#     txdb <<- loadDb(txDBPath)
#     seqlevelsStyle(txdb) <- "UCSC"
#     txDBPath<<-txDBPath
#   } else if(exists("txDBPath") & !file.exists("txDBPath")){
#     txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
#     seqlevelsStyle(txdb) <- "UCSC"
#     txdbName <- tools::file_path_sans_ext(gff_File)
#     txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
#     saveDb(txdb, file=txDBPath)
# 
#   } else if(! exists("txDBPath") ){
#     txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
#     seqlevelsStyle(txdb) <- "UCSC"
#     txdbName <- tools::file_path_sans_ext(gff_File)
#     txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
#     saveDb(txdb, file=txDBPath)
#   }
#   return(txDBPath)
# }
# ###Saving and Loading a TxDb Object
# #gff_File = "~/Downloads/GIT/Tmp_Data/Homo_sapiens.GRCh38.99.gff3.gz"
# 
# TxDB_generate_from_URI <- function(GFFURI,save_path){
#   gff_File<-download_GFF_fromURI(GFFURI,save_path)
#   txdb <<- makeTxDbFromGFF(gff_File,organism="Homo sapiens")
#   seqlevelsStyle(txdb) <- "UCSC"
#   txDBtempPath<-dirname(gff_File)
#   txdbName <- tools::file_path_sans_ext(gff_File)
#   txDBPath <<-paste(txdbName,"txdb.sqlite",sep = "_")
#   saveDb(txdb, file=txDBPath)
#   #  txdb <- loadDb(txDBPath)
#   return(txDBPath)
# }

#txDBPath = "~/Downloads/GIT/Tmp_Data/Homo_sapiens.GRCh38.99.gff3_txdb.sqlite"
#txDBPath<-TxDB_generate(gff_File)

################################################
VCF_to_SNV <- function(Input_Vcf,txDBPath){
    #txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  #txDBPath ="gff_annotation.sqlite"
  txdb <- loadDb(txDBPath)
  vcf <- readVcf(Input_Vcf)

  print("t3")
  seqlevelsStyle(vcf) <- "UCSC"
  seqlevelsStyle(txdb) <- "UCSC"

  #print("check seq stuff")
  #print(intersect(seqlevels(vcf), seqlevels(txdb))) #Returns chr1 through chr22 and chrX, chrY, chrM. This is as it should be. However, sarek-files are empty in predictcoding.

  seqlevels(vcf,pruning.mode="coarse")<-intersect(seqlevels(vcf), seqlevels(txdb))

  coding <- predictCoding(vcf, txdb, seqSource=Hsapiens)
  coding$GENEID <- gsub("gene:","",coding$GENEID) #Removes gene: which for some reson appears in GENEID after input
  
  #print("Check result of predictCoding() in VCF_to_SNV colCoding")  #Checked, says unspecified genome..
  #print(head(coding))

  colCoding <-mcols(coding)

  colCoding$PROTEINLOC <- sapply(colCoding$PROTEINLOC, function(x) paste(unlist(x), collapse=", ")) #PROTEINLOC is of type IntegerList - must be split

  Protein_Change <- paste(colCoding$REFAA,colCoding$PROTEINLOC,colCoding$VARAA,sep="")
  #print("Protein_change paste runs ok")
  #replace item last char as '=' which amino acid code did not changed.
  Protein_Change[colCoding$REFAA==colCoding$VARAA] <- paste(colCoding$REFAA,colCoding$PROTEINLOC,"=",sep="")[colCoding$REFAA==colCoding$VARAA]
  #print("Protein_change subset runs ok")
  Hugo_Symbol <- colCoding$GENEID #These appear to actually be ensembl_ids?

  annots <- select(org.Hs.eg.db, keys=Hugo_Symbol, columns="SYMBOL", keytype="ENSEMBL") #Added code to convert ensembl to gene symbols
  annots <- annots[match(Hugo_Symbol, annots$ENSEMBL),]
  Hugo_Symbol <- annots$SYMBOL

  #print("Check Hugo_symbol in VCF_to_SNV") #19.4. issue between here and check coding
  #print(head(Hugo_Symbol))

  Variant_Classification<-colCoding$CONSEQUENCE

  #calculate insertion or deletion by length difference of varAllele and REF
  Indel_status<- lengths(colCoding$varAllele)-lengths(colCoding$REF)
  SNVTable <- data.frame(Hugo_Symbol,Variant_Classification,Protein_Change,Indel_status)
  #print("Check SNVTable")
  #print(head(SNVTable))
  #print(nrow(SNVTable))                          #At this point the table has data for mutect2, but only 5 rows
  #                                               This could be explained by the data having only 14 mutations...
  ####convert variant_classification
  print("t4")                                    # !!!!  When using mutect2 output this is the last output before erroring out  !!!!
    ###remove silent mutation
  SNVTable<-SNVTable[!(SNVTable[,"Variant_Classification"]=="synonymous"),]

  if (nrow(SNVTable) == 0) {
    stop("SNVTable in variant_annotation.r is empty after filtering silent mutations. Check input data.")
  }

  SNVTable$Variant_Classification<- as.character(SNVTable$Variant_Classification)
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="nonsense"] <-  "Nonsense_Mutation"
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="nonsynonymous"] <-  "Missense_Mutation"
  SNVTable$Variant_Classification[SNVTable$Indel_status > 0 & SNVTable$Variant_Classification=="frameshift"] <- "Frame_Shift_Ins"
  SNVTable$Variant_Classification[SNVTable$Indel_status < 0 & SNVTable$Variant_Classification=="frameshift"] <- "Frame_Shift_Del"
  SNVTable$Variant_Classification[SNVTable$Indel_status > 0 & SNVTable$Variant_Classification=="missense"] <- "In_Frame_Ins"
  SNVTable$Variant_Classification[SNVTable$Indel_status < 0 & SNVTable$Variant_Classification=="missense"] <- "In_Frame_Del"
  SNVTable$Indel_status <- NULL

  ##################################
  SNVTable = unique(SNVTable)

}


Coordinate_to_SNV <- function(Coordinate_input,txDBPath){
  #txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
 # txDBPath ="gff_annotation.sqlite"
  txdb <- loadDb(txDBPath)
#  Coordinate_input_path

    Coordinate_input_Table <- tryCatch(read.table(Coordinate_input,sep="\t", header=TRUE),error = function(error_condition) {
      return(data.frame())
    })


  #Coordinate_input_Table <- read.table(Coordinate_input,sep="\t", header=TRUE)
  emptyCoordTable <-data.frame(Hugo_Symbol = character(),
                               Variant_Classification = character(),
                               Protein_Change = character())
  if(nrow(Coordinate_input_Table)==0){
    return(emptyCoordTable)
  }
  if(!("Chromosome" %in% names(Coordinate_input_Table)) ){
    print("Input file missing chromosome column")
    return(emptyCoordTable)
  }
  if(!("start" %in% names(Coordinate_input_Table))){
    print("Input file missing start column")
    return(emptyCoordTable)
  }
  if(!("end" %in% names(Coordinate_input_Table))){
    print("Input file missing end column")
    return(emptyCoordTable)
  }
  if(!("refAllele" %in% names(Coordinate_input_Table))){
    print("Input file missing refAllele column")
    return(emptyCoordTable)
  }
  if(!("varAllele" %in% names(Coordinate_input_Table))){
    print("Input file missing varAllele column")
    return(emptyCoordTable)
  }


  if("strand" %in%  colnames(Coordinate_input_Table)){
    Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                  ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                  strand = Coordinate_input_Table$strand,
                                  REF = Coordinate_input_Table$refAllele)
  } else {
    Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                  ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                  REF = Coordinate_input_Table$refAllele)
  }
  seqlevelsStyle(txdb) <- "UCSC"
  seqlevelsStyle(Coordinate_Granges) <- "UCSC"

  VarAllele<-DNAStringSet(Coordinate_input_Table$varAllele)

  coding <- predictCoding(Coordinate_Granges, txdb, Hsapiens,VarAllele)
  colCoding <-mcols(coding)

  Protein_Change <- paste(colCoding$REFAA,colCoding$PROTEINLOC,colCoding$VARAA,sep="")
  #replace item last char as '=' which amino acid code did not changed.
  Protein_Change[colCoding$REFAA==colCoding$VARAA] <- paste(colCoding$REFAA,colCoding$PROTEINLOC,"=",sep="")[colCoding$REFAA==colCoding$VARAA]

  Hugo_Symbol <- colCoding$GENEID
  Variant_Classification<-colCoding$CONSEQUENCE
  #calculate insertion or deletion by length difference of varAllele and REF
  Indel_status<- lengths(colCoding$varAllele)-lengths(colCoding$REF)
  SNVTable<- data.frame(Hugo_Symbol,Variant_Classification,Protein_Change,Indel_status)
  ####convert variant_classification

    ###remove silent mutation
  SNVTable<-SNVTable[!(SNVTable[,"Variant_Classification"]=="synonymous"),]

    #replace consequence name
  SNVTable$Variant_Classification<- as.character(SNVTable$Variant_Classification)
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="nonsynonymous"] <-  "missense"
  SNVTable$Variant_Classification[SNVTable$Variant_Classification=="frameshift"] <- "frame_shift"
  SNVTable$Variant_Classification[SNVTable$Indel_status > 0 & SNVTable$Variant_Classification=="missense"] <- "ins"
  SNVTable$Variant_Classification[SNVTable$Indel_status < 0 & SNVTable$Variant_Classification=="missense"] <- "del"
  SNVTable$Indel_status <- NULL

  ##################################
  SNVTable = unique(SNVTable)
  return(SNVTable)

}

#snvTest<-Coordinate_to_SNV(Coordinate_input_path,txDBPath)
#snvTest<-VCF_to_SNV(vcfTest,txDBPath)


getSNVfromMaf <- function(mafInput) {
  maftable<-mafInput@data


  ######### unify synonyms############
  if("HGVSp_Short" %in% names(maftable)&"Protein_Change" %in% names(maftable))
  {
    print("Problem: duplicate protein change column ")
  } else if("HGVSp_Short" %in% names(maftable)&!("Protein_Change" %in% names(maftable)))
  {
    protein_Change_Index <- which(names(maftable)=="HGVSp_Short")
    names(maftable)[protein_Change_Index] <- "Protein_Change"
  }
  ##########################
  SNV<-maftable

  #  SNV<-maftable[maftable$Variant_Type!="SNP",]    #remove SNP
  SNV<-as.data.frame(SNV)                                     #convert data.table to data.frame
  SNV                = SNV[SNV$Variant_Classification!='Silent',]
  SNV$Protein_Change = gsub("p.","",SNV[,"Protein_Change"])
  SNV                = SNV[,c("Hugo_Symbol","Variant_Classification","Protein_Change")]
  SNV                = SNV[which(SNV[,"Protein_Change"] != ""),]
  SNV                = unique(SNV)
  return(SNV)
}

###### generate SNV_table from maf file
MAF_to_SNVTable<-function(MafInputFile,txDBPath){
  mafInput = read.maf(MafInputFile)
  mafFields<-getFields(mafInput)
  if("HGVSp_Short" %in% mafFields |"Protein_Change" %in% mafFields  )
  {
    SNVTable<-getSNVfromMaf(mafInput)
    SNVTable = unique(SNVTable)
    return(SNVTable)
  } else {
    maftable<-mafInput@data
    Chromosome <- sub("chr","", maftable$Chromosome )
    #maftable[,"Start_Position"]

    Coordinate_input_Table<-data.frame(Chromosome,maftable$Start_Position,maftable$End_Position,maftable$Reference_Allele,
                              maftable$Tumor_Seq_Allele2,maftable$Strand)
    Coordinate_input_Table <-unique(Coordinate_input_Table)
    colnames(Coordinate_input_Table) <- c("Chromosome",	"start","end","refAllele","varAllele","strand")
    txdb <- loadDb(txDBPath)
    if("strand" %in%  colnames(Coordinate_input_Table)){
      Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                    ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                    strand = Coordinate_input_Table$strand,
                                    REF = Coordinate_input_Table$refAllele)
    } else {
      Coordinate_Granges <- GRanges(seqnames = Rle(Coordinate_input_Table$Chromosome),
                                    ranges = IRanges(start = Coordinate_input_Table$start, end = Coordinate_input_Table$end),
                                    REF = Coordinate_input_Table$refAllele)
    }
    seqlevelsStyle(txdb) <- "UCSC"
    seqlevelsStyle(Coordinate_Granges) <- "UCSC"

    VarAllele<-DNAStringSet(Coordinate_input_Table$varAllele)

    coding <- predictCoding(Coordinate_Granges, txdb, Hsapiens,VarAllele)
    colCoding <-mcols(coding)

    Protein_Change <- paste(colCoding$REFAA,colCoding$PROTEINLOC,colCoding$VARAA,sep="")
    #replace item last char as '=' which amino acid code did not changed.
    Protein_Change[colCoding$REFAA==colCoding$VARAA] <- paste(colCoding$REFAA,colCoding$PROTEINLOC,"=",sep="")[colCoding$REFAA==colCoding$VARAA]

    SNVTable<- data.frame(maftable$Hugo_Symbol,maftable$Variant_Classification,Protein_Change)
    SNVTable = unique(SNVTable)
    return(SNVTable)
  }
}


####### detect the extension and call the corresponding functions ##########
metadata_to_SNV <- function(metafile,txDBPath){
  if(endsWith(metafile,".maf")){
    SNVTable <- MAF_to_SNVTable(metafile,txDBPath)
    return(SNVTable)
  }
  else if (endsWith(metafile,".txt")  ) {
    SNVTable<-Coordinate_to_SNV(metafile,txDBPath)
    return(SNVTable)
  }
  else if(endsWith(metafile,".vcf")){
    SNVTable<-VCF_to_SNV(metafile,txDBPath)
    return(SNVTable)
  }
  else {
    print("incorrect file type, the input file can only be '.maf, .vcf and coordinate .txt files'")
    return(NULL)
  }
}
#####################
# TestDataload<-function(testDataFilesDir){
#   path="~/Downloads"
#   subDir = "variant_annotation_testData"
#   testDataPath<-file.path(path, subDir)
#   dir.create(testDataPath, showWarnings = FALSE)
# 
#   for(i in seq_along(testDataFilesDir)){
#     print(basename(testDataFilesDir[i]))
#     testDataFile<-paste(testDataPath,basename(testDataFilesDir[i]),sep = "/")
#     print(testDataFile)
#     download.file(testDataFilesDir[i], testDataFile, mode="wb")
#   }
# 
# }
