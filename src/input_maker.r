### INPUT MAKER ###
### Charlotte Höltermann ###
# makes: CNV.csv and SNV.csv files that can be given to localmaster.r as they are, as arguments from the command line.
# SPECIFICALLY MADE FOR the supplementary data files from:
#Nintog, S.Y., Yoshida, N., Christie, A.L. et al. Targetable vulnerabilities in T- and NK-cell lymphomas identified through preclinical models. Nat Commun 9, 2024 (2018) doi:10.1038/s41467-018-04356-9
#########
#Gene names must be Hugo Symbols.
#Variant Classification comprises the following levels: 
#Frame_Shift_Del, 
#Frame_Shift_Ins, 
#In_Frame_Del, 
#In_Frame_Ins, 
#Missense_Mutation, 
#Nonsense_Mutation, 
#Silent, 
#Splice_Site, 
#Translation_Start_Site, 
#Nonstop_Mutation, 
#RNA, 
#Targeted_Region
#
#Protein Change must be e.g. T790M
#CN alteration must be one of the foollowing levels: {amplification, deletion} #bug here?
#install.packages("openxlsx")
library(openxlsx)

#args <- commandArgs(trailingOnly=TRUE)
# arg1: 2,7,8,9
# arg2: path to file to be read
#arg_suppnum = as.integer(args[1])
#arg_file = toString(args[2])
#destination = toString(args[3])

setwd("/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw")
#destination = "Input"
destination = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Input"
#destination = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Input"
#arg_file = "Paper_research_nicole/MTB_Paper_Summary_v6.xlsx" 
#arg_suppnum = 1
#arg_file = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Ng_et_al_2018_Sequencing_Data_T-NHL/SupplementaryDataFile2-RNASeqfusions.xlsx" 
#arg_suppnum = 2
#arg_file = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Ng_et_al_2018_Sequencing_Data_T-NHL/SupplementaryDataFile7-CNV.xlsx" 
#arg_suppnum = 7
#arg_file = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Ng_et_al_2018_Sequencing_Data_T-NHL/SupplementaryDataFile8-mutationscelllines.xlsx" 
#arg_suppnum = 8
#arg_file = "/sybig/home/chh/Dokumente/TCL_Project_Koch_Dönitz/Input_MTB_raw/Ng_et_al_2018_Sequencing_Data_T-NHL/SupplementaryDataFile9-mutationsPDXs.xlsx" 
#arg_suppnum = 9
###
make_variant_classification_ok = function(c1, fr, ra){
  fr = fr
  range = 1:ra
  for (r in range){
    ro <- toString(fr[r, c1])
    # Convert strings meaning the same into the desired format
    if (ro == "Frame_Shift_Del" | ro == "Frame_Shift_Ins" | ro == "In_Frame_Del" | ro == "In_Frame_Ins" | ro == "Missense_Mutation" |ro == "Nonsense_Mutation" |ro == "Silent" |ro == "Splice_Site" |ro == "Translation_Start_Site" |ro == "Nonstop_Mutation" |ro == "RNA" |ro == "Targeted_Region" | ro == "De_novo_Start_InFrame" | ro == "De_novo_Start_OutOfFrame"){
    } else if (ro == "Missense"){
      fr[r, c1] <- gsub(ro, "Missense_Mutation", fr[r, c1])
    } else if (ro == "Frameshift"){
      fr[r, c1] <- gsub(ro, "Frame_Shift", fr[r, c1])
    } else if (ro == "Nonsense"){
      fr[r, c1] <- gsub(ro, "Nonsense_Mutation", fr[r, c1])
    } else {
      #warning(paste("There is an unknown description: ", toString(ro)))
      }
  }
  return(fr)
}
### functions for the different types of .xlsx files in the supplementary data ###
# Supplementary file 2 # FUSIONS
r_supp2 <- function(){
  dest = paste(destination ,"/Ng_et_al_2018", sep ="")
  dir.create(dest, showWarnings = FALSE)
  fr = readWorkbook(arg_file, sheet = 1, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
  fr = fr[c(4:nrow(fr)), c(1,5,8)]
  fr = setNames(fr, c("Sample","Gene1","Gene2")) #Gene1  = left, Gene2= right gene
  for (i in 1:nrow(fr)){
    fr[i,1] = gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", fr[i,1])
  }
  s = unique(fr[, 1]) #character vector with unique "Sample" names
  for (e in s){
    e = toString(e)
    fr_sub <- subset(fr, fr[,1]==e)
    fr_sub = fr_sub[, c(2,3)] #subset the frame as required by MTB
    e = gsub(" ", "", e)
    fr_sub = setNames(fr_sub, c(paste(e, "_Gene1", sep = ""), paste(e, "_Gene2", sep = "")))
    fr_sub = unique(fr_sub)
    write.csv(fr_sub, paste(dest, "/", "Ng_et_al_2018_", "Supp_2_Fusion_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
  }
}
# Supplementary file 7 # CNV
r_supp7 <- function(){
  dest = paste(destination ,"/Ng_et_al_2018", sep = "")
  dir.create(dest, showWarnings = FALSE)
  fr = readWorkbook(arg_file, sheet = 1, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
  fr = fr[c(2:nrow(fr)), c(14, 2, 12)]
  fr = setNames(fr, c("Sample","Gene","CN_type"))
  fr = subset.data.frame(fr, fr$CN_type == "+" | fr$CN_type == "-") #which are not 0 in the segment_call
  for (n in 1:nrow(fr)){
    el = toString(fr[n, 3])
    if (el == "+"){fr[n, 3] = "amplification"}
    else if (el == "-"){fr[n, 3] = "deletion"}
  }
  for (i in 1:nrow(fr)){
    fr[i,1] = gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", fr[i,1])
  }
  s = unique(fr[,1]) #character vector with unique "Sample" names
  for (e in s){
    e = toString(e)
    fr_sub = subset(fr, fr$Sample==e)
    fr_sub = fr_sub[, c(2,3)] #subset the frame as required by MTB
    e = gsub(" ", "", e)
    fr_sub = setNames(fr_sub, c(paste(e, "_Gene", sep = ""), "CN_type"))
    fr_sub = unique(fr_sub)
    write.csv(fr_sub, paste(dest, "/", "Ng_et_al_2018_", "Supp_7_CNV_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
  }
}
# Supplementary file 8 # SNV cell lines
#load the supplementary file as a workbook. Gives access to the excel sheets names which are the cell lines.
r_supp8 <- function(){
  dest = paste(destination ,"/Ng_et_al_2018", sep ="")
  dir.create(dest, showWarnings = FALSE)
  excelfile <- loadWorkbook(arg_file)
  samples = names(excelfile)
  big_frame = data.frame(matrix(ncol = 7, nrow = 0))
  big_frame = setNames(big_frame, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type", "Protein_Change", "Ref_allele", "Tumor_allele"))
  for (i in c(1,2,3,4,5,20,21,22,33)) { #type1
    sa = toString(samples[i])
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    #modifing the frame here to fit the required format
    fr = fr[2:nrow(fr), c(1,8,9,19,10,11)]
    firstcol = c()
    #append the name of the sample, meaning the name of the cell line as a column
    for (n in 1:nrow(fr)){
      firstcol = append(firstcol, sa)
      }
    fr = cbind(firstcol, fr)
    fr = setNames(fr, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type","Protein_Change", "Ref_allele", "Tumor_allele"))
    big_frame = rbind(big_frame, fr)
  }
  for (i in c(6,7,8,9,10,14,15,16,17,18,19,23,26)){ #type 2
    sa = toString(samples[i])
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    #modifing the frame here to fit the required format
    fr = fr[2:nrow(fr), c(1,9,10,20,11,13)]
    firstcol = c()
    for (n in 1:nrow(fr)){
      firstcol = append(firstcol, sa)
    }
    fr = cbind(firstcol, fr)
    fr = setNames(fr, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type","Protein_Change", "Ref_allele", "Tumor_allele"))
    big_frame = rbind(big_frame, fr)
  }
  for (i in c(11)){ #type 3
    sa = toString(samples[i])
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    #modifing the frame here to fit the required format
    fr = fr[2:nrow(fr), c(7,6,6,9,3,4)]
    firstcol = c()
    for (n in 1:nrow(fr)){
      firstcol = append(firstcol, sa)
    }
    fr = cbind(firstcol, fr)
    fr = setNames(fr, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type","Protein_Change", "Ref_allele", "Tumor_allele"))
    big_frame = rbind(big_frame, fr)
  }  
  for (i in c(12,13)){ #type 4
    sa = toString(samples[i])
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    #modifing the frame here to fit the required format
    fr = fr[2:nrow(fr), c(8,7,7,10,4,5)]
    firstcol = c()
    for (n in 1:nrow(fr)){
      firstcol = append(firstcol, sa)
    }
    fr = cbind(firstcol, fr)
    fr = setNames(fr, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type","Protein_Change", "Ref_allele", "Tumor_allele"))
    big_frame = rbind(big_frame, fr)
  } 
  for (i in c(24,25,27,28,29,30,31,32)){ #type 5 
    sa = toString(samples[i])
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    fr = fr[2:nrow(fr), c(10,11,11,11,4,5)]
    #make adjustments for one row containing more than one mutation
    nr = nrow(fr)
    for (r in 1:nr){
      fr[r, ] = gsub("GI=", "", fr[r, ])
      fr[r, ] = gsub("FC=", "", fr[r, ])
      hugo = strsplit(fr[r, 1], ",", fixed = TRUE)
      variant = strsplit(fr[r, 3], ",", fixed = TRUE)
      proteinchange = strsplit(fr[r, 4], ",", fixed = TRUE)
      hugo = hugo[[1]]
      variant = variant[[1]]
      proteinchange = proteinchange[[1]]
      for (e in 1:length(variant)){
        if (grepl("Missense|Nonsense|Deletion|Frameshift|Silent", variant[e]) == TRUE){
        variant[e] = gsub(".*(Missense).*", "Missense", variant[e]) 
        variant[e] = gsub(".*(Nonsense).*", "Nonsense", variant[e]) 
        variant[e] = gsub(".*(Deletion).*", "Deletion", variant[e])
        variant[e] = gsub(".*(Frameshift).*", "Frameshift", variant[e])
        variant[e] = gsub(".*(Silent).*", "Silent", variant[e])
        } else {warning(paste("There is more than Missense, Nonsense, Deletion, Frameshift, Silent in type 5", toString(variant[e])))}
      }
      for (e in 1:length(proteinchange)){
       proteinchange[e] = gsub("_?(Missense|Nonsense|Deletion|Frameshift|Silent)_?", "", proteinchange[e]) 
      }
      for (i in 1:length(hugo)){
        vec = c(toString(hugo[i]), toString(variant[i]), "", toString(proteinchange[i]), toString(fr[r, 5]), toString(fr[r, 6]))
        fr = rbind(fr, vec)
      }
    }
    fr = fr[-c(1:nr), ]
    # tumor allele = fr[r, 6]
    # ref allele = 5
    for (r in 1:nrow(fr)){
      if (nchar(toString(fr[r, 5])) == 1 & nchar(toString(fr[r, 6])) == 1){
        
      } else if (nchar(toString(fr[r, 6])) < nchar(toString(fr[r, 5])) & nchar(toString(fr[r, 6])) == 1){
        #adjust Ins/del synatx to the one used in the other sheets
        fr[r, 5] = gsub("^.", "", fr[r, 5])
        fr[r, 6] = gsub("^.", "-", fr[r, 6])
      } else if (nchar(toString(fr[r, 6])) > nchar(toString(fr[r, 5])) & nchar(toString(fr[r, 5])) == 1){
        fr[r, 6] = gsub("^.", "", fr[r, 6])
        fr[r, 5] = gsub("^.", "-", fr[r, 5])
      } else if (nchar(toString(fr[r, 6])) == nchar(toString(fr[r, 5]))){
          warning("There is something wrong with supp8 nt type 5")
        }
    }
    firstcol = c()
    for (n in 1:nrow(fr)){
      firstcol = append(firstcol, sa)
    }
    fr = cbind(firstcol, fr)
    fr = setNames(fr, c("Sample", "Hugo_Symbol", "Variant_Classification", "Variant_Type","Protein_Change", "Ref_allele", "Tumor_allele"))
    big_frame = rbind(big_frame, fr)
  }
    ######
  for (r in 1:nrow(big_frame)){
    big_frame[r, 5] <- gsub("p.", "", big_frame[r, 5])
    big_frame[r, 1] = gsub("_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "", big_frame[r, 1])
  }
  big_frame = make_variant_classification_ok(3, big_frame, nrow(big_frame))
  #delete the not needed columns from the dataframe by subsetting
  s = unique(big_frame[, 1])
  for (e in s){
    e = toString(e)
    fr_sub = subset(big_frame, big_frame$Sample==e)
    fr_sub = fr_sub[, c(2,3,5)] #subset the frame as required by MTB
    e = gsub(" ", "", e)
    fr_sub = setNames(fr_sub, c(paste(e, "_Hugo_Symbol", sep = ""), "Variant_Classification", "Protein_Change"))
    fr_sub = unique(fr_sub)
    write.csv(fr_sub, paste(dest, "/", "Ng_et_al_2018_", "Supp_8_SNV_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
  }
}

# Supplementary file 9 # SNV PDX
r_supp9 <- function(){
  for (num in c(1,2)){
    dest = paste(destination ,"/Ng_et_al_2018", sep ="")
    dir.create(dest, showWarnings = FALSE)
    fr = readWorkbook(arg_file, sheet = num, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    fr = fr[c(2:nrow(fr)), c(1, 7, 6, 6, 8, 4, 5)]
    fr = setNames(fr, c("PDX", "Hugo_Symbol","Variant_Classification","Variant_Type", "Protein_Change", "Ref_allele", "Tumor_allele"))
    #modifing the frame here to fit the required format
    for (n in 1:nrow(fr)){
      fr[n, 4] <- gsub("p.", "", fr[n, 4])
    }
    fr = make_variant_classification_ok(3, fr, nrow(fr))
    # Separate by PDX/sample name
    # write a .csv for each PDX
    pdx = unique(fr[, 1]) #character vector with unique "Sample" names
    for (e in pdx){
      e = toString(e)
      fr_sub <- subset(fr, fr$PDX==e)
      fr_sub = fr_sub[, c(2,3,5)] #subset the frame as required by MTB
      e = gsub(" ", "", e)
      fr_sub = setNames(fr_sub, c(paste(e, "_Hugo_Symbol", sep = ""), "Variant_Classification", "Protein_Change"))
      fr_sub = unique(fr_sub)
      if (num == 1){
        write.csv(fr_sub, paste(dest, "/", "Ng_et_al_2018_", "Supp_9_SNV_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
      }
      if (num == 2){
        write.csv(fr_sub, paste(dest, "/", "Ng_et_al_2018_", "Supp_9_SNV_", e, "_SNP_filtered", ".csv", sep = ""), row.names = FALSE, quote=FALSE)
      }
    }
  }
}
###
# dir = ? make dir for every publication, include a file with metadata
# destination = paste(destination, dir)

### Paper research nicole (Data received on 16.01.2020)
# contact: Nicole
# paper: 
#Sakata-Yanagimoto et al. 14 
#Yoo et al. 14
#Odejide et al. 14
#Palomero et al. 14
#Crescenzo et al. 15 
#Schatz et al. 15 
#McKinney et al. 17
#Abate et al. 17 
#Manso et al. 17
#Moffitti et al. 17
#Ji et al. 18
#Song et al. 18
#Laginestra et al. 19
#Watatani et al. 19
#Heavican et al. 19 
#Vallois et al. 18
#Nguyen et al. 17
#Yoo et al. 16
#Boddicker et al. 16
r_nicole <- function(){ 
  excelfile <- loadWorkbook(arg_file)
  publications = names(excelfile)
  for (i in c(1,2,4,6,7,9,10,11,12,13,14,15,18,19,23,24,26,27,30,31)) { #type1 + type3 -- SNV
    #make directory for each publication
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    pub = toString(publications[i])
    pub = gsub(" ", "_", pub)
    pub = gsub("\\.", "", pub)
    dir.create(paste("Input/", pub, sep = ""), showWarnings = FALSE)
    #subsetting the frame here to fit the required format
    fr = fr[2:nrow(fr), c(1, 2, 3, 6, 7)]
    fr = setNames(fr, c("Sample", "Subtype", "Hugo_Symbol", "Variant_Classification", "Protein_Change"))
    fr = make_variant_classification_ok(4, fr, nrow(fr))
    for (num in 1:nrow(fr)){
      fr[num, 5] <- gsub("p.", "", fr[num, 5], fixed = TRUE)
    }
    # teilweise mehrere protein changes pro zelle in Excel! das muss ich lösen:
    fr_exp <- data.frame(matrix(ncol = 6, nrow = 0))
    fr_exp <- setNames(fr_exp, c("Sample", "Subtype", "Hugo_Symbol", "Variant_Classification", "Protein_Change", "subset"))
    fr$subset <- rep("NO" , nrow(fr))
    for (num in 1:nrow(fr)){
      p_changes = toString(fr[num, 5])
      p_changes = unlist(strsplit(p_changes, "(\\s|\\c)"))
      if (length(p_changes) > 1){
        fr[num, 6] <- "YES"
        for (ch in p_changes){
          #create new dataframe from current row containing just this change
          #append new line with just this change to df
          newrow = fr[num,]
          newrow$Protein_Change = rep(ch, nrow(newrow))
          fr_exp = rbind(fr_exp, newrow)
        }
      }
    }
    fr <- fr[which(fr$subset == "NO"), ]
    fr <- rbind(fr, fr_exp)
    fr <- fr[, 1:5]
    samples = unique(fr[, 1])
    for (e in samples){
      e = toString(e)
      fr_sub = subset(fr, fr$Sample==e)
      e = gsub(" ", "", e)
      e = gsub("/", "_", e)
      subtype = toString(fr_sub[1,2])
      subtype = gsub(" ", "_", subtype)
      fr_sub = fr_sub[, c(3,4,5)] #subset the frame as required by MTB
      fr_sub = setNames(fr_sub, c(paste(e, "_Hugo_Symbol", sep = ""), "Variant_Classification", "Protein_Change"))
      fr_sub = unique(fr_sub)
      write.csv(fr_sub, paste("Input/", pub, "/", pub, "_SNV_", subtype, "_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
    }
  }
  for (i in c(8,20,21,22,25,28)){ #type2 -- Fusions
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    pub = toString(publications[i])
    pub = gsub(" ", "_", pub)
    pub = gsub("\\.", "", pub)
    dir.create(paste("Input/", pub, sep = ""), showWarnings = FALSE)
    #subsetting the frame here to fit the required format
    fr = fr[2:nrow(fr), c(1, 2, 3, 4)]
    fr = setNames(fr, c("Sample", "Subtype", "Gene1", "Gene2"))
    samples = unique(fr[, 1])
    for (e in samples){
      e = toString(e)
      fr_sub = subset(fr, fr$Sample==e)
      e = gsub(" ", "", e)
      e = gsub("/", "_", e)
      subtype = toString(fr_sub[1,2])
      subtype = gsub(" ", "_", subtype)
      fr_sub = fr_sub[, c(3,4)] #subset the frame as required by MTB
      fr_sub = setNames(fr_sub, c(paste(e, "_Gene1", sep = ""), "Gene2"))
      fr_sub = unique(fr_sub)
      write.csv(fr_sub, paste("Input/", pub, "/", pub, "_Fusion_", subtype, "_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
    }
  }
  for (i in c(29)){ #type4 -- CNV
    fr = readWorkbook(arg_file, sheet = i, startRow = 1, colNames = FALSE, rowNames = FALSE, detectDates = FALSE, skipEmptyRows = FALSE ,skipEmptyCols = TRUE, check.names = FALSE, sep.names = "",namedRegion = NULL, na.strings = "NA", fillMergedCells = FALSE)
    pub = toString(publications[i])
    pub = gsub(" ", "_", pub)
    pub = gsub("\\.", "", pub)
    dir.create(paste("Input/", pub, sep = ""), showWarnings = FALSE)
    #subsetting the frame here to fit the required format
    fr = fr[2:nrow(fr), c(1, 2, 5, 12)]
    fr = setNames(fr, c("Sample", "Subtype", "Hugo_Symbol", "CN_type"))
    samples = unique(fr[, 1])
    for (e in samples){
      e = toString(e)
      fr_sub = subset(fr, fr$Sample==e)
      e = gsub(" ", "", e)
      e = gsub("/", "_", e)
      subtype = toString(fr_sub[1,2])
      subtype = gsub(" ", "_", subtype)
      fr_sub = fr_sub[, c(3,4)] #subset the frame as required by MTB
      fr_sub = setNames(fr_sub, c(paste(e, "_Hugo_Symbol", sep = ""), "CN_type"))
      fr_sub = unique(fr_sub)
      write.csv(fr_sub, paste("Input/",pub, "/", pub, "_CNV_", subtype, "_", e, ".csv", sep = ""), row.names = FALSE, quote=FALSE)
    }
  }
  # 3,5
}

### execute ###
if (arg_suppnum == 2){r_supp2()}
if (arg_suppnum == 7){r_supp7()}
if (arg_suppnum == 8){r_supp8()}
if (arg_suppnum == 9){r_supp9()} 
if (arg_suppnum == 1){r_nicole()}
#if (arg_suppnum != 1|2|7|8|9){warning("You need to specify the kind of input! type 2, 7, 8 or 9 for the respective supplement files in Ng et al 2018 or 1 for nicoles summary file!")}