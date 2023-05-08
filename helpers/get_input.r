
##============================================================================##
##                     METHODS FOR INPUT PARSING                              ##
##----------------------------------------------------------------------------##
##                                                                            ##
## Set of functions that help with getting input for the tool from different  ##
## file formats, specifically:                                                ##                                                                           ##
## 1) Variant Call Files (.vcf)                                               ##
## 2) Mutation Annotation Format (.maf)                                       ##
##                                                                            ##
## Also contains functions for input formatting and checking.                 ##
##============================================================================##

###############
### PARSING ###
###############

helpers.path <- "./helpers"
cc_sourced  <- FALSE

#' Function that parses a MAF file and creates a dataframe for the tool from it.
#' @param file - The according MAF file
#' @return A dataframe that can be processed by the MTB tool.
parse_maf <- function(file) {
  if (!cc_sourced) {
    source(paste(helpers.path,"/coordinate_converter.r",sep=""))
    cc_sourced = TRUE
  }
  options(stringsAsFactors = FALSE)
  
  # Parse the MAF file and create dataframe:
  full_maf <- read.delim(file, sep="\t", comment.char = "#", header=TRUE, stringsAsFactors = FALSE)
  maf_snv = data.frame(matrix(ncol=3, nrow=0), stringsAsFactors = FALSE)

  # Gather information from the main tool column
  for (i in 1:nrow(full_maf)) {
    if (full_maf[i,"Variant_Classification"] != "Silent") {
      maf_pc <- ""
      maf_gene = as.character(full_maf[i,"Hugo_Symbol"])
      maf_vc = as.character(full_maf[i,"Variant_Classification"])
      
      # Get provided protein change
      if (full_maf[i, "Protein_Change"] == "" || is.na(full_maf[i, "Protein_Change"])) {
        # Calculate protein change by positions
        query_frame <- data.frame(matrix(ncol=5, nrow=0), stringsAsFactors = FALSE)
        colnames(query_frame) = c("seq_name", "seq_start", "seq_end", "refAllele", "varAllele")
        t_allele <- full_maf[i, "Tumor_Seq_Allele1"]

        if (t_allele == full_maf[i, "Reference_Allele"]) {
          t_allele = full_maf[i, "Tumor_Seq_Allele2"]
        }

        query_frame = rbind(query_frame, data.frame(seq_name=full_maf[i, "Chromosome"], seq_start=full_maf[i, "Start_Position"], seq_end=full_maf[i, "End_Position"],
                                                    refAllele = full_maf[i, "Reference_Allele"], varAllele=t_allele))
        print(query_frame)
        cc_out <- Coordinate_Converter(query_frame, edbx_hg19)
        cc_out <- subset(cc_out, cc_out$start != -1 & cc_out$end != -1 & cc_out$names != "")
        if (is.null(cc_out) || nrow(cc_out) == 0) {
          maf_pc = "-"
        } else {
          if (cc_out[1,"aa_change"] != "") {
            print(cc_out)
            cc_split <- str_split(cc_out[1,"aa_change"], ">")[[1]]
            if (cc_out[1, "start"] == cc_out[1, "end"]) {
              maf_pc = paste(cc_split[1], cc_out[1, "start"], cc_split[2], sep="")
            } else {
              maf_pc = paste(cc_split[1], cc_out[1, "start"], "-", cc_out[1,"end"], cc_split[2], sep="")
            }
          } else { maf_pc == "" }
        }
        message(paste("Finished CCS calculation for a '", maf_vc, "' variant on gene '", maf_gene, "' for which no protein change was provided.", sep=""))
      } else {
        # Take the provided protein change
	maf_pc = as.character(gsub("p.","",full_maf[i,"Protein_Change"]))
      }
      maf_snv = rbind(maf_snv, c(maf_gene, maf_vc, maf_pc)) 
    }
  }
  colnames(maf_snv) = c("Hugo_Symbol", "Variant_Classification", "Protein_Change")
  print(maf_snv)
  return(maf_snv)
  
}

#' Function that parses a VCF file and creates a dataframe for the tool from it.
#' @param file - The according VCF file
#' @return A dataframe that can be processed by the MTB tool.
parse_vcf <- function(file) {
  if (!cc_sourced) {
    source(paste(helpers.path,"/coordinate_converter.r",sep=""))
    cc_sourced = TRUE
  }
  options(stringsAsFactors = FALSE)
  
  # Parse the VCF file and create dataframe:
  full_vcf <- read.delim(file, header=FALSE, comment.char = "#", stringsAsFactors = FALSE)[1:8] # The "format" column (9) is optional and not at all interesting for the tool
  colnames(full_vcf) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  vcf_snv = data.frame(stringsAsFactors = FALSE)
  
  for (i in 1:nrow(full_vcf)) {
    print("t1")
    temp_pc = "unknown"
    temp_vc = "unknown"
    temp_gene = "unknown"
    # Standards for this have to be checked. For now, use "HS", "VC" and "PC" for Hugo Symbol and Protein Change, respectively.
    # If Hugo Symbol, Variant Classification and Protein Change are given, just use those. Check the info column:
    additional_information <- split(full_vcf[i, "INFO"], ";")
    for (kv_pair in additional_information) {
      kvs <- split(kv_pair, "=")
      if (grepl("^HS$",kvs[1])) {
        temp_gene = kvs[2]
      } else if (grepl("^VC$", kvs[1])) {
        temp_vc = kvs[2]
      } else if (grepl("^PC$", kvs[1])) {
        temp_pc = kvs[2]
      }
    }
    
    if (temp_pc != "unknown" && temp_vc != "unknown" && temp_gene != "unknown" && temp_vc != "Silent") {
      vcf_snv = rbind(vcf_snv, c(temp_gene, temp_vc, temp_hs))
    } else {
      #This is where the genomic location mapping service has to come in
      #This is not yet supported, though!
    }
    print("t2")
    #colnames(vcf_snv) = c("Hugo_Symbol", "Variant_Classification", "Protein_Change")
    print("t2-2")
    return(vcf_snv)
  }
}

#######################
### OTHER FUNCTIONS ###
#######################

#' Function for making a string uppercase (from "Help" for 'toupper()'):
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# Function to convert 3 letter amino acid code to 1 amino acid letter code:
AAcode_3to1 <- function(aa, aacode_file){
  aa_code <- aacode_file
  if (any(grepl(aa, aa_code$Three.Letter.Code, ignore.case = TRUE))){
    aa <- gsub("[[:space:]]", "", aa)
    aa_1 <- toString(aa_code[grepl(aa, aa_code$Three.Letter.Code, ignore.case = TRUE), 3])
  } else {aa_1 <- aa}
  return(aa_1)
}

#Function to check if a cancer synonym is known or not
check_cancer <- function(cancer, c_synonyms) {
  for (i in 1:ncol(c_synonyms)) {
    for (j in 1:nrow(c_synonyms)) {
      if (c_synonyms[j,i] != "" && grepl(tolower(cancer),paste("^",tolower(c_synonyms[j,i]),"$",sep=""), fixed=TRUE)) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}
