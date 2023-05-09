#==============================================================#
#       GENERATION OF MOLECULAR TUMOR BOARD (MTB) REPORT       #
#--------------------------------------------------------------#
#   This script filters SNVs and CNVs using gene-drug public   #
#   databases. Then classifies the variants into levels of     #
#   evidence and finally presents the results in a report      #
#==============================================================#
# contact: nadine.kurz@bioinf.med.uni-goettingen.de
# contact: charlotte-hoeltermann@web.de
# contact: tim.tucholski@stud.uni-goettingen.de

#####################
## ADDITIONAL INFO ##
#####################
# Is to be run from the docker environment or, when the working directory is correctly set beforehand, may run on the command line on a local system.
# RUN THIS SCRIPT WITH: Rscript --vanilla localmaster.R [.csv file to be analyzed]
# IMPORTANT: The input files NEED to have "Fusion", "CNV" or "SNV" in them to correctly identify them.
#### REPS AND DOCS: ####
# https://github.com/jperera-bel/iMTB-Report
# https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-018-0529-2

##################
## Dependencies ##
##################
#mtbreport_packages <- c("knitr","stringr","stringi","magrittr","xtable","pander")

#mtbreport_packages <- c('DT', 'ggplot2', 'xtable', 'timeDate','timeSeries', 'pander', 'knitr', 'RCurl', 'httr', 'openssl', 'curl', 'XML', 'covr', 'stringr', 'stringi', 'plyr','dplyr', 'BiocManager', 'parallel', 'RestRserve', 'stats4', 'testthat', 'jsonlite', 'tools','devtools','gdtools','xml2','shinyjs','openxlsx','maftools')
#mtbreport_bm_packages <- c('shinydashboard','S4Vectors', 'AnnotationDbi', 'XVector', 'Biostrings', 'Biobase', 'GenomeInfoDb', 'VariantAnnotation', 'BiocGenerics', 'EnsDb.Hsapiens.v86', 'IRanges', 'GenomicRanges', 'GenomicFeatures', 'gwascat', 'biomaRt', 'liftOver', 'EnsDb.Hsapiens.v75','maftools','EnsDb.Hsapiens.v75','gwascat','liftOver','AnnotationHub','interactiveDisplayBase')

# not_installed <- mtbreport_packages[!(mtbreport_packages %in% installed.packages()[ , "Package"])]
# if(length(not_installed)) install.packages(not_installed, repos='http://cran.us.r-project.org', dependencies=TRUE)


mtbreport_bm_packages = c('VariantAnnotation')
not_installed_bm <- mtbreport_bm_packages[!(mtbreport_bm_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed_bm)) BiocManager::install(not_installed_bm)


# install_github('mariodeng/FirebrowseR')

library(devtools)

suppressPackageStartupMessages(library(knitr, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(stringr, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(xtable, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(pander, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(ggplot2, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(timeSeries, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(FirebrowseR, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(yaml, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(httr, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(tools, quietly=TRUE, warn.conflicts = FALSE))
suppressPackageStartupMessages(library(latexpdf, quietly=TRUE, warn.conflicts = FALSE))


######################################################
## Load everything necessary for downstream analysis #
######################################################

setwd('./')
helpers.path <- "./helpers"
data.path <- "./data"
getwd()

## Load helper functions:
source(paste(helpers.path,"/get_druggable.r",sep=""))
source(paste(helpers.path,"/get_levels.r",sep=""))
source(paste(helpers.path,"/get_input.r",sep=""))
source(paste(helpers.path,"/variant_annotation.r",sep=""))
# source("helpers/txDBdata_preProcess.r")
invisible(Sys.setlocale(locale = "C"))

## Load databases:
db_GDKD <- read.table(paste(data.path,"/GDKD.csv",sep=""), sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE, fill=FALSE) 
db_CIVIC = read.table(paste(data.path,"/CIViC.csv",sep=""), sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE) 
civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
db_CIVIC <- civic_del(db_CIVIC)  
colnames(db_CIVIC) <- gsub("X.","",colnames(db_CIVIC))
colnames(db_CIVIC) <- gsub("\\.","",colnames(db_CIVIC))

db_TARGETMERIC = read.table(paste(data.path,"/TARGET_MERIC.csv",sep=""), sep="\t", header=TRUE)
#db_OncoKB = read.csv2(paste(data.path,"/OncoKB.csv",sep=""), header=T, sep=";")

# Cancer type synonyms between databases:
synonyms = read.csv(paste(data.path,'/cancer_types.csv',sep=""), header = TRUE,sep="\t")

# Checked db files, they seem to be ok 

#################
### GET INPUT ###
#################

# Input is now preferably read from a single metadata file.
# The metadata file must be in YAML format and adhere to the given standards.
# It must at least contain paths to the files and the file formats (maf, csv, vcf)
# VCF and MAF files are assumed to contain SNV data. For CSV files, the variation type must be specified (snv, cnv, fusion).

args <- commandArgs(trailingOnly=TRUE)
path_keep <- data.frame()
cancer <- "unspecified"
print_yaml <- FALSE
ref_genome <- "unspecified"
var_control_sign <- ""

if (length(args)==0) {
    stop("At least one file path argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1 && grepl(".yaml|.yml", args[1], ignore.case = TRUE)) { # Case 1: Metadata file is given

    yaml_input <- read_yaml(args[1])
    print_yaml = TRUE
    
    #Check if the ref genome is provided, if not default set the ref genome as hg38
    if(is.null(yaml_input$data[["ref_genome"]])){
        ref_genome <- "hg38"
        warning(paste("The reference genome version has not been set in ",args[1],". The hg38 has been chosen as the default reference genome.", sep=""))
    } else if(!grepl("\\bhg38\\b|\\bGRCh38\\b|\\bhg19\\b|\\bGRCh37\\b", yaml_input$data[["ref_genome"]],ignore.case = TRUE)){
        ref_genome <- "hg38"
        warning(paste("The reference genome version set using ", yaml_input$data[["ref_genome"]], " is invalid, The hg38 has been chosen as the default reference genome.", sep=""))
    } else{
        ref_genome <- yaml_input$data[["ref_genome"]]
    }
    # Get the cancer type (if set). This must always be the first sub-list in the YAML file.
    #Check if the cancer matches any known cancer synonym
    if (check_cancer(yaml_input$data$cancer_type, synonyms)) {
        #Cancer is known - Use as is
        cancer = yaml_input$data$cancer_type
    }

    temp <- yaml_input$data[names(yaml_input$data) != "cancer_type"]
    temp <- temp[names(temp) != "ref_genome"]
    
    path_keep <- data.frame(matrix(unlist(temp), nrow=length(temp), byrow=T))
    path_keep <- path_keep[,c(2,1,3)]
    colnames(path_keep) = c("path", "format", "var_type")
    
    ### get variants extract control sign ###
    if(!is.null(yaml_input[["control"]][["var_control"]])){
        var_control_sign <- yaml_input[["control"]][["var_control"]]
    }

  # Alternatively, the "old" way of passing parameters can be used.
  # In the case of CSV files, there MUST be a source file and an input type. Giving a cancer type is optional.
  # In the case of MAF or VCF files, the input type is assumed to be SNV. A cancer type is optional.
  # The order of the arguments does not matter.
  
  } else { # Case 2: Old parameter style - One file plus input type and optionally the cancer type
    
      temp_inp_type <- "unknown"
      temp_filename <- "missing"
      temp_format <- "csv"
      
      for (i in 1:length(args)) {
          if (grepl("^snv$|^cnv$|^fusion$", args[i], ignore.case = TRUE)) {
              temp_inp_type = args[i]
          } else if (grepl(".csv", args[i], ignore.case = TRUE)) {
              temp_filename = args[i]
              temp_format = "csv"
          } else if (grepl(".maf", args[i], ignore.case = TRUE)) {
              temp_filename = args[i]
              temp_format = "maf"
              temp_inp_type = "snv"
          } else if (grepl(".vcf", args[i], ignore.case = TRUE)) {
              temp_filename = args[i]
              print(temp_filename)
              temp_format = "vcf"
              temp_inp_type = "snv"
          } else if (grepl(".vcf.gz", args[i], ignore.case = TRUE)) { #Began adding support for gzipped .vcf, needs work as of 04/23
              temp_filename = args[i]
              temp_format = "vcf"
              temp_inp_type = "snv"    
          } else if(grepl("\\bhg38\\b|\\bGRCh38\\b|\\bhg19\\b|\\bGRCh37\\b", args[i],ignore.case = TRUE) ){
            ref_genome <- args[i]
          } else {
            cancer = args[i]
          }
      }
      
      if(ref_genome == "unspecified"){
          ref_genome <- "hg38"
          warning(paste("The reference genome version was not set - hg38 has been chosen as the default reference genome.", sep=""))
      }
      if (temp_filename != "missing" && temp_inp_type != "unknown") {
          path_keep = rbind(path_keep, c(temp_filename, temp_format, temp_inp_type))
          colnames(path_keep) = c("path", "format", "var_type")
      } else {
          message(paste("Filename:", temp_filename))
          message(paste("Input Type:", temp_inp_type))
          stop("No (correctly formatted) file was submitted OR the input type is missing. We accept data in the following formats: [CSV, MAF, VCF] and the following input types: [SNV, CNV, Fusion]")
      }
}
#print("Check path_keep object")
#print(head(path_keep))      #path_keep checked and OK
#### generate txDB object and save it to specified path(see doceker and data/txDBdata_preProcess.r)
#GFFURI<-"ftp://ftp.ensembl.org/pub/release-100/gff3/homo_sapiens/Homo_sapiens.GRCh38.100.gff3.gz"
## generate TxDB and get the TxDB object path
#save_path <- "data"
#txDBPath <- TxDB_generate_from_URI(GFFURI,save_path)
txDBPath_hg38 <- paste(data.path,"/Homo_sapiens.GRCh38.109.gff3_txdb.sqlite",sep="") #change from 100 to 109
txDBPath_hg19 <- paste(data.path,"/Homo_sapiens.GRCh37.87.gff3_txdb.sqlite",sep="")
if(grepl("\\bhg38\\b|\\bGRCh38\\b", ref_genome,ignore.case = TRUE)){
    txDBPath <- txDBPath_hg38
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg38
} else if(grepl("\\bhg19\\b|\\bGRCh37\\b", ref_genome,ignore.case = TRUE)){
    txDBPath <- txDBPath_hg19
    Hsapiens <- BSgenome.Hsapiens.UCSC.hg19
}

###################
### COMPUTATION ###
###################

for (i in 1:nrow(path_keep)) {
    tc_out <- tryCatch({  
        inp_file = as.character(path_keep[i,1])
        inp_format = as.character(path_keep[i,2])
        inp_type = as.character(path_keep[i,3])
        
        message(paste("Processing file",as.character(i-1)))
        message(paste("Input path:", inp_file, "with input type", inp_type, "in the format",inp_format))
        
        if (grepl("csv", inp_format, ignore.case = TRUE)) {
            inputfile = read.csv(inp_file, header = TRUE, row.names=NULL, stringsAsFactors = FALSE, na.strings = c("", "NA", "-"))
        } else if (grepl("maf", inp_format, ignore.case = TRUE)) {
            #inputfile = parse_maf(inp_file)
            inputfile <- MAF_to_SNVTable(inp_file,txDBPath)       #### maf parse function from variant_annotation.r
        } else {
            #stop("VCF is not yet supported")
            #inputfile = parse_vcf(inp_file)
            inputfile <- VCF_to_SNV(inp_file,txDBPath)          #### vcf parse function from variant_annotation.r
        }

        SNV = setNames(data.frame(matrix(ncol = 3, nrow = 0), stringsAsFactors = FALSE),c("Hugo_Symbol","Variant_Classification","Protein_Change"))
        CNV = setNames(data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE),c("Gene","cn_alteration"))
        TX_file = setNames(data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE),c("Gene1","Gene2"))

        if (grepl("SNV", inp_type, ignore.case = TRUE)){
            inputfile = setNames(inputfile, c("Hugo_Symbol","Variant_Classification","Protein_Change"))
            SNV = rbind(SNV, inputfile)                                                                # This part is relevant for vcf. Basically "inputfile" is turned into "SNV"
        } else if (grepl("CNV", inp_type, ignore.case = TRUE)){
            inputfile = setNames(inputfile, c("Gene","cn_alteration"))
            CNV = rbind(CNV, inputfile)
        } else if (grepl("Fusion", inp_type, ignore.case = TRUE)){
            inputfile = setNames(inputfile, c("Gene1","Gene2"))
            TX_file = rbind(TX_file, inputfile)
        } else {
            stop("Error(s) occurred when processing the input type - Make sure it is either 'snv', 'cnv' or 'fusion'!")
        }
        if (nrow(SNV) == 0 & nrow(CNV) == 0 & nrow(TX_file) == 0){
            stop('The patient has no SNVs, CNVs or gene fusions. Not worth continuing with the analysis!')
        }

        cancer_GDKD      = unique(as.character(synonyms[grep(cancer,synonyms$tcga_cancer,ignore.case = T),"knowledge"]))
        cancer_CIVIC     = sapply(cancer,function(x) as.character(na.omit(synonyms[grep(x,synonyms$tcga_cancer,ignore.case = T),"civic"])[1]))
        cancer_CIVIC     = paste(cancer_CIVIC,collapse = ",")
        cancer_CIVIC     = unique(strsplit(cancer_CIVIC,",")[[1]])
        #cancer_ONCO      = unique(as.character(synonyms[grep(cancer,synonyms$oncokb,ignore.case = T),"oncokb"]))

        ########################################
        ### FORMATTING OF GENE FUSIONS INPUT ###
        ########################################
        
        # tx is a vector with gene rearrangements. Each gene rearrangement is defined by two "=" separated genes "gene1=gene2" ###
        # separate genes with "="
        TX = c()
        if (nrow(TX_file) > 0){
            for (r in 1:nrow(TX_file)){
                TX1 = toString(TX_file[r,1])
                TX2 = toString(TX_file[r,2])
                if (TX1 != "" & TX2 != ""){
                    v = paste(TX1, "=", TX2)
                    v = gsub(" ", "", v)
                    TX = append(TX, v, after = length(TX))
                }
            }
        }
        
        #############################################
        ### MAKE SURE INPUT IS IN REQUIRED FORMAT ###
        #############################################
        # aa_code <- read.csv("data/aa_code.csv", stringsAsFactors = FALSE)
        
        ### FORMATTING OF AA CODE IN "aa change" IN SNV ###                                       

        # Exchange 3 letter amino acid code to 1 letter amino acid code in p.XXX of SNV Input:
        if (nrow(SNV)>=1){
        ### CHECK "Variant_type" IN SNV ###
            for (num in 1:nrow(SNV)){
                if (!is.na(SNV[num, 2]) & !is.na(SNV[num, 2]) & !grepl("Frame_Shift_Del|Frame_Shift_Ins|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Nonsense_Mutation|Silent|Splice_Site|Translation_Start_Site|Nonstop_Mutation|RNA|Targeted_Region|De_novo_Start_OutOfFrame|De-novo_Start_InFrame|Start_Codon_SNP|Stop_Codon_Del", toString(SNV[num, 2]), ignore.case = TRUE)){
                    # SNV[num, 2] <- NA
                    warning(paste("This Variant type: ", toString(SNV[num, 2]), " is not accepted by MTB Report. Please refer to our documentation using -h.", toString(SNV[num, 2]), sep=""))
                }
            }
        }
        ### CHECK "CN type" IN CNV ###
        if (nrow(CNV)>=1){
            for (num in 1:nrow(CNV)){
                if (!grepl("amplification|deletion", toString(CNV[num, 2]), ignore.case = TRUE)){
                    stop(paste("This CN type: ", toString(CNV[num, 2]), " is not accepted by MTB Report. It must be \"deletion\" or \"amplification\". Please refer to our documentation using -h.", sep=""))
                }
            }
        }
        
        #####################################
        ## Filter SNVs and CNVs by database #  
        #####################################       
        print("locwtcivic0")                         

        #### GDKD DB
        druggableGDKD = data.frame()
        if (nrow(SNV) >=1) {
            druggableGDKD = match_SNV_GDKD(SNV,db_GDKD)
        }
        if (nrow(CNV) >=1) {
            druggableGDKD = rbind(druggableGDKD,match_CNV_GDKD(CNV,db_GDKD))
        }
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
            druggableGDKD = rbind(druggableGDKD,match_WT_GDKD(SNV,CNV,cancer_GDKD,db_GDKD))
        }
        if (length(TX) >= 1) {
            druggableGDKD = rbind(druggableGDKD,match_TX_GDKD(TX,db_GDKD))
        }
        rownames(druggableGDKD) = NULL

        print("locwtcivic1")
        #### CIVIC
        druggableCIVIC = data.frame()
        if (nrow(SNV) >=1) {
            druggableCIVIC = match_SNV_CIVIC(SNV, db_CIVIC)
        }
        print("locwtcivic1-1")
        if (nrow(CNV) >= 1) {
            druggableCIVIC = unique(rbind(druggableCIVIC, match_CNV_CIVIC(CNV, db_CIVIC)))
        }
        print("locwtcivic1-2")
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
            druggableCIVIC = unique(rbind(druggableCIVIC, match_WT_CIVIC(SNV, CNV, cancer_CIVIC, db_CIVIC)))
        }
        if (length(TX) >= 1) {
            druggableCIVIC = rbind(druggableCIVIC,match_TX_CIVIC(TX,db_CIVIC))
        }
        rownames(druggableCIVIC) = NULL
        print("locwtcivic1-3")

        #### TARGET DB
        druggableTARGET = data.frame()
        if (nrow(SNV) >= 1) {
            druggableTARGET = match_TARGET_MERIC(SNV, db_TARGETMERIC) #added 27.4.23
        } 
        if (nrow(SNV) >= 1 || nrow(CNV) >= 1 || length(TX) >= 1) {
            druggableTARGET = match_TARGET_MERIC(SNV, CNV, TX, db_TARGETMERIC)
        }       
        print("wtcivic 2")

        #### OncoKB
        #druggableONCO = data.frame()
        #if (nrow(SNV) >=1) {
        #    druggableONCO = match_SNV_OncoKB(SNV,db_OncoKB)
        #}
        #if (nrow(CNV) >=1) {
        #    druggableONCO = rbind(druggableONCO,match_CNV_OncoKB(CNV,db_OncoKB))
        #}
        #if (length(TX) >= 1) {
        #    druggableONCO = rbind(druggableONCO,match_TX_OncoKB(TX,db_OncoKB))
        #}
        #rownames(druggableONCO) = NULL

        #if( nrow(druggableGDKD) == 0 && nrow(druggableCIVIC) == 0 && nrow(druggableTARGET)==0 && nrow(druggableONCO) == 0) {
        #    warning('No gene-drug interations were found. Not worth continuing with the analysis!')
        #}


        #########################################################
        ## Classify filtered variants by Levels of Evidence    ##
        #########################################################
        sort = "levels"
        print("levels 0: Classify filtered variants by Levels of Evidence")

        ### KNOwLEDGE
        levelsGDKD = c()
        if (nrow(druggableGDKD) > 0) {
            levelsGDKD = get_levels_GDKD(druggableGDKD, cancer_GDKD)
        }

        ### CIVIC
        levelsCIVIC = c()
        if (nrow(druggableCIVIC) > 0) {
            levelsCIVIC = get_levels_CIVIC(druggableCIVIC, cancer_CIVIC)
        }
        ### ONCOKB
        #levelsONCO = c()
        #if (nrow(druggableONCO > 0)) {
        #    levelsONCO = get_levels_ONCO(druggableONCO, cancer_ONCO)
        #}

        # levelsCIVIC = return_list
        # levelsGDKD = c()

        if (length(levelsCIVIC) > 0 && length(levelsGDKD) > 0) {
            levels = merge_levels(levelsGDKD, levelsCIVIC)
        } else if (length(levelsCIVIC) > 0) {
            levels = levelsCIVIC
        } else if (length(levelsGDKD) > 0) {
            levels = levelsGDKD
        } else {
            levels = c()
        }

        # Homogenize/clean the final table
        table = clean_levels(levels, synonyms=synonyms, sort_by=sort) #removed "", levelsONCO ""     <-----  removed O argument from function in get_levels.r

        #print("Check cleaned levels")                                                                                        #levels is empty even after gene symbol mod, check code above
        #print(head(table))

        # table = A
        #################################################
        ## Generate report from table and patient data ##
        #################################################
        print("inp0") 

        # Create a name for the output file:
        inp_nice = strsplit(inp_file, "/", fixed = TRUE)    #when testing using mutect2 output: "inp_name could not be found"
        inp_name = tail(inp_nice[[1]], n=1)
        inp_name = gsub(" ", "", inp_name)
        # Parse the input file location:
        inp_loc = inp_nice[[1]]
        inp_loc = head(inp_loc, -1)
        loc = ""
        for (i in 1:length(inp_loc)){
            loc = paste(loc, inp_loc[i], "/", sep = "")
        }
        print("LOC:")
        print(loc)
        
        # Fix formatting issues:
        # Reformat the table to Gene1-Gene2 for better legibility
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub("NEW", "", x))
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub("([A-Z]) ([A-Z])", "\\1, \\2", x))
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub(" ", "", x))
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub("=", "-", x))
        table$`Sample`<- rep(inp_name, nrow(table))
        table = table[, c(10,1,2,3,4,5,6,7,8,9)]
        table$Drugs <- sapply(table$Drugs, function(x) capwords(toString(x)))
        table$Drugs <- sapply(table$Drugs, function(x) paste(sort(unlist(strsplit(x, ", |,"))), collapse=", "))

        #Reformat variants containing "*" for better legibility
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub("([A-Z])\\*", "\\1,\\*", x))
        table$Pat_Var <- sapply(table$Pat_Var , function(x) gsub("\\*([A-Z])", "\\*,\\1", x))

        ### Extract table by variants type(any Mut, sepcific or wildtype) ###
        Specific_MutTB_extract<-function(table,var_control_sign){
          if(var_control_sign == "all" || var_control_sign == ""){
              return(table)
          } else if(var_control_sign == "without_wildtype"){
              temp <- table[- grep("wild", table$Known_Var,ignore.case = TRUE), ]
              return(temp)
          } else if(var_control_sign == "wildtype_only"){
              temp <- table[grep("wild", table$Known_Var,ignore.case = TRUE), ]
              return(temp)
          } else if(var_control_sign == "specific_var"){
              temp <- table[- grep("wild|any", table$Known_Var,ignore.case = TRUE), ]
              return(temp)
          } else if(var_control_sign == "any_var_only"){
              temp <- table[grep("any", table$Known_Var,ignore.case = TRUE), ]
              return(temp)
          } else {
              warning(paste("The control sign set in yaml metadata file with ", yaml_input[["control"]][["var_control"]], " is invalid, all the known variants will be extracted", sep=""))
              return(table)
          }
        }
        table <- Specific_MutTB_extract(table,var_control_sign)
        table <- table[!duplicated(table), ]
        table$Predicts <- gsub("sensitivity/response", "sensitivity/ response", table$Predicts)
        
        ### SAVE AS CSV: ###
        # Create a directory:
        print("create dir")
        dir.create(paste(loc, "MTB_Reports", sep=""), showWarnings = FALSE)
        print("dir create ok")

        # Do not write any empty reports, but report those inputs that were empty:
        if (nrow(table) == 0){
            warning(paste("There were no results for this input: ", inp_name, " No report file was written."))
        } else {
            inp_name = gsub(".csv|.maf|.vcf", "", inp_name)
            write.csv(table, paste(loc, "MTB_Reports/", inp_name, "_report.csv", sep = ""), row.names = F)
            # message for command line:
            print(paste("MTB_Reports/", inp_name, "_report.csv", " has successfully been written. ", sep = ""))
        }

    },
    error = function(cond) {
        message(paste("Error encountered while processing file:", inp_file, sep=" "))
        message("Here's the original error message:")
        message(paste(cond,"\n",sep=""))
        return(NA)
    })
  # End TryCatch
}
# End For

#####################
###  FINAL STEPS  ###
#####################

# Append some tool information to the passed metadata file and save it to output as well
if (!print_yaml) { yaml_input=yaml.load(paste("data:\n path: ", inp_name, "\n cancer_type: ", cancer, "\n variation_type: ", inp_type, "\n ref_genome: ", ref_genome, "\n", sep="")) }
metafile_addendum <- yaml.load(paste("tool:\n git_checksum: undefined\n date: ", Sys.Date(), "\ndatabases:\n gdkd: v20.0\n civic: 15-jan-01\n target-meric: v3"))#\n oncokb: v1"))
invisible(write_yaml(append(yaml_input, metafile_addendum), paste(loc, "MTB_Reports/", inp_name, "_metadata.yml", sep="")))

# Copy file with database information to output:
invisible(file.copy(paste(data.path,"/MTB_Databases_versions.csv",sep=""), "loc", overwrite = FALSE, recursive = FALSE))

# Show warnings:
# summary(warnings())

#####################
# OUTPUT PDF REPORT #
#####################

print("Create .pdf report")

setwd('./')
wd <- getwd()
f = paste(wd,"/src/tmp/",sep="")
report_file <- paste(loc, "MTB_Reports/",inp_name,".tex",sep="")
report = file.path(wd, "/src/Rnw/Report_local-knitr.Rnw")
knit2pdf(input = report, output = report_file)
print("Knitting done")
