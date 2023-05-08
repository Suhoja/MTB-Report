##============================================================================##
##                     MATCHING GENE-DRUG INTERACTIONS                        ##
##----------------------------------------------------------------------------##
##                                                                            ##
## Set of functions that search for matching variants (SNVs, CNVs and fusion  ##
## genes) in specialized  gene-drug interaction databases, specifically:      ##
##                                                                            ##
## 1) GDKD: Dienstmann et al., Cancer discovery 5.2 (2015): 118-123, v19      ##
## 2) CIVic: Griffith et al., bioRxiv (2016): 072892                          ##
## 3) TARGET: Van Allen et al., Nature medicine 20.6 (2014): 682-688, v3      ##
## 4) Meric-Bernstam et al., J Natl Cancer Inst. 107(7) (2015)                ##
##                                                                            ##
##============================================================================##

library(stringr)
library(plyr)

#############
## 1) GDKD ##
#############
data.path <- './data'


# Debugging Switch
debug = FALSE

################
## 2) GDKD    ##
################

match_SNV_GDKD = function(snv,db){
    ## Given a list of SNVs, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
    ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
    ## output: returns those Gene Drug Knowledge Database rows matching to the input SNVs.
    if (nrow(snv) == 0) {
        stop('Provided SNV file was empty!')
    }
    else if (ncol(snv) < 3) {
        stop('The input SNV file should have at least three columns!')
    }
    if (missing(db)){
        gdkd <- read.table(paste(data.path,"/GDKD.csv",sep=""), sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,fill=FALSE)
    } else {
        gdkd <- db
    }
    
    gdkd$Patient_variant = ""
    druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
    colnames(druggable)       = colnames(gdkd)
    colnames(druggable)       = colnames(gdkd)

    ## We apply a narrowing proceedure
    ## Starting with all variants of the database
    ## We filter out genes in the database not mutated in the patient
    gdkd$Gene = gsub(" ", "",gdkd$Gene)
    index               = as.character(gdkd$Gene) %in% as.character(snv[,1])  #Database genes mutated in patient
    druggable           = subset(gdkd,index)                                  #Subset from complete DB
    rownames(druggable) = NULL

    #print("Check druggable GDKD after first narrowing") #Empty. Looks like gdkd has gene names, whereas SNV has ensembl_id -- fix'd in variant_annotation.r
    #print(head(druggable))

    ## From the subset previously filtered by gene
    ## We annotate only matching variants or repurposed ones (NEW)
    if (nrow(druggable) != 0){
        for (i in 1:nrow(druggable)){
            gene       =  druggable[i,"Gene"]
            p_var_type = as.character(snv[which(as.character(snv[,1])==gene),2]) # Patient's variant type on current gene
            p_variants = unique(as.character(snv[which(as.character(snv[,1])==gene),3])) # Patient's protein change on the current gene
            for (p_variant in p_variants){
                ### Match Missense mutations with variant not specified in the database
                if (grepl("missense|mutation|any",druggable[i,"Description"], ignore.case = TRUE)
                  && grepl("mut|any|unknown",druggable[i,"Variant"], ignore.case = TRUE)
                  && grepl("missense|frame_shift",p_var_type,ignore.case = TRUE)){
                    druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                }
                if (grepl("loss|lof|fs",druggable[i,"Variant"], ignore.case = TRUE)
                  && grepl("splice|shift",p_var_type,ignore.case = TRUE)){
                    druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                }
                ### Match Missense mutations with specific variant specified in the database
                if (grepl("missense",druggable[i,"Description"], ignore.case = TRUE) 
                    && grepl("missense",p_var_type,ignore.case = TRUE) 
                    && grepl("mut|any|unknown",druggable[i,"Variant"],ignore.case = TRUE)==FALSE){
                    if (grepl(p_variant,druggable[i,"Variant"])){
                        druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant) 
                    }
                    ## Repurposing - New variants
                    else {
                        # if (grepl(gene,"BRAF,KIT,FGFR3,KRAS,PDGFRA,EGFR,AKT1,MTOR,ALK",ignore.case = TRUE)){
                        #     if (gene=="EGFR" &&  druggable[i,"Variant"] != "T790M" && druggable[i,"Variant"] != "C797S"){
                        #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
                        #     if (gene=="MTOR" && druggable[i,"Variant"] != "F2108L"){
                        #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
                        #     if (gene=="AKT1" && druggable[i,"Variant"] != "Q79K"){
                        #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
                        #     if (gene=="BRAF" && grepl("V600",druggable[i,"Variant"])==F){
                        #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")}
                        # }
                        # else {
                        #     druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                        # }
                        if (grepl('^(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)*(?:[a-zA-Z]*\\d+[a-zA-Z]*)$', druggable[i,"Variant"], ignore.case = TRUE)){
                              druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                        }
                    }
                }
                if (grepl("exon|delins|indel",druggable[i,"Variant"], ignore.case = TRUE)
                    && grepl("missense|frame_shift",p_var_type,ignore.case = TRUE)){
                      druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                }
                ### Match indels, small deletions and insertions - repurposing rule (position is not checked)
                if (grepl("indel|deletion mutation|insertion mutation",druggable[i,"Description"], ignore.case = TRUE) &&
                    grepl("ins|del",p_var_type,ignore.case = TRUE) ){
                    druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                }

                ### Match Nonsense mut - repurposing rule (deletion, LoF)
                if (grepl("loss-of-function",druggable[i,"Effect"], ignore.case = TRUE) &&
                    grepl("nonsense",p_var_type,ignore.case = TRUE) ){
                    druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                }

                ### Find matches in variants annotated at exon level (exon X p.XX-XX)
                if (grepl("-",druggable[i,"Variant"])){
                    aa_pos = as.numeric(str_match(p_variant,"[0-9]+"))
                    range  = as.character(str_match_all(as.character(gdkd[i,"Variant"]),"[0-9]+-[0-9]+")[[1]])
                    range  = strsplit(range,"-")
                    if (length(range) != 0){
                        for (x in 1:length(range)){
                            if (aa_pos >= range[[x]][1] && aa_pos <= range[[x]][2]){
                                druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant) 
                            }
                        }
                    }
                }
            }
        }
    } # This is indentation hell!

    index= which(druggable$Patient_variant != "")     ## Keep all those entries with "Patient_variant"
    druggable=druggable[index,]   #Subset from first selection
    return(druggable)
  #print("Check druggable GDKD after annotation")
  #print(head(druggable))
}

match_CNV_GDKD = function(cnv,db){
  ## Given a list of CNVs, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## output: returns those Gene Drug Knowledge Database rows matching to the input CNVs.
  if (nrow(cnv) == 0) {
    stop('Provided CNV file was empty!')
  }
  else if (ncol(cnv) < 2) {
    stop('The input CNV file should have atleast two columns!')
  }

  if (missing(db)){
    gdkd=read.delim(paste(data.path,"/GDKD.csv",sep=""),sep="\t")
  } else {
    gdkd=db
  }
  gdkd$Patient_variant = ""
  druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)

  ## We apply a narrowing proceedure
  ## Starting with all variants of the database
  ## We filter out genes in the database with no CNVs in the patient
  index               =gdkd$Gene %in% as.character(cnv[,1]) #Database genes mutated in patient
  druggable           =subset(gdkd,index)
  rownames(druggable) = NULL

  # Retain "copy number gain/loss" variants in the database
  index               =grep("copy number|promoter methylation or deletion",druggable[,"Description"])
  druggable           =druggable[sort(unique(index)),]
  rownames(druggable) = NULL

  # Matching patient variants
  index=character()
  if (nrow(druggable) != 0){
    for (i in 1:nrow(druggable)){
      gene = as.character(druggable[i,"Gene"])
      pvar = as.character(cnv[which(cnv[,1]==gene),2])[1]

      #  "mutation or copy number loss"
      if (grepl("copy number loss|promoter methylation or deletion",druggable[i,"Description"]) && pvar == "deletion"){
       druggable[i,"Patient_variant"] = pvar
      }
      if (grepl("copy number gain",druggable[i,"Description"]) && pvar == "amplification"){
        druggable[i,"Patient_variant"] = pvar
      }
    }
  }
  index               = which(druggable$Patient_variant != "")
  druggable           = druggable[index,]
  rownames(druggable) = NULL
  return(druggable)
}

match_TX_GDKD = function(tx,db){
  ## Given a list of gene rearrangements, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two "="-separated genes "gene1=gene2"
  # tx = fusion
  if (length(tx) == 0) {
    stop('Provided fusion file was empty!')
  }
  for (element in tx) {
    if (!grepl('=', element)) {
      warning('Fusion file might contain wrongly formatted lines.')
    }
  }

  if (missing(db)){
    gdkd=GDKD <- read.table(paste(data.path,"/GDKD.csv",sep=""), sep="\t", header=TRUE, comment.char="#",
                            na.strings=".", stringsAsFactors=FALSE,
                            fill=FALSE)
  } else {
    gdkd=db
  }
  gdkd$Patient_variant = ""
  druggable                 = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  colnames(druggable)       = colnames(gdkd)

  # Retain "fusion gene" or "rearrangement" variants
  index               = grep("fusion gene|rearrangement",gdkd[,"Description"])
  druggable           = gdkd[index,]

  druggable$Patient_variant = "" #added

  rownames(druggable) = NULL

  # grep patient single gene to gdkd$Gene
  for (fusion in tx){
    genes  = strsplit(as.character(fusion),"=")[[1]]
    for (gene in genes){
      i  = grep(gene,druggable$Gene,ignore.case=TRUE)
      druggable[i,"Patient_variant"] = fusion
    }
  }

  ## Retain rows with "Patient_variant"
  index = character()
  index = which(druggable$Patient_variant != "")
  druggable = druggable[index,]
  rownames(druggable) = NULL

  return(druggable)
}

match_WT_GDKD = function(snv, cnv, cancer_gdkd,db){
  ## Searches for wild type variants in Gene Drug Knowledge Database not matching any gene in the inputs SNVs or CNVs.
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: cancer_gdkd is a character vector with a cancer type used by GDKD.
  ## output: returns those GDKD wild type rows (annotated for disease=cancer_gdkd) not matching to the input SNVs and CNVs.
  if (nrow(snv) == 0 && nrow(cnv) == 0) {
    stop('Both of the provided variant files were empty!')
  }
  else if (ncol(snv) < 3 && ncol(cnv) < 2) {
    stop('At least one of the input variant files is wrongly formatted!')
  }

  if (missing(db)){
    gdkd <- read.table(paste(data.path,"/GDKD.csv",sep=""), sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE, fill=FALSE)
  } else {
    gdkd <- db
  }
  print("druggable matchwt 0")
  gdkd$Patient_variant <- ""
  druggable = data.frame(matrix(ncol=ncol(gdkd),nrow=0))
  print("druggable match 0-0-1")
  colnames(druggable) <- colnames(gdkd)
  print("druggable match 0-1")

  wt <- grep("wild type",gdkd[,"Variant"],ignore.case=TRUE)
  if (!is.null(cnv) && nrow(cnv) >= 1) {
    wt  = wt[!(gdkd[wt,"Gene"] %in% as.character(cnv[,1]))] # keep those that are not in CNV
  }
  if (!is.null(snv) && nrow(snv) >= 1) {
    wt  = wt[!(gdkd[wt,"Gene"] %in% as.character(snv[,1]))] # keep those that are not in SNV
  }
  print("druggable match1")
  #check cancer type=cancer_gdkd
  same_cancer=c()
  for (ck in cancer_gdkd){
    same_cancer = c(which(gdkd[wt,"Disease"]==ck),same_cancer)
  }
  druggable  = gdkd[wt[unique(same_cancer)],]
  if (nrow(druggable) > 0) {
    druggable$Patient_variant="WILDTYPE"
  }
  return(druggable)
}


match_SNV_CIVIC = function(snv,db){
    ## Given a list of SNVs searches for matching variants in CIViC database.
    ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
    ## output: returns those CIViC rows matching to the input SNVs.
    if (nrow(snv) == 0) {
        stop('Provided SNV file was empty!')
    }
    else if (ncol(snv) < 3) {
        stop('The input SNV file should have atleast three columns!')
    }

    if (missing(db)){
        civic <- read.csv(paste(data.path,"/CIViC.csv",sep=""), sep="\t", header=TRUE, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=TRUE)
        civic = as.data.frame(sapply(civic, function(x) gsub("\"", "", x)))
        # civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
        # civic <- civic_del(civic)
        
        colnames(civic) <- gsub("X.","",colnames(civic))
        colnames(civic) <- gsub("\\.","",colnames(civic))
    } else {
        civic <- db
    }

    civic$Patient_variant = ""
    druggable             = data.frame(matrix(ncol=ncol(civic),nrow=0))
    colnames(druggable)   = colnames(civic)

    # Gene match
    index               = as.character(civic$gene) %in% as.character(snv[,1])  #Database genes mutated in patient
    druggable           = subset(civic,index)                                  #Subset from complete DB
    rownames(druggable) = NULL

  # Keep matching variants
    if (nrow(druggable) != 0){
        for (i in 1:nrow(druggable)){
            gene       = druggable[i,"gene"]
            p_var_type = as.character(snv[which(as.character(snv[,1])==gene),2])
            p_variants = unique(as.character(snv[which(as.character(snv[,1])==gene),3]))
            for (p_variant in p_variants){
                ### Any mutation
                #if (grepl("MUTATION|LOSS OF FUNCTION|UNDEREXPRESSION|OVEREXPRESSION|[a-zA-Z]\\d+[a-zA-Z]", druggable[i,"variant"], ignore.case = TRUE)){
                    if (grepl("MUTATION",druggable[i,"variant"], ignore.case = TRUE)){
                        druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"], p_variant)
                    }
                    if (grepl("loss|lof|deletion|fs",druggable[i,"variant"], ignore.case = TRUE)
                      && grepl("splice|shift",p_var_type,ignore.case = TRUE)){
                        druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                    }
                    ### frameshift,ins,del
                    if (grepl("FRAMESHIFT|FRAME SHIFT|FS",druggable[i,"variant"], ignore.case = TRUE) 
                      && grepl("frame_shift", p_var_type, ignore.case = TRUE)){
                        druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"], p_variant)
                    }
                    ### stopgain ~ deletion
                    if (grepl("DELETION|LOSS",druggable[i,"variant"], ignore.case = TRUE) 
                      && grepl("nonsense", p_var_type, ignore.case = TRUE)){
                        druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW", p_variant,sep=" "))
                    }
                    ### Missense & Specific variant
                    if (grepl("MISSENSE", p_var_type[match(p_variant,p_variants)], ignore.case = TRUE)){
                        if (grepl(p_variant, druggable[i,"variant"])){
                            druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"], p_variant)
                        }
                        # Repurposing - New variants
                        else {
                            # if (grepl(gene,"BRAF,KIT,FGFR3,KRAS,PDGFRA,EGFR,AKT1,MTOR,ALK",ignore.case = TRUE)){
                            #     if (gene=="EGFR" &&  druggable[i,"variant"] != "T790M" && druggable[i,"variant"] != "C797S"){
                            #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                            #     }
                            #     if (gene=="MTOR" && druggable[i,"variant"] != "F2108L"){
                            #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                            #     }
                            #     if (gene=="AKT1" && druggable[i,"variant"] != "Q79K"){
                            #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                            #     }
                            #     if (gene=="BRAF" && grepl("V600",druggable[i,"variant"])==F){
                            #         druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                            #     }
                            # }
                            # else {
                            #     druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],paste("NEW",p_variant,sep=""),sep=" ")
                            # }

                            if (grepl('^(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)*(?:[a-zA-Z]*\\d+[a-zA-Z]*)$|delins|indel|exon', druggable[i,"variant"], ignore.case = TRUE)){
                              druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
                            }
                        }
                    }

                #}
            }
        }
    }

    ## Keep all those entries with "Patient_variant" and WT
    index      = which(druggable$Patient_variant != "")
    druggable  = druggable[index,]   #Subset from first selection

    return(druggable)
  #print("Check druggable CIVIC after annotation")
  #print(head(druggable))  
}

match_CNV_CIVIC = function(cnv,db){
  ## Given a list of CNVs searches for matching variants in CIViC database.
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## output: returns those CIViC rows matching to the input CNVs.
  if (nrow(cnv) == 0) {
    stop('Provided CNV file was empty!')
  }
  else if (ncol(cnv) < 2) {
    stop('The input CNV file should have atleast two columns!')
  }

  if (missing(db)){
    civic <- read.table(paste(data.path,"/CIViC.csv",sep=""), sep="\t", header=TRUE, comment.char="#", na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)
    civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
    civic <- civic_del(civic)
    
    colnames(civic) <- gsub("X.","",colnames(civic))
    colnames(civic) <- gsub("\\.","",colnames(civic))
  } else {
    civic <- db
  }
  civic$Patient_variant = ""
  druggable             = data.frame(matrix(ncol=ncol(civic),nrow=0))
  colnames(druggable)   = colnames(civic)

  # Gene match
  index     = as.character(civic$gene) %in% as.character(cnv[,1])  #Database genes mutated in patient
  druggable = subset(civic,index)                                  #Subset from complete DB
  rownames(druggable) = NULL

 # WT
 # wt  = grep("wild type",civic[,"variant"],ignore.case=TRUE)
 # wt  = wt[!(civic[wt,"gene"] %in% as.character(cnv[,1]))] # keep those that are not in cnv

  # Keep "amplification" or "deletion" variants
  index     = grepl("amplification|^deletion$|loss",druggable[,"variant"],ignore.case=TRUE)
  druggable = druggable[sort(unique(index)),]
  rownames(druggable) = NULL

  # Keep matching patient variants
  index=character()
  if (nrow(druggable) != 0){
    for (i in 1:nrow(druggable)){
      gene = as.character(druggable[i,"gene"])
      pvar = as.character(cnv[which(cnv[,1]==gene),2])[1]
      # Get "mutation or copy number loss"
      if (grepl("loss|deletion",druggable[i,"variant"],ignore.case=TRUE) && pvar == "deletion"){
        druggable[i,"Patient_variant"] = pvar
      }
      if (grepl("amplification",druggable[i,"variant"],ignore.case=TRUE) && pvar == "amplification"){
        druggable[i,"Patient_variant"] = pvar
      }
    }
  }
  index               = which(druggable$Patient_variant != "")
  druggable           = druggable[index,]
 # druggable           = rbind(druggable,civic[wt,])
  rownames(druggable) = NULL

  return(druggable)
}

match_TX_CIVIC = function(tx,db){
  ## Given a list of gene rearrangements, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
  ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two "="-separated genes "gene1=gene2"
  if (length(tx) == 0) {
    stop('Provided fusion file was empty!')
  }
  for (element in tx) {
    if (!grepl('=', element)) {
      warning('Fusion file might contain wrongly formatted lines.')
    }
  }

  if (missing(db)){
    print("LOAD CIVIC")
    civic <- read.csv(paste(data.path,"/CIViC.csv",sep=""), sep="\t", header=TRUE, comment.char="#", na.strings=".", stringsAsFactors=FALSE, quote="", fill=TRUE)
    civiv = colwise(function(x) str_replace_all(x, '\"', ""))
    #civic = as.data.frame(sapply(civic, function(x) gsub("\"", "", x)))
  } else {
    print("CIVIC ALREADY LOADED")
    civic=db
    #print(colnames(civic))
  }

  druggable                 = data.frame(matrix(ncol=ncol(civic),nrow=0))
  colnames(druggable)       = colnames(civic)

  # Retain "fusion gene" or "rearrangement" variants
  index               = grepl("fusion|::", civic[,"variant"], ignore.case=TRUE)
  druggable           = civic[index,]

  druggable$Patient_variant = "" #added

  rownames(druggable) = NULL

  print("This is really loaded")

  # grep patient single gene to civic$Gene
  for (fusion in tx){
    genes  = strsplit(as.character(fusion),"=")[[1]]
    for (gene in genes){
      i  = grepl(gene, druggable$"gene", ignore.case=TRUE)
      druggable[i,"Patient_variant"] = fusion
    }
  }

  ## Retain rows with "Patient_variant"
  index = character()
  index = which(druggable$Patient_variant != "")
  druggable = druggable[index,]
  rownames(druggable) = NULL

  return(druggable)
}


match_WT_CIVIC  = function(snv,cnv,cancer_civic,db){
  ## Searches for wild type variants in CIViC not matching any gene in the inputs SNVs or CNVs.
  ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: cancer_civic is a character vector with a cancer type used by CIViC.
  ## output: returns those CIViC wild type rows (annotated for disease=cancer_civic) not matching to the input SNVs and CNVs.
  if (nrow(snv) == 0 && nrow(cnv) == 0) {
    stop('Both of the provided variant files were empty!')
  }
  else if (ncol(snv) < 3 && ncol(cnv) < 2) {
    stop('At least one of the input variant files is wrongly formatted!')
  }

  if (missing(db)){
    civic <- read.table(paste(data.path,"/CIViC.csv",sep=""), sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)
    civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
    civic <- civic_del(civic)
    
    colnames(civic) <- gsub("X.","",colnames(civic))
    colnames(civic) <- gsub("\\.","",colnames(civic))
  } else {
    civic <- db
  }
  print("match civic wt0")
  wt  = grep("wild type",civic[,"variant"],ignore.case=TRUE)
  if (!is.null(cnv) && nrow(cnv) >= 1) {
    wt  = wt[!(civic[wt,"gene"] %in% as.character(cnv[,1]))] # keep those that are not in CNV
  }
  if (!is.null(snv) && nrow(snv) >= 1) {
    wt  = wt[!(civic[wt,"gene"] %in% as.character(snv[,1]))] # keep those that are not in SNV
  }

  #check cancer type=cancer_civic
  same_cancer=c()
  for (ck in cancer_civic){
    same_cancer = c(which(civic[wt,"disease"]==ck),same_cancer)
  }

  druggableCIVIC_wt  = civic[same_cancer,]
  if (nrow(druggableCIVIC_wt) > 0) {
    druggableCIVIC_wt$Patient_variant="WILDTYPE"
  }
  print("match civic wt1")
  return(druggableCIVIC_wt)
}

############################
## 3) TARGET              ##
## 4) MERIC-BERNSTAM DB   ##
############################

match_TARGET_MERIC = function(snv=c(),cnv=c(),tx=c(),db){
  ## Given a list of SNVs, CNVs and gene rearrangements, searches for matching variants in TARGET database and gene list from Meric-Bernstam et al 2015.
   ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
  ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
  ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two "="-separated genes "gene1=gene2"
  ## output: returns those TARGET + MERIC-BERNSTAM DB rows matching to the input SNVs, CNVs and gene rearrangements.
  if (nrow(snv) == 0 && nrow(cnv) == 0 && length(tx) == 0) {
    stop('All of the provided variant files were empty!')
  }

  if (missing(db)){
    target_meric <- read.table(paste(data.path,"/TARGET_MERIC.csv",sep=""), sep="\t", header=TRUE)
  } else {
    target_meric <- db
  }
  #target_meric = read.delim("../data/TARGET_MERIC.csv",sep="\t") # localmaster.r loads databases now.
  druggable  = data.frame(matrix(ncol=5,nrow=0))
  druggable2 = data.frame(matrix(ncol=5,nrow=0))
  druggable3 = data.frame(matrix(ncol=5,nrow=0))
  druggable4 = data.frame(matrix(ncol=5,nrow=0))

  colnames(druggable)  = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable2) = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable3) = c(colnames(target_meric)[1:4],"Patient_variant")
  colnames(druggable4) = c(colnames(target_meric)[1:4],"Patient_variant")

  print("targetmeric0: data frames created")

  # Match SNVs
  if(!is.null(snv) && !nrow(snv)==0) {
    index     = target_meric$Gene %in% as.character(snv[grep("missense",snv[,2],ignore.case=TRUE),1])
    druggable = subset(target_meric,index,1:4)
    druggable$Patient_variant = snv[na.omit(match(target_meric$Gene,as.character(snv[grep("missense",snv[,2],ignore.case=TRUE),1]))),3]
    index     = grepl("Mutation",druggable[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable = subset(druggable,index)
    print("targetmeric 1: SNVs matched")
    #print("Check targetmeric 1")
    #print(head(druggable))
  }


  # Match CNVs
  if(!is.null(cnv) && !nrow(cnv)==0) {
    index      = target_meric$Gene %in% as.character(cnv[,1])
    druggable2 = subset(target_meric,index,1:4)
    druggable2$Patient_variant = cnv[na.omit(match(target_meric$Gene,as.character(cnv[,1]))),2]
    index      = grepl("amplification|deletion|biallelic",druggable2[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable2 = subset(druggable2,index)
    # Match specific structural alteration
    if (nrow(druggable2)!=0){
      druggable2[,"Types_of_recurrent_alterations"] = gsub("[Ii]nactivation","deletion",druggable2[,"Types_of_recurrent_alterations"])
      index      = paste(lapply(1:nrow(druggable2), function(x) grep(druggable2[x,5],druggable2[x,3], ignore.case = TRUE)))
      index      = index==1
      druggable2 = druggable2[index,]
    }
    print("targetmeric 2: CNVs matched")
  }
  
  # Match TX
  if(!is.null(tx) && !length(tx)==0) {
    druggable3  = subset(target_meric,grepl("Rearrangement",target_meric[,"Types_of_recurrent_alterations"]),1:4)
    for (fusion in tx){
      genes  = strsplit(fusion,"=")[[1]]
      for (gene in genes){
        i  = grep(gene,druggable3$Gene,ignore.case=TRUE)
        druggable3$Patient_variant = c(NA)
        druggable3[i,"Patient_variant"] = paste(fusion,"fusion",sep=" ")
      }
    }
    index                = character()
    index                = which(druggable3$Patient_variant != "")
    druggable3           = druggable3[index,]  #Subset those with Patient variant
    rownames(druggable3) = NULL
    print("targetmeric 3: TX matched")
  }
  

  # Match Nonsense mut - repurposing rule (nonsense ~ deletion, biallelic inactivation)
  if(!is.null(snv) && nrow(snv)!=0){
    index      = target_meric$Gene %in% as.character(snv[grep("nonsense",snv[,2],ignore.case=TRUE),1])
    druggable4 = subset(target_meric,index,1:4)
    druggable4$Patient_variant = snv[na.omit(match(target_meric$Gene,as.character(snv[grep("nonsense",snv[,2],ignore.case=TRUE),1]))),3]
    index      = grepl("deletion|inactivation",druggable4[,"Types_of_recurrent_alterations"],ignore.case = TRUE)
    druggable4 = subset(druggable4,index)
  }

  # rbind snv + cnv + tx druggable
  #if(nrow(druggable)==0 || nrow(druggable2)==0 || nrow(druggable3)==0 || nrow(druggable4)==0){
  #  warning("No data for one of the following: snv, cnv, fusions") #  This warning is always true even if all input is provided - Seems to be broken?
  #}
  print("targetmeric 4: Nonsense mutations matched")
  druggable = rbind(druggable,druggable2, druggable3, druggable4)
  druggable = druggable[,c(1,5,2,3,4)]
  if (nrow(druggable) >= 1) {
  	druggable["database_origin"] = "T/M"
  }
  return(druggable)
  #print("Check targetmeric druggable")
  #print(head(druggable))
}

################
## 5) OncoKB  ##
################

# match_SNV_OncoKB = function(snv,db) {
#     ## Given a list of SNVs searches for matching variants in OncoKB database.
#     ## input: snv is a dataframe with SNVs. Must have three columns: gene symbol, variant classification and amino_acid_change (AXXXB).
#     ## output: returns those OncoKB rows matching to the input SNVs.
#     if (nrow(snv) == 0) {
#         stop('Provided SNV file was empty!')
#     }
#     else if (ncol(snv) < 3) {
#         stop('The input SNV file should have atleast three columns!')
#     }
    
#     if (missing(db)){
#         onco <- read.csv2(paste(data.path,"/OncoKB.csv",sep=""), header=T, sep=";")
#     } else {
#         onco <- db
#     }
    
#     druggable             = data.frame(matrix(ncol=ncol(onco),nrow=0))
#     colnames(druggable)   = colnames(onco)  

#     # Gene match
#     index               = as.character(onco$Gene) %in% as.character(snv[,1])  #Database genes mutated in patient
#     druggable           = subset(onco,index)       
#                                #Subset from complete DB
#     if (nrow(druggable) > 0) {
#         druggable$Patient_variant = ""
#     } else {
#         return(data.frame())
#     }

#     rownames(druggable) = NULL
#     # Output DF where actionable variants are stored
#     output              = data.frame(matrix(ncol = ncol(druggable), nrow=0))
#     rownames(output)    = NULL
#     colnames(output)    = colnames(druggable)
#     output_counter = 1
    
#     # Keep matching variants
#     if (nrow(druggable) != 0){
#         for (i in 1:nrow(druggable)){
          
#             gene = druggable[i,"Gene"]
#             p_var_type = as.character(snv[which(as.character(snv[,1])==gene),2])         #p_var_type denotes the variant type (missense, frame shift etc) for gene at i
#             p_variants = unique(as.character(snv[which(as.character(snv[,1])==gene),3])) #p_variants denotes the amino acid change(s) for gene at i
#             for (p_variant in p_variants){
              
#               ### Any mutation
#                 if (grepl("mutation|mutations", druggable[i,"Variant"], ignore.case = TRUE)){
#                     output= rbind(output, druggable[i,])
#                     output[output_counter,"Patient_variant"] = p_variant
#                     output_counter = output_counter + 1
#                 }
                

#                 if (grepl("loss|lof|fs",druggable[i,"variant"], ignore.case = TRUE)
#                   && grepl("splice|shift",p_var_type,ignore.case = TRUE)){
#                     druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
#                 }
#                 ### Positional SNV
#                 #if (grepl("A-Z][0-9]+[A-Z]",druggable[i,"Variant"], ignore.case = TRUE) && grepl("ins|del",p_var_type,ignore.case = TRUE) ){
#                 #  druggable[i,"Patient_variant"] = paste(druggable[i,"Patient_variant"],p_variant)
#                 #}
#                 # New Regex to find also comma-separated SNVs in a cell
#                 if (grepl('^(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)*(?:[a-zA-Z]*\\d+[a-zA-Z]*)$', druggable[i,"Variant"], ignore.case = TRUE)){
#                     output= rbind(output, druggable[i,])
#                     output[output_counter,"Patient_variant"] = paste(druggable[i,"Variant"],p_variant)
#                     output_counter = output_counter + 1
#                 }
                
#                 ### Stopgain ~ Deletion
#                 if (grepl("deletion|loss",druggable[i,"Variant"], ignore.case = TRUE) && grepl("nonsense",p_var_type,ignore.case = TRUE) ){
#                     output= rbind(output, druggable[i,])
#                     output[output_counter,"Patient_variant"] = p_variant
#                     output_counter = output_counter + 1
#                 }
#                 ### Missense & Specific variant
#                 if (grepl("missense",p_var_type[match(p_variant,p_variants)],ignore.case = TRUE)){
#                     if (grepl(p_variant,druggable[i,"Variant"])){
#                         output= rbind(output, druggable[i,])
#                         output[output_counter,"Patient_variant"] = p_variant
#                         output_counter = output_counter + 1
#                     }
#                     else if (grepl("[A-Z][0-9]+[A-Z]",druggable[i,"Variant"])) { #Unknown, but specific variant
#                         output= rbind(output, druggable[i,])
#                         output[output_counter,"Patient_variant"] = paste(druggable[i,"Variant"], paste("NEW", p_variant))
#                         output_counter = output_counter + 1
#                     }
#                     else if (grepl('^(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)(?:[a-zA-Z]*\\d+[a-zA-Z]*,? ?)*(?:[a-zA-Z]*\\d+[a-zA-Z]*)$|delins|indel|exon', druggable[i,"Variant"], ignore.case = TRUE)){
#                         output= rbind(output, druggable[i,])
#                         output[output_counter,"Patient_variant"] = paste(druggable[i,"Variant"], paste("NEW", p_variant))
#                         output_counter = output_counter + 1
#                     }
#                 }
#             }
#         }
#     }

#     if (nrow(output)>0) {
#         output$database_origin ="OncoKB"
#     }

#     return(output)
# }

# match_CNV_OncoKB = function(cnv, db) {
#   ## Given a list of CNVs searches for matching variants in CIViC database.
#   ## input: cnv is a dataframe with CNVs. Must have two columns: gene symbol, variant (amplification or deletion).
#   ## output: returns those CIViC rows matching to the input CNVs.
#   if (nrow(cnv) == 0) {
#     stop('Provided CNV file was empty!')
#   }
#   else if (ncol(cnv) < 2) {
#     stop('The input CNV file should have atleast two columns!')
#   }
  
#   if (missing(db)){
#     onco=read.csv2(paste(data.path,"/OncoKB.csv",sep=""),header=T, sep=";")
#   } else {
#     onco=db
#   }
  
#   druggable             = data.frame(matrix(ncol=ncol(onco),nrow=0))
#   colnames(druggable)   = colnames(onco)
  
#   # Gene match
#   index               = as.character(onco$Gene) %in% as.character(cnv[,1])  #Database genes mutated in patient
#   druggable           = subset(onco,index)                                  #Subset from complete DB
#   if (nrow(druggable) > 0) {
#     druggable$Patient_variant = ""
#   } else {
#     return(data.frame())
#   }
#   rownames(druggable) = NULL
#   # Output DF where actionable variants are stored
#   output              = data.frame(matrix(ncol = ncol(druggable), nrow=0))
#   rownames(output)    = NULL
#   colnames(output)    = colnames(druggable)
#   output_counter = 1
  
#   # Keep "amplification" or "deletion" variants
#   index     = grepl("amplification|deletion",druggable[,"Variant"],ignore.case=TRUE)
#   druggable = druggable[sort(unique(index)),]
#   rownames(druggable) = NULL
  
#   # Keep matching patient variants
#   index=character()
#   if (nrow(druggable) != 0){
#     for (i in 1:nrow(druggable)){
#       gene = as.character(druggable[i,"Gene"])
#       pvar = as.character(cnv[which(cnv[,1]==gene),2])[1]
      
#       # Get "mutation or copy number loss"
#       if (grepl("deletion",druggable[i,"Variant"],ignore.case=TRUE) && pvar == "deletion"){
#         output = rbind(output, druggable[i,])
#         output[output_counter,"Patient_variant"] = pvar
#         output_counter = output_counter + 1
#       }
#       if (grepl("amplification",druggable[i,"Variant"],ignore.case=TRUE) && pvar == "amplification"){
#         output = rbind(output,druggable[i,])
#         output[output_counter,"Patient_variant"] = pvar
#         output_counter = output_counter + 1
#       }
#     }
#   }
#   if (nrow(output)>0) {
#     output$database_origin ="OncoKB"
#   }
#   return(output)
# }

# match_TX_OncoKB = function(tx, db) {
#   ## Given a list of gene rearrangements, searches for matching variants in Gene Drug Knowledge Database (Dienstmann et al., 2015).
#   ## input: tx is a vector with gene rearrangements. Each gene rearrangement is defined by two "="-separated genes "gene1=gene2"
#   if (length(tx) == 0) {
#     stop('Provided fusion file was empty!')
#   }
#   for (element in tx) {
#     if (!grepl('=', element)) {
#       warning('Fusion file might contain wrongly formatted lines.')
#     }
#   }
  
#   if (missing(db)){
#     onco=read.csv2(paste(data.path,"/OncoKB.csv",sep=""),header=T, sep=";")
#   } else {
#     onco=db
#   }
  
#   druggable             = data.frame(matrix(ncol=ncol(onco),nrow=0))
#   colnames(druggable)   = colnames(onco)
  
#   # Retain "fusion gene" or "rearrangement" variants
#   index               = grepl("fusion",onco[,"Variant"], ignore.case = TRUE)
#   druggable           = onco[index,]
#   if (nrow(druggable) > 0) {
#     if (nrow(druggable) > 0) {
#     druggable$Patient_variant = ""
#   } else {
#     return(data.frame())
#   }
#     rownames(druggable) = NULL
    
#     output              = data.frame(matrix(ncol = ncol(druggable), nrow=0))
#     rownames(output)    = NULL
#     colnames(output)    = colnames(druggable)
#     output_counter = 1
  
#   for (fusion in tx){
#     genes  = strsplit(as.character(fusion),"=")[[1]]
#     for (i in 1:nrow(druggable)) {
#       if (grepl(genes[1],druggable[i, "Gene"]) || grepl(genes[2],druggable[i, "Gene"])) {
#         if (druggable[i,"Variant"] == "Fusions") {
#           output = rbind(output,druggable[i,])
#           output[output_counter,"Patient_variant"] = fusion
#           output_counter = output_counter + 1
#         }
#         else if (grepl(".*-.*", druggable[i,"Variant"])) { #Precise fusion
#           precise_fusions = strsplit(as.character(druggable[i,"Variant"]), "-")[[1]]
#           if ((grepl(genes[1],precise_fusions[1],ignore.case=TRUE) && grepl(genes[2],precise_fusions[2],ignore.case=TRUE)) || (grepl(genes[2],precise_fusions[1],ignore.case=TRUE) && grepl(genes[1],precise_fusions[2],ignore.case=TRUE))) {
#             output = rbind(output,druggable[i,])
#             output[output_counter,"Patient_variant"] = fusion
#             output_counter = output_counter + 1
#           }
#         }
#       }
#     }
#   }
#   }
#   if (nrow(output)>0) {
#     output$database_origin ="OncoKB"
#   }
#   return(output)
# }


# if(debug == TRUE){
#     data.path <- "/mnt/c/Users/kevin/OneDrive/Desktop/Git-Repositorys/mtb-report/data"

#     snv = matrix(c('JAK2','ALK', 'EGFR', 'BRAF', "BRAF", 'BCL2', 'PDGFRA', 'Missense_Mutation','Missense_Mutation', 'Missense_Mutation', 'Missense_Mutation',
#     'Missense_Mutation', 'Missense_Mutation', 'Missense_Mutation','V617F', 'F1174L', 'T790M', 'V600', 'V600E', 'P59T', 'L747S'), ncol=3)
#     colnames(snv) = c('Gene', 'Variant', 'Protein_Change')
#     snv = as.table(snv)

#     snv = matrix(c('ATM', 'BRAF', 'CREBBP', 'KMT2D', 'MAPK15', 
#                     'frame_shift_del', 'missense_mutation', 'frame_shift_ins', 'frame_shift_ins', 'frame_shift_del', 
#                     'L1328fs', 'R726H', 'T848fs', 'P1110fs', 'L285fs'), ncol=3)
#     colnames(snv) = c('Gene', 'Variant', 'Protein_Change')
#     snv = as.table(snv)

#     snv = matrix(c('B2M',
#                     'missense',
#                     'L1328fs'), ncol=3)



#     colnames(snv) = c('Gene', 'Variant', 'Protein_Change')
#     snv = as.table(snv)


#     #fusion = matrix(c('JAK2=GOLGA5','ALK=NPM'), ncol=1)

#     snv_gdkd = match_SNV_GDKD(snv)
#     snv_civic = match_SNV_CIVIC(snv)
#     snv_oncokb = match_SNV_OncoKB(snv)

#     cnv = matrix(c('CD274',
#                 'amplification'), ncol=2)
    
#     cnv_gdkd = match_CNV_GDKD(cnv)
#     cnv_civic = match_CNV_CIVIC(cnv)
#     cnv_oncokb = match_CNV_OncoKB(cnv)
# }


