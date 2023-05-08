#==============================================================#
#              CCS FUNCTIONALITY FOR THE MTB TOOL              #
#--------------------------------------------------------------#
#   These helper functions are designed to work exclusively    #
#         with the functionality of the MTB tool               #
#==============================================================#

suppressPackageStartupMessages(library(EnsDb.Hsapiens.v75, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(Biostrings, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(stringr, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(stringi, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(httr, quietly=TRUE, warn.conflicts=FALSE))
suppressPackageStartupMessages(library(ensembldb, quietly=TRUE, warn.conflicts=FALSE))

edbx_hg19 <- EnsDb.Hsapiens.v75
#edbx_hg38 <- EnsDb.Hsapiens.v86

#' Creates a genomic range object from a set of given genomic locations
#' @param range_dataframe A dataframe containing one or more sets of genomic coordinates
#' @return The according GRanges object
createGRange <- function(range_dataframe) {
  
  ### Sort dataframe ###
  range_dataframe <- range_dataframe[order(range_dataframe$seq_name, range_dataframe$seq_start),]
  
  ### Create GRange lists ###
  df_seqnames <- c()
  df_seqquant <- c()
  df_rangestart <- c()
  df_rangeend <- c()
  
  ### Fill GRange lists ###
  seq_counter <- 0
  for (i in 1:nrow(range_dataframe)) {
    
    ### Chromosome
    if (range_dataframe[i, "seq_name"] %in% df_seqnames) {
      df_seqquant[seq_counter] = df_seqquant[seq_counter] + 1
    }
    else {
      df_seqnames = append(df_seqnames, as.character(range_dataframe[i, "seq_name"]))
      seq_counter = seq_counter + 1
      df_seqquant = append(df_seqquant, 1)
    }
    
    ### Start & End
    df_rangestart = append(df_rangestart, range_dataframe[i, "seq_start"])
    df_rangeend = append(df_rangeend, range_dataframe[i, "seq_end"])
  }
  
  return(GRanges(seqnames = Rle(df_seqnames, df_seqquant), ranges = IRanges(start = df_rangestart, end = df_rangeend)))
}

#' Creates an IRange object from a set of given genomic locations
#' @param range_dataframe A dataframe containing one or more sets of genomic coordinates
#' @return The according IRanges object
createIRange <- function(range_dataframe) {
  
  ### Sort dataframe ###
  range_dataframe <- range_dataframe[order(range_dataframe$tx_id, as.integer(range_dataframe$var_start)),]
  
  ### Create IRange lists ###
  df_txnames <- c()
  df_txquant <- c()
  df_rangestart <- c()
  df_rangeend <- c()
  
  ### Fill IRange lists ###
  seq_counter <- 0
  for (i in 1:nrow(range_dataframe)) {
    
    ### Transcript ID
    if (range_dataframe[i, "tx_id"] %in% df_txnames) {
      df_txquant[seq_counter] = df_txquant[seq_counter] + 1
    }
    else {
      df_txnames = append(df_txnames, as.character(range_dataframe[i, "tx_id"]))
      seq_counter = seq_counter + 1
      df_txquant = append(df_txquant, 1)
    }
    
    ### Start & End
    df_rangestart = append(df_rangestart, as.integer(range_dataframe[i, "var_start"]))
    df_rangeend = append(df_rangeend, as.integer(range_dataframe[i, "var_end"]))
  }
  
  return(IRanges(names = df_txnames, start = df_rangestart, end = df_rangeend))
}

#' Creates an XML query in the format that ENSEMBL/biomaRt requires
#' @param dataset: The name of the biomaRt dataset the data should come from. The default is 'hsapiens_gene_ensembl'.
#' @param filter_name: The ENSEMBL term for the provided value that is given. For a transcript ID, for example, it is 'ensembl_transcript_id'
#' @param filter_value: The actual value that is to be used as a search basis. In the case of a transcript ID, for example, it could be 'ENST00000378816'
#' @param attributes: The ENSEMBL terms of the attributes that are to be queried from biomaRt
#' @return A XML Query that can be POSTed to biomaRt
createBMQuery <- function(dataset, filter_name, filter_value, attributes) {
  attrstring <- ""
  for (element in attributes) {
    attrstring <- paste(attrstring, '<Attribute name="', element, '"/>', sep="")
  }
  return(paste('query=<Query client="UMG" processor="JSON" limit="-1" header="1"><Dataset name="', dataset, '"><Filter name="', filter_name, '" value="', filter_value, '"/>"', attrstring, '</Dataset></Query>', sep="" ))
}

#' Finds the triplet start position for a given position in a sequence
#' @param pos: The sequence position
#' @return The according triplet start position
findSeqStartPosition <- function(pos) {
  if (pos %% 3 == 1) {
    return(pos)
  } else if (pos %% 3 == 2) {
    return(pos-1)
  } else {
    return(pos-2)
  }
}

#' Finds the triplet end position for a given position in a sequence
#' @param pos: The sequence position
#' @return The according triplet end position
findSeqEndPosition <- function(pos) {
  if (pos %% 3 == 1) {
    return(pos+2)
  } else if (pos %% 3 == 2) {
    return(pos+1)
  } else {
    return(pos)
  }
}

############ MAIN FUNCTIONS ###############

#' Takes a genomic location + a nucleotide exchange and returns the respective HUGO symbol plus all transcripts,
#' the exchanged amino acid and the position of the exchange in the primary transcript.
#' @param range_dataframe - A dataframe containing all the relevant data.
#' @param genome - An EnsemblDb object, either v75 or v86
#' @return HUGO symbol, primary transcript ID, position and exchanged amino acid
Coordinate_Converter<-function(range_dataframe, genome) {
  
  ### Create GRanges object ###
  genomic_range <- createGRange(range_dataframe)
  
  ### Find transcripts for GRange ###
  suppressWarnings({transcripts <- as.data.frame(genomeToTranscript(genomic_range, genome))})
  if (!"names" %in% colnames(transcripts)) {
    return(NULL)
  } 
  transcripts <- merge(range_dataframe, transcripts)
  
  ### Extract transcript information and sequence ###
  trsc_dataframe <- subset(transcripts, transcripts$start != -1 & transcripts$end != -1 & transcripts$names != "", select=c("names", "start", "end", "refAllele", "varAllele"))
  rownames(trsc_dataframe) = NULL
  if (nrow(trsc_dataframe) == 0) { return(NULL)}

  for (row in 1:nrow(trsc_dataframe)) {
    ### Convert given coordinates to CDS coordinates ###
    cds_ir <- IRanges(start=trsc_dataframe[row, "start"], end = trsc_dataframe[row, "end"], name = trsc_dataframe[row, "names"])
    cds_mapped <- transcriptToCds(cds_ir, edbx_hg19)
    trsc_dataframe[row, "start"] = start(cds_mapped)
    trsc_dataframe[row, "end"] = end(cds_mapped)
    
    ### Download sequence, strand and HGNC symbol from ensembl ###
    current_sequence_r <- POST("http://useast.ensembl.org/biomart/martservice/results",
                               body=createBMQuery('hsapiens_gene_ensembl', 'ensembl_transcript_id', trsc_dataframe[row, "names"], c("coding", "hgnc_symbol", "strand"))
    )
    information <- str_split(httr::content(current_sequence_r, "text"), "\t|\n")[[1]]
    current_sequence <- information[4]
    if (current_sequence == "" || current_sequence == "Sequence unavailable" || trsc_dataframe[row, "start"] == -1 || trsc_dataframe[row, "end"] == -1) {
      #trsc_dataframe[row, "comment"] <- "No sequence found for this transcript OR no coordinates found"
      trsc_dataframe[row, "aa_change"] <- ""
      trsc_dataframe[row, "hgnc_symbol"] <- ""
      trsc_dataframe[row, "strand"] <- ""
      next 
    }
    alt_sequence <- information[4]
    gene <- information[5]
    
    if (length(information) >= 6) { # Strand is contained
      if (information[6] == "-1") {
        strand <- "-"
        variation <- as.character(complement(DNAString(trsc_dataframe[row, "varAllele"]))) # On minus strand, alleles have to be their complements
      } else {
        strand <- "+"
        variation <- trsc_dataframe[row, "varAllele"]
      }
    } else {
      strand <- "+"
      variation <- trsc_dataframe[row, "varAllele"]
    }
    
    substr(alt_sequence, trsc_dataframe[row, "start"], trsc_dataframe[row, "end"]) <- variation
    suppressWarnings({
      seq1 <- translate(DNAString(substr(current_sequence, findSeqStartPosition(as.integer(trsc_dataframe[row, "start"])), findSeqEndPosition(as.integer(trsc_dataframe[row, "end"])))))
      seq2 <- translate(DNAString(substr(alt_sequence, findSeqStartPosition(as.integer(trsc_dataframe[row, "start"])), findSeqEndPosition(as.integer(trsc_dataframe[row, "end"])))))
    })
    
    ### Calculate amino acid exchange ###
    from <- ""
    difference <- ""
    for (i in 1:nchar(as.character(seq1))) {
      #message(paste("Comparing ", substr(as.character(seq1), i, i), " to ", substr(as.character(seq2), i, i), sep="" ))
      from = paste(from, substr(as.character(seq1), i, i), sep="")
      difference = paste(difference, substr(as.character(seq2), i, i), sep="")
    }
    trsc_dataframe[row, "aa_change"] <- paste(from, ">", difference, sep="")
    trsc_dataframe[row, "hgnc_symbol"] <- gene
    trsc_dataframe[row, "strand"] <- strand
  }
  return(trsc_dataframe)
}

