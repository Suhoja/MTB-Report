txDBPath <- "data/Homo_sapiens.GRCh38.100.gff3_txdb.sqlite"
inp_file <- "example_data/tcga_laml.maf"
cancer <- "unspecified"
inputfile <- MAF_to_SNVTable(inp_file,txDBPath) 
 
var_control_sign <- "any_var_only"
var_control_sign <- "without_wildtype"
var_control_sign <- "only_wildtype"
var_control_sign <- "specific_var"
var_control_sign <- "all"

Specific_MutTB_extract<-function(table,var_control_sign){
  if(var_control_sign == "all"){
    return(table)
  } else if(var_control_sign == "without_wildtype"){
    temp <- table[- grep("wild", table$`Known Var`,ignore.case = TRUE), ]
    return(temp)
  } else if(var_control_sign == "wildtype_only"){
    temp <- table[grep("wild", table$`Known Var`,ignore.case = TRUE), ]
    return(temp)
  } else if(var_control_sign == "specific_var"){
    temp <- table[- grep("wild|any", table$`Known Var`,ignore.case = TRUE), ]
    return(temp)
  } else if(var_control_sign == "any_var_only"){
    temp <- table[grep("any", table$`Known Var`,ignore.case = TRUE), ]
    return(temp)
  } else {
    warning(paste("The control sign set in yaml metadata file with ", yaml_input[["control"]][["var_control"]], " is invalid, all the known variants will be extracted", sep=""))
    return(table)
  }
}

yaml_input <- read_yaml("example_data/metadata_example.yml")
