
setwd('./')
helpers.path <- "../helpers"
data.path <- "../data"

source(paste(helpers.path,"/get_druggable.r",sep=""))
source(paste(helpers.path,"/get_levels.r",sep=""))
source(paste(helpers.path,"/get_input.r",sep=""))
source(paste(helpers.path,"/variant_annotation.r",sep=""))

# GDKD
gdkd_db_path = paste(data.path,"/GDKD.csv",sep="")
if (file.exists(gdkd_db_path)){
  print("gdkd loaded")
  #GDKD = read.delim(gdkd_db_path,sep="\t")
  GDKD <- read.table(gdkd_db_path, sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,fill=FALSE)
  #used_db_data[["GDKD"]] <- GDKD
} else {
  print("gdkd databases not found")
  #used_db_data[["GDKD"]] <- NULL
}

# CIVIC
civic_db_path = paste(data.path,"/CIViC.csv",sep="")
if (file.exists(civic_db_path)){
  print("civic loaded")
  #CIVIC = read.delim(civic_db_path, sep="\t")
  #used_db_data[["CIVIC"]] <- CIVIC
  CIVIC <- read.table(civic_db_path, sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)
  civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
  CIVIC <- civic_del(CIVIC)
  
  colnames(CIVIC) <- gsub("X.","",colnames(CIVIC))
  colnames(CIVIC) <- gsub("\\.","",colnames(CIVIC))
  print("load complete")
} else{
  print("civic db not found")
  #used_db_data[["CIVIC"]] <- NULL
}

# OncoKB
oncokb_db_path = paste(data.path,"/OncoKB.csv",sep="")
if (file.exists(oncokb_db_path)){
  print("oncokb loaded")
  ONCOKB =  read.csv2(oncokb_db_path, header=T, sep=";") #read.delim(oncokb_db_path, sep="\t")
  #used_db_data[["ONCOKB"]] <- ONCOKB
} else{
  print("oncokb not found")
  #used_db_data[["ONCOKB"]] <- NULL
}

# TARGET
target_db_path = paste(data.path,"/TARGET_MERIC.csv",sep="")
if (file.exists(target_db_path)){
  print("target loaded")
  TARGET = read.delim(target_db_path, sep="\t")
  #used_db_data[["TARGET"]] <- ONCOKB
} else{
  #used_db_data[["TARGET"]] <- NULL
  print("target not found")
  TARGET <- NULL
}


shinyServer(function(input, output, session) {

  ########################################################
  ######### TUTORIAL TAB #################################
  ########################################################
  rv <- reactiveValues(page = 1)
  observe({
    toggleState(id = "prevBtn", condition = rv$page > 1)
    toggleState(id = "nextBtn", condition = rv$page < NUM_PAGES)
    hide(selector = ".page")
    show(paste0("step", rv$page))
  })
  
  navPage <- function(direction) {
    rv$page <- rv$page + direction
  }
  
  observeEvent(input$prevBtn, navPage(-1))
  observeEvent(input$nextBtn, navPage(1))
  ########################################################
  ######### UPLOAD TAB ###################################
  ########################################################
  
  # Button to go to "upload" tab
  observeEvent(input$upload_button, {
    updateTabItems(session, "tab", "upload")
  })
  
  # Button to go to "tcga" tab
  observeEvent(input$tcga_button, {
    updateTabItems(session, "tab", "tcga")
  })
  
  ###############
  ## REACTIVES ##
  ###############
  
UPLOAD1 <- reactive({
  req(input$inputSNV)
  typeof(input$inputSNV)
  txDBPath <- paste(data.path,"/Homo_sapiens.GRCh38.100.gff3_txdb.sqlite",sep="")     # like in localmaster.r
  inFile <- input$inputSNV
  if (is.null(inFile)){print("input null");return(NULL)}
  if (grepl(".vcf$", inFile$datapath, ignore.case = TRUE)) {
    parse_vcf(inFile$datapath)
    print("parse vcf")
    VCF_to_SNV(inFile$datapath,txDBPath)          #### vcf parse function from variant_annotation.r
  } else if (grepl(".maf$", inFile$datapath, ignore.case = TRUE)) {
    #parse_maf(inFile$datapath)
    MAF_to_SNVTable(inFile$datapath,txDBPath)       #### maf parse function from variant_annotation.r
  } else {
    read.csv(inFile$datapath, sep=input$sep,  header=input$header1)
  }
})

UPLOAD2 <- reactive({
  inFile2 <- input$inputCNV
  if (is.null(inFile2)){return(NULL)}
  read.csv(inFile2$datapath,  sep=input$sep2,header=input$header2)
})

UPLOAD3 <- reactive({
  inFile3 <- input$inputTX
  if (is.null(inFile3)){return(NULL)}
  read.csv(inFile3$datapath,sep=input$sep3,col.names = c("Gene1","Gene2"), header=input$header3)
})

# Load database data
UPLOAD_DB <- reactive({
  used_db_data <- list()
  
  v_gdkd <- "gdkd"
  v_civic <- "civic"
  v_oncokb <- "onco"
  v_target <- "target"
  
  #GDKD
  gdkd_db_path = paste(data.path,"/GDKD.csv",sep="")
  if (file.exists(gdkd_db_path)){
    print("gdkd loaded")
    #GDKD = read.delim(gdkd_db_path,sep="\t")
    GDKD <- read.table(gdkd_db_path, sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,fill=FALSE)
    used_db_data[["GDKD"]] <- GDKD
  } else {
    print("gdkd databases not found")
    used_db_data[["GDKD"]] <- NULL
  }
  
  #if (input$gdkd=='no'){GDKD = 'no';v_gdkd="not used"}
  #if (input$gdkd=='default'){GDKD = read.delim("../data/GDKD.csv",sep="\t");v_gdkd="v20.0"}
  #if (input$gdkd=='other'){
  #  if (is.null(input$file_gdkd)){return(NULL)}
  #  GDKD=read.xlsx2(input$file_gdkd$datapath, header=TRUE, sheetIndex = 1,colIndex = c(1:45))
  #  # Small changes
  print("uploaddb 0")
  GDKD$Gene  = gsub(" ", "", GDKD$X.Gene.)
  for (c in grep("PMID",colnames(GDKD))){
    GDKD[,c] = gsub("\\^|_","",GDKD[,c])
    gdkd_agg   = aggregate(GDKD, by = list(GDKD$X.Disease.,GDKD$X.Gene.,GDKD$X.Description.,GDKD$X.Association_1.,GDKD$X.Therapeutic.context_1.),
                           FUN = function(X) paste(unique(X), collapse=", "))[,6:50]
    #    GDKD=gdkd_agg
    #    v_gdkd=input$file_gdkd$name
    #  }
    #}
    
    #CIVIC
    civic_db_path = paste(data.path,"/CIViC.csv",sep="")
    if (file.exists(civic_db_path)){
      #CIVIC = read.delim(civic_db_path, sep="\t")
      CIVIC <- read.table(civic_db_path, sep="\t", header=TRUE, comment.char="#",na.strings=".", stringsAsFactors=FALSE,quote="", fill=FALSE)
      civic_del <- colwise(function(x) str_replace_all(x, '\"', ""))
      CIVIC <- civic_del(CIVIC)
      
      colnames(CIVIC) <- gsub("X.","",colnames(CIVIC))
      colnames(CIVIC) <- gsub("\\.","",colnames(CIVIC))
      
      used_db_data[["CIVIC"]] <- CIVIC
    } else{
      used_db_data[["CIVIC"]] <- NULL
    }
    
    print("uploaddb 1")
    #if (input$civic=='no'){CIVIC = 'no';v_civic="not used"}
    #if (input$civic=='default'){CIVIC = read.delim("../data/CIViC.csv",sep="\t");v_civic="1 June 2019"}
    #if (input$civic=='other'){
    #  if (is.null(input$file_civic)){return(NULL)}
    #  CIVIC = read.delim(input$file_civic$datapath, header=T,stringsAsFactors = F,sep="\t",quote = "")
    #  CIVIC <- CIVIC[which(CIVIC$evidence_type=="Predictive"),]
    civic_agg  = aggregate(CIVIC, by = list(CIVIC$X.gene.,CIVIC$X.disease.,CIVIC$X.drugs.,CIVIC$X.evidence_level.,CIVIC$X.clinical_significance.),
                           FUN = function(X) paste(unique(X), collapse=", "))[,6:41]
    #  CIVIC=civic_agg
    #  v_civic=input$file_civic$name
    
    #}
    print("uploaddb 1-1")
    
    #OncoKB
    oncokb_db_path = paste(data.path,"/OncoKB.csv",sep="")
    if (file.exists(oncokb_db_path)){
      ONCOKB = read.delim(oncokb_db_path, sep="\t")
      used_db_data[["ONCOKB"]] <- ONCOKB
    } else{
      used_db_data[["ONCOKB"]] <- NULL
    }
    
    # if (input$oncokb=='no'){ONCOKB = 'no';v_oncokb="not used"}
    # if (input$oncokb=='default'){ONCOKB = read.csv2("../data/OncoKB.csv",sep=";",header=T);v_oncokb="v1"}
    # if (input$oncokb=='other'){
    #   if (is.null(input$file_oncokb)){return(NULL)}
    #   ONCOKB = read.csv2(input$file_oncokb$datapath, header=T,stringsAsFactors = F,sep=";",quote = "")
    #   ONCOKB_agg  = aggregate(ONCOKB, by = list(ONCOKB$Gene,ONCOKB$Disease,ONCOKB$Drugs,ONCOKB$Evidence_level),
    #                          FUN = function(X) paste(unique(X), collapse=", "))
    #   ONCOKB = ONCOKB_agg
    #   v_oncokb=input$file_oncokb$name
    #  
    # }
    
    #TARGET-MERIC
    #if (input$target=='no'){TARGET = 'no';v_target="not used"}
    #if (input$target=='default'){TARGET = read.delim("../data/TARGET_MERIC.csv",sep="\t");v_target="v3"}
    print("uploaddb2")
    target_db_path = paste(data.path,"/TARGET_MERIC.csv",sep="")
    if (file.exists(target_db_path)){
      TARGET = read.delim(target_db_path, sep="\t")
      used_db_data[["TARGET"]] <- ONCOKB
    } else{
      used_db_data[["TARGET"]] <- NULL
      TARGET <- NULL
    }
    
    # returns a list of used databases with the keys: "GDKD, CIVIC ,ONCOKB, TARGET"
    return(list(GDKD,CIVIC,ONCOKB,TARGET,v_gdkd,v_civic,v_oncokb,v_target))
  }
})

## Druggable
 DruggableUPLOAD = eventReactive({
    input$action
    input$sortA
  },{
    if (input$action > 0 && (!is.null(UPLOAD1()) || !is.null(UPLOAD2())|| !is.null(UPLOAD3()) ) ) {
      #Cancer type
      print("daction 0")
      cancer <- input$cancer
      cancer_GDKD <- unique(as.character(na.omit(synonyms[grep(cancer,synonyms$Shiny),"knowledge"])))
      cancer_CIVIC<- sapply(cancer,function(x) as.character(na.omit(synonyms[grep(x,synonyms$Shiny,ignore.case = T),"civic"])[1]))
      cancer_CIVIC <- paste(cancer_CIVIC,collapse = ",")
      cancer_CIVIC <- unique(strsplit(cancer_CIVIC,",")[[1]])
      cancer_ONCO <- unique(as.character(synonyms[grep(cancer,synonyms$oncokb,ignore.case = T),"oncokb"]))
      
      # Sort by
      sort <- input$sortA
      print("daction 1")
      
      #### Parse Uploads
      SNV <- UPLOAD1()#[,1:3]
      CNV <- UPLOAD2()#[,1:2]
      
      if (!is.null(UPLOAD2()[,1:2])) {
        TX <- UPLOAD3()[,1:2]
        if (!is.null(TX)){TX <- apply(TX,1,paste,collapse="=")}
      }
      print("daction 2")
      
      #used_db =  UPLOAD_DB()
      #print(used_db)
      #if(!is.null(used_db[["GDKD"]])){
      #  GDKD <- used_db[["GDKD"]]
      #}
      
      #if(!is.null(used_db[["CIVIC"]])){
      #  CIVIC <- used_db[["CIVIC"]]
      #}
      print("daction 2-1")
      #GDKD <- UPLOAD_DB()[[1]]
      #CIVIC <- UPLOAD_DB()[[2]]
      #ONCOKB <- UPLOAD_DB()[[3]]
      #TARGET <- UPLOAD_DB()[[4]]
      
      print("daction 3")
      
      if (is.null(SNV)) {
        SNV <- data.frame()
      }
      if (is.null(CNV)) {
        CNV <- data.frame()
      }
      if (!exists("TX")) {
        TX <- c()
      }
      print("daction 4")
      
      #### GDKD DB
      druggableGDKD = data.frame()
      #if(!is.null(used_db[["GDKD"]])){
        if (nrow(SNV) >=1) {
          print("daction4-1 match snv")
          druggableGDKD = match_SNV_GDKD(SNV,GDKD)
        }
        if (nrow(CNV) >=1) {
          druggableGDKD = rbind(druggableGDKD,match_CNV_GDKD(CNV,GDKD))
        }
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
          druggableGDKD = rbind(druggableGDKD,match_WT_GDKD(SNV,CNV,cancer_GDKD,GDKD))
        }
        if (length(TX) >= 1) {
          druggableGDKD = rbind(druggableGDKD,match_TX_GDKD(TX,GDKD))
        }
        rownames(druggableGDKD) = NULL
      #}
      print("result gdkd: ")
      print(druggableGDKD)
      print("daction 5")
      
      #### CIVIC
      #if(!is.null(used_db[["CIVIC"]])){
      druggableCIVIC = data.frame()
        if (nrow(SNV) >=1) {
          druggableCIVIC = match_SNV_CIVIC(SNV, CIVIC)
        }
        if (nrow(CNV) >= 1) {
          druggableCIVIC = unique(rbind(druggableCIVIC, match_CNV_CIVIC(CNV, CIVIC)))
        }
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
          druggableCIVIC = unique(rbind(druggableCIVIC, match_WT_CIVIC(SNV, CNV, cancer_CIVIC, CIVIC)))
        }
        rownames(druggableCIVIC) = NULL
      #}
      
      print("daction 5-1")
      #### OncoKB
      druggableONCO = data.frame()
      #if (!is.null(used_db[["ONCOKB"]])) {
        if (nrow(SNV) >=1) {
          druggableONCO = match_SNV_OncoKB(SNV,ONCOKB)
        }
        if (nrow(CNV) >=1) {
          druggableONCO = rbind(druggableONCO,match_CNV_OncoKB(CNV,ONCOKB))
        }
        if (length(TX) >= 1) {
          druggableONCO = rbind(druggableONCO,match_TX_OncoKB(TX,ONCOKB))
        }
        rownames(druggableONCO) = NULL
      #}
      
      print("daction 6")
      #### TARGET DB
      druggableTARGET = data.frame()
      #if (input$target != "no") {
         if (nrow(SNV) >= 1 || nrow(CNV) >= 1 || length(TX) >= 1) {
           druggableTARGET = match_TARGET_MERIC(SNV, CNV, TX, TARGET)
         }
      #}
      
      #########################################################
      ## Classify filtered variants by Levels of Evidence    ##
      #########################################################
      
      ### KNOwLEDGE
      levelsGDKD = c()
      if (nrow(druggableGDKD) > 0) {
        levelsGDKD = get_levels_GDKD(druggableGDKD, cancer_GDKD)
      }
      print("daction 6-1")
      ### CIVIC
      levelsCIVIC = c()
      if (nrow(druggableCIVIC) > 0) {
        levelsCIVIC = get_levels_CIVIC(druggableCIVIC, cancer_CIVIC)
      }
      print("daction 6-2")
      ### ONCOKB
      levelsONCO = c()
      if (nrow(druggableONCO > 0)) {
        levelsONCO = get_levels_ONCO(druggableONCO, cancer_ONCO)
      }
      
      if (length(levelsCIVIC) > 0 && length(levelsGDKD) > 0) {
        levels = merge_levels(levelsGDKD, levelsCIVIC)
      } else if (length(levelsCIVIC) > 0) {
        levels = levelsCIVIC
      } else if (length(levelsGDKD) > 0) {
        levels = levelsGDKD
      } else {
        levels = c()
      }
      print("daction 7")
      
      # Homogeneize/clean the final table
      # If OncoKB is to be taken out, replace levelsONCO with an empty vector here. This will not break the code and instead just ignore OncoKB input.
      table = clean_levels(levels,levelsONCO,synonyms,sort_by=sort)
      #print(table)
      print("dim target")
      print(dim(druggableTARGET))
      
      print("daction 8")
      return(list(druggableTARGET,druggableGDKD,druggableCIVIC,druggableONCO, table, cancer_GDKD))
    } 
    else{return()}
})
  
## Formatting output tables
outUPLOAD <- reactive({
    print("outupload 0")
    L       = c(input$LevelAUPLOAD,input$LevelBUPLOAD)
    G       = input$actionableSNVs_rows_selected    
    G2      = input$actionableCNVs_rows_selected   
    gs      = c(G,G2,L)
    
    validate(need(try(gs),'Please select gene(s) or level(s)'))
    g=c()
    g2=c()
    rowsL=c()
    rowsG=c()
    rowsG2=c()
    
    
    #gets a list with: (1 druggableTarget), 1 druggableGDKD, 2 druggableCIVIC, 3 druggableONCO, 4 table, 5 cancerGDKD
    druggableUpload <- DruggableUPLOAD()
    
    
    if(length(G)!=0){
      g       = summary_actionable_genes(DruggableUPLOAD()[[5]],UPLOAD1(),UPLOAD2(),DruggableUPLOAD()[[1]])[[1]][G,1]
      patternG= paste(g,collapse = "|")
      rowsG   = grep(patternG,DruggableUPLOAD()[[5]]$Gene)
    }
    if(length(G2)!=0){
      g2       = summary_actionable_genes(DruggableUPLOAD()[[5]],UPLOAD1(),UPLOAD2(),DruggableUPLOAD()[[1]])[[2]][G2,1]
      patternG2= paste(g2,collapse = "|")
      rowsG2   = grep(patternG2,DruggableUPLOAD()[[5]]$Gene)
    }
    if(length(L)!=0){
      pattern = paste(L,collapse = "|")
      rowsL    = grep(pattern,DruggableUPLOAD()[[5]]$level)
    }
    
    rows    = unique(c(rowsL,rowsG,rowsG2))
    if(length(rows)!=0){
      df           = DruggableUPLOAD()[[5]][rows,]
      df$Pat_Var = gsub("NEW","<font color='red'>",df$Pat_Var)
      Links        = geneLink(df$Gene)
      df$Gene      = paste(df$Gene," (",Links,")",sep="")
      Ref          = refLink(df$Ref)
      df           = cbind(df[,1:7],Ref,df[9],df[,"level"])
      colnames(df) = colnames(DruggableUPLOAD()[[5]])
      res          = grep("resistance|no|decrease",unique(df[,"Predicts"]),value=TRUE)
      sens         = grep("resistance|no|decrease",unique(df[,"Predicts"]),value=TRUE,invert = T)
      
      datatable(df,selection="single",rownames=F,escape=F,filter = 'top',options=list(scrollX = TRUE)) %>%
        formatStyle('Predicts',
                    backgroundColor = styleEqual(c(sens,res), c(rep('#66b3ff',length(sens)),rep('#FF704D',length(res))))
        ) %>%
        formatStyle('level',
                    backgroundColor = styleEqual(c("A1","A2","A3","B1","B2","B3"), c('#1A237E','#1976D2','#81D4FA','#1B5E20','#4CAF50','#AED581'))
        )
    }
    else{return()}
})
  
  
  #############
  ## OUTPUTS ##
  #############
  output$cancerGDKD <- renderText(DruggableUPLOAD()[[6]])
  output$actionableSNVs <- DT::renderDataTable({
    print("render0")
    datatable(summary_actionable_genes(DruggableUPLOAD()[[5]],UPLOAD1(),UPLOAD2(),DruggableUPLOAD()[[1]])[[1]],
              caption = "Single Nucleotide Variants ",
              selection  = "multiple",  
              escape     = FALSE ,
              rownames   = F,
              options    = list(lengthMenu = list(c( 10, 20, 30, -1), c("10", "20", "30", "All")), 
                                pageLength = 5,sDom  = '<"top">lrt<"bottom">ip'))
  })
  output$actionableCNVs <- DT::renderDataTable({
    datatable(summary_actionable_genes(DruggableUPLOAD()[[5]],UPLOAD1(),UPLOAD2(),DruggableUPLOAD()[[1]])[[2]],
              caption = "Copy Number Variants ",
              selection  = "multiple",  
              escape     = FALSE ,
              rownames   = F,
              options    = list(lengthMenu = list(c( 10, 20, 30, -1), c("10", "20", "30", "All")), 
                                pageLength = 5,sDom  = '<"top">lrt<"bottom">ip'))
  })
  
  output$mafUPLOAD <- DT::renderDataTable({
    # if(is.null(UPLOAD1())){return(datatable(data.frame(c("No SNVs have been provided yet")),rownames=F,colnames=""))}
    print("outputmaf 0")
    validate(need(try(!is.null(UPLOAD1())),"No SNVs have been provided yet"))
    
    datatable(cbind(UPLOAD1()[,c(1,2)],
                    cosmicLink(UPLOAD1()[,1],UPLOAD1()[,3]),
                    UPLOAD1()[,setdiff(1:ncol(UPLOAD1()),1:3)]),
              colnames   = colnames(UPLOAD1()),
              selection  = "single",  
              escape     = FALSE ,
              rownames   = F,
              options    = list(lengthMenu = list(c( 10, 20, 30, -1), c("10", "20", "30", "All")), 
                                pageLength = 5))
  })
  output$cnvUPLOAD <- DT::renderDataTable({
    # if(is.null(UPLOAD2())){return(datatable(data.frame(c("No CNVs have been provided yet")),rownames=F,colnames=""))}
    validate(need(try(!is.null(UPLOAD2())),"No CNVs have been provided yet"))
    
    datatable(UPLOAD2(),
              selection  = "single",
              rownames   = F,
              options    = list(lengthMenu = list(c(5, 10, 25, -1), c("5", "10", "25", "All")),
                                pageLength = 5))
  })
  output$txUPLOAD <- DT::renderDataTable({
    # if(is.null(UPLOAD3())){return(datatable(data.frame(c("No fusions have been provided yet")),rownames=F,colnames=""))}
    validate(need(try(!is.null(UPLOAD3())),"No fusions have been provided yet"))
    
    datatable(UPLOAD3(),
              selection  = "single",
              rownames   = F,
              options    = list(lengthMenu = list(c(5, 10, 25, -1), c("5", "10", "25", "All")), 
                                pageLength = 5))
  })
  
  ## Figure
  output$figureUPLOAD <- renderPlot({
    print("output figure")
    #print(DruggableUPLOAD()[[1]])
    #print(DruggableUPLOAD()[[2]])
    #print(DruggableUPLOAD()[[3]])
    #print(DruggableUPLOAD()[[4]])
    #print(DruggableUPLOAD()[[5]])
    #print(DruggableUPLOAD()[[6]])
    validate(need(try(!is.null(DruggableUPLOAD()[[5]]) &&  nrow(DruggableUPLOAD()[[5]])!=0),
                  "Sorry, no actionable variants were found."))
    
    A                     = DruggableUPLOAD()[[5]]
    df                    = expand.grid(c("1.Approved","2.Clinical","3.Preclinical"), c("B.Other\nCancers","A.Same\nCancer"))
    df$value              = c(nrow(A[A$level=="B1",]),
                              nrow(A[A$level=="B2",]),
                              nrow(A[A$level=="B3",]),
                              nrow(A[A$level=="A1",]),
                              nrow(A[A$level=="A2",]),
                              nrow(A[A$level=="A3",]))
    df$value[df$value==0] = NA
    df$value              = as.numeric(df$value)
    
    g = ggplot(df, aes(Var1, Var2)) +
      geom_point(aes(size = value), colour = c("#1B5E20","#4CAF50","#AED581","#1A237E","#1976D2","#81D4FA")) +
      theme_classic() + 
      xlab("") + 
      ylab("")
    g = g + 
      labs(title="     FINDINGS BY LEVEL")+
      scale_size_continuous(range=c(5,10)) + 
      geom_text(
        data = df, 
        aes(x = Var1, y = Var2, label = value),
        size = 5, 
        vjust = 0.5, 
        hjust = 0.5,
        colour="white",
        fontface="bold") + 
      theme(
        plot.title = element_text(size=10),
        legend.position="none",
        axis.text.x = element_text(size=10,angle = 45,hjust = 1),
        axis.text.y.right = element_text(size=10,angle = 45,vjust = 1), 
        plot.margin = margin(0.5, 0, 0, 1, "cm"),
        plot.background = element_rect(fill = "#f4f4f4",colour = "#8c8c8c",size =1))+
      scale_y_discrete(position = "right") 
    coord_fixed(ratio = 0.75)
    g
    
  })
  
  ## Target targetUPLOAD
  # output$othersUPLOAD = DT::renderDataTable({
  #   
  #   validate(need(try(!is.null(DruggableUPLOAD()[[1]]) && nrow(DruggableUPLOAD()[[1]])!=0 ), "No other actionable variants found"))
  #   
  #   datatable(DruggableUPLOAD()[[1]],
  #     selection = "single",
  #     rownames  = F,
  #     escape    = F,
  #     options=list(scrollX = TRUE)) 
  # })
  
  ## By Levels
  output$UPLOAD = DT::renderDataTable({
    print("redner level 0")
    # validate(need(try(outUPLOAD() != "" ), ""))
    validate(need(try(DruggableUPLOAD()[[5]] != ""), "No actionable variants found"))
    
    outUPLOAD()
  })
  
  ## Download PDF (this could be utilized in localmaster.r to skip the web intervace thing)
  output$reportUPLOAD = downloadHandler(
    filename = function(){
      paste(input$ID,'_MTBReport.pdf',sep='')
    },
    content  = function(file) {
      wd = getwd()
      print(wd)
      f = paste(wd,"/tmp/",sep="")
      #print(f)
      s = myFun(1)
      #print(s)
      print(dim(DruggableUPLOAD()[[5]]))
      print(dim(DruggableUPLOAD()[[1]]))
      
      report = file.path(wd, "Rnw/Report_UPLOAD-knitr.Rnw")
      print("Report thingy created")     
      invisible(knit(report, output =paste(f,s,".tex",sep="")))
      print("Knitting done")
      setwd(f)
      print("WD switched to tmp")
      tryCatch(texi2pdf(paste(f,s,".tex",sep=""),clean=T), error=function(e){write.csv(paste(e,paste(f,s,".tex",sep="")), file)})
      print("Trycatch done")
      file.copy(paste(f,s,".pdf",sep = ""), file, overwrite = TRUE)
      unlink(paste(f,s,".Rnw",sep = ""))
      unlink(paste(f,s,".pdf",sep = ""))
      setwd(wd)
    }
  )  
  
  ## Download CSV
  output$csvUPLOAD = downloadHandler(
    filename = function() {
      paste(input$ID,'_MTBReport.csv',sep='')
    },
    content = function(file) {
      write.csv(DruggableUPLOAD()[[5]], file,row.names = F)
    })
  
  ##Download Metadata
  output$metaUPLOAD = downloadHandler(
    filename = function() {
      paste(input$ID,'_Metadata.yml',sep='')
    },
    content = function(file) {
      snv_path <- "not used"
      cnv_path <- "not used"
      tx_path <- "not used"
      if (!is.null(input$inputSNV)) {
        snv_path <- input$inputSNV$datapath
      }
      if (!is.null(input$inputCNV)) {
        cnv_path <- input$inputCNV$datapath
      }
      if (!is.null(input$inputTX)) {
        tx_path <- input$inputTX$datapath
      }
      yaml_input=yaml.load(paste("data:\n cancer_type: ", input$cancer, "\n snv_path: ", snv_path, "\n cnv_path: ", cnv_path, "\n tx_path: ", tx_path, "\n date: ", Sys.Date(), sep=""))
      metafile_addendum <- yaml.load("tool:\n version: alpha\n git_checksum: undefined\ndatabases:\n gdkd: v20.0\n civic: 15-jan-01\n target-meric: v3\n oncokb: v1")
      write_yaml(append(yaml_input, metafile_addendum), file)
    })
  
  ######################################################
  ######### TCGA TAB ###################################
  ######################################################
  output$x1 = DT::renderDataTable({
    datatable(cohort_table(),
              rownames = F,
              selection = 'single', 
              options=list(lengthMenu = list(5, c("5")),pageLength=5,scrollX = TRUE)) %>% 
      formatStyle(names(cohort_table()), cursor = 'pointer')
  })
  
  ###############
  ## REACTIVES ##
  ###############
  ## Table cohort
  cohort_table = reactive({
    Samples.Clinical(format = "csv",  cohort=strsplit(input$cohort," - ")[[1]][1],page_size = 25000)[,grep("tcga_participant_barcode|histological_type|pathologic_stage|clinical_stage|acute_myeloid_leukemia_calgb_cytogenetics_risk_category",
                                                                                                           colnames(Samples.Clinical(format = "csv",  cohort=strsplit(input$cohort," - ")[[1]][1],page_size = 25000)))]
  })
  
  ## Clinical
  clinical = reactive({
    r = input$x1_rows_selected
    if (!is.null(r)) {
      patient = cohort_table()[r,1]
      clinical_data=Samples.Clinical(format = "csv",tcga_participant_barcode = patient)
      return(as.vector(clinical_data[1,]))
    }
  })
  
  ## Genomic
  TCGA = reactive({q
    print("genomic 0")
    r = input$x1_rows_selected
    if (!is.null(r)) {
      patient = cohort_table()[r,1]
      maf=try(Analyses.Mutation.MAF(format = "csv",
                                    tcga_participant_barcode = patient)[,c("Hugo_Symbol","Variant_Classification","Protein_Change")],silent=T)
      cnv = try(Analyses.CopyNumber.Genes.Thresholded(format = "csv",
                                                      tcga_participant_barcode = patient)[,c("gene","cn_alteration")],silent=T)
      if(is.data.frame(maf)){
        maf$Protein_Change=gsub("p.","",maf[,"Protein_Change"])
      } else {maf=NULL}
      if(is.data.frame(cnv)){
        cnv[cnv[,2]==-1,2]<-"low level deletion"
        cnv[cnv[,2]==1,2]<-"low level amplification"
        cnv[cnv[,2]==-2,2]<-"biallelic inactivation"
        cnv[cnv[,2]==2,2]<-"biallelic amplification"
        cnv                 = cnv[which(cnv[,2] != 0),]
        cnv                 = unique(cnv)
      } else {cnv=NULL}
      return(list(maf,cnv))
    }
  })
  
  ## Druggable
  Druggable = eventReactive(input$action2,{
    print("druggabl 0")
    r = input$x1_rows_selected
    if (!is.null(r) && (!is.null(TCGA()[[1]]) || !is.null(TCGA()[[2]]) )) {
      
      #Cancer type
      cancer      <- strsplit(input$cohort," - ")[[1]][1]
      cancer_GDKD      = unique(as.character(synonyms[grep(cancer,synonyms$tcga_cancer,ignore.case = T),"knowledge"]))
      cancer_CIVIC     = sapply(cancer,function(x) as.character(na.omit(synonyms[grep(x,synonyms$tcga_cancer,ignore.case = T),"civic"])[1]))
      cancer_CIVIC     = paste(cancer_CIVIC,collapse = ",")
      cancer_CIVIC     = unique(strsplit(cancer_CIVIC,",")[[1]])
      cancer_ONCO      = unique(as.character(synonyms[grep(cancer,synonyms$oncokb,ignore.case = T),"oncokb"]))
      
      
      #### Parse variations
      SNV <- TCGA()[[1]]
      CNV <- TCGA()[[2]]
      
      #### If any file is null, just revert it to a normal data.frame. This prevents a lot of pesky errors
      if (!exists("SNV")) {
        SNV <- data.frame()
      }
      if (!exists("CNV")) {
        CNV <- data.frame()
      }
      if (!exists("TX")) {
        TX <- c()
      }
      print("taction 0")
      
      #### GDKD DB
      druggableGDKD <- data.frame()
      if(exists("GDKD")){
        if (nrow(SNV) >=1) {
          print("taction 0-1 snv")
          druggableGDKD = match_SNV_GDKD(SNV,GDKD)
        }
        if (nrow(CNV) >=1) {
          print("taction 0-2 cnv")
          druggableGDKD = rbind(druggableGDKD,match_CNV_GDKD(CNV,GDKD))
        }
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
          print("taction 0-3 snv")
          druggableGDKD = rbind(druggableGDKD,match_WT_GDKD(SNV,CNV,cancer_GDKD,GDKD))
        }
        if (length(TX) >= 1) {
          print("taction 0-4 tx")
          druggableGDKD = rbind(druggableGDKD,match_TX_GDKD(TX,GDKD))
        }
        rownames(druggableGDKD) = NULL
      }
      print("taction 1")
      
      #### CIVIC
      druggableCIVIC = data.frame()
      if(exists("CIVIC")){
        druggableCIVIC = data.frame()
        if (nrow(SNV) >=1) {
          druggableCIVIC = match_SNV_CIVIC(SNV, CIVIC)
        }
        if (nrow(CNV) >= 1) {
          druggableCIVIC = unique(rbind(druggableCIVIC, match_CNV_CIVIC(CNV, CIVIC)))
        }
        if (nrow(SNV) >=1 || nrow(CNV) >= 1) {
          druggableCIVIC = unique(rbind(druggableCIVIC, match_WT_CIVIC(SNV, CNV, cancer_CIVIC, CIVIC))) #? julia did not use this for her analysis. Why?
        }
        rownames(druggableCIVIC) = NULL
      }
      
      print("taction 2")
      
      #### OncoKB
      druggableONCO = data.frame()
      if (nrow(SNV) >=1) {
        druggableONCO = match_SNV_OncoKB(SNV,ONCOKB)
      }
      if (nrow(CNV) >=1) {
        druggableONCO = rbind(druggableONCO,match_CNV_OncoKB(CNV,ONCOKB))
      }
      if (length(TX) >= 1) {
        druggableONCO = rbind(druggableONCO,match_TX_OncoKB(TX,ONCOKB))
      }
      rownames(druggableONCO) = NULL

      ### TARGET
      druggableTARGET = data.frame()
      #if (input$target != "no") {
      if (nrow(SNV) >= 1 || nrow(CNV) >= 1 || length(TX) >= 1) {
        druggableTARGET = match_TARGET_MERIC(SNV, CNV, TX, TARGET)
      }
      
      #### TARGET DB
      # druggableTARGET = data.frame()
      # if (nrow(SNV) >= 1 || nrow(CNV) >= 1 || length(TX) >= 1) {
      #   druggableTARGET = match_TARGET_MERIC(SNV, CNV, TX)
      # }
      # 
      # if( nrow(druggableGDKD) == 0 & nrow(druggableCIVIC) == 0 & nrow(druggableONCO) == 0 ) return(list(druggableTARGET,druggableGDKD,druggableCIVIC, druggableONCO, data.frame()))
      # 
      #########################################################
      ## Classify filtered variants by Levels of Evidence    ##
      #########################################################
      
      print("taction 3")
      
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
      levelsONCO = c()
      if (nrow(druggableONCO > 0)) {
        levelsONCO = get_levels_ONCO(druggableONCO, cancer_ONCO)
      }
      
      if (length(levelsCIVIC) > 0 && length(levelsGDKD) > 0) {
        levels = merge_levels(levelsGDKD, levelsCIVIC)
      } else if (length(levelsCIVIC) > 0) {
        levels = levelsCIVIC
      } else if (length(levelsGDKD) > 0) {
        levels = levelsGDKD
      } else {
        levels = c()
      }
      print("taction 4")
      
      # Homogenize/clean the final table
      # If OncoKB is to be taken out, replace levelsONCO with an empty vector just like above!
      table = clean_levels(levels, levelsONCO, synonyms,sort_by="drug_freq")
      print("return")
      #print(druggableGDKD)
      #print(druggableCIVIC)
      #print(druggableONCO)
      #print(table)
      #return(list(druggableGDKD,druggableCIVIC, druggableONCO, table))

      #return(list(druggableTARGET,druggableGDKD,druggableCIVIC,druggableONCO, table, cancer_GDKD))
      cancer_GDKD <- unique(as.character(na.omit(synonyms[grep(cancer,synonyms$Shiny),"knowledge"])))
      #return(list(druggableTARGET,druggableGDKD,druggableCIVIC,druggableONCO, table, cancer_GDKD))
      return(list(druggableGDKD,druggableCIVIC,druggableONCO,druggableONCO, table, cancer_GDKD))
    }
    else{return()}
  })
  
  ## Formatting output tables
  outTCGA <- reactive({
    L       = c(input$LevelA,input$LevelB)
    validate(need(try(!is.null(L)),"Select one level"))
    pattern = paste(L,collapse = "|")
    print("tcga 0")
    rows    = grep(pattern,Druggable()[[5]]$level)
    if(length(rows)!=0){
      df           = Druggable()[[5]][rows,]
      print("tcga 1")
      df$Pat_Var = gsub("NEW","<font color='red'>",df$Pat_Var)
      Links        = geneLink(df$Gene)
      df$Gene      = paste(df$Gene," (",Links,")",sep="")
      Ref          = refLink(df$Ref)
      df           = cbind(df[,1:7],Ref,df[9],df[,"level"])
      colnames(df) = colnames(Druggable()[[5]])
      res          = grep("resistance|no|decrease",unique(df[,"Predicts"]),value=TRUE)
      sens         = grep("resistance|no|decrease",unique(df[,"Predicts"]),value=TRUE,invert = T)
      
      print("tcga2")
      datatable(df,selection="single",rownames=F,escape=F,options=list(scrollX = TRUE)) %>%
        formatStyle('Predicts',
                    backgroundColor = styleEqual(c(sens,res), c(rep('#66b3ff',length(sens)),rep('#FF704D',length(res))))
        ) %>%
        formatStyle('level',
                    backgroundColor = styleEqual(c("A1","A2","A3","B1","B2","B3"), c('#1A237E','#1976D2','#81D4FA','#1B5E20','#4CAF50','#AED581'))
        )
    }
    else{return()}
  })
  
  
  #############
  ## Outputs ##
  #############
  
  output$clinical= renderUI({
    HTML(
      paste("<br><b>",names(clinical())[grep("tcga_participant_barcode|age_at|gender|metastasis|pathologic_|histological_type|pathologic_stage|clinical_stage|vital_status|acute_myeloid_leukemia_calgb_cytogenetics_risk_category|KRAS|NRAS|BRAF|PIK3",names(clinical()))],": </b>",
            clinical()[,grep("tcga_participant_barcode|age_at|gender|metastasis|pathologic_|histological_type|pathologic_stage|clinical_stage|vital_status|acute_myeloid_leukemia_calgb_cytogenetics_risk_category|KRAS|NRAS|BRAF|PIK3",names(clinical()))]      )
      
    )
  })
  print("renderoutput 01")
  output$maf= DT::renderDataTable({
    validate(need(try(!is.null(TCGA()[[1]])),"No SNVs available for this patient"))
    datatable(cbind(TCGA()[[1]][,c(1,2)],cosmicLink(TCGA()[[1]][,1],TCGA()[[1]][,3])),colnames = c("Gene","Variant Type","aa change"),
              selection="single",  escape = FALSE ,rownames=F,options=list(lengthMenu = list(c(5, 10, 25, -1), c("5", "10", "25", "All")), pageLength=5))
  })
  output$cnv= DT::renderDataTable({
    validate(need(try(!is.null(TCGA()[[2]])),"No SNVs available for this patient"))
    datatable(TCGA()[[2]],
              selection="single",rownames=F,options=list(lengthMenu = list(c(5, 10, 25, -1), c("5", "10", "25", "All")), pageLength=5))
  })
  
  ## Figure
  output$figure <- renderPlot({
    print("renderplot 0")
    validate(need(try(!is.null(Druggable()[[5]]) &&  nrow(Druggable()[[5]])!=0), "Sorry, no actionable variants were found."))
    A                     = Druggable()[[5]]
    df                    = expand.grid(c("1.Approved","2.Clinical","3.Preclinical"), c("B.Other\nCancers","A.Same\nCancer"))
    df$value              = c(nrow(A[A$level=="B1",]),nrow(A[A$level=="B2",]),nrow(A[A$level=="B3",]),nrow(A[A$level=="A1",]),nrow(A[A$level=="A2",]),nrow(A[A$level=="A3",]))
    
    df$value[df$value==0] = NA
    df$value              = as.numeric(df$value)
    
    g = ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("#1B5E20","#4CAF50","#AED581","#1A237E","#1976D2","#81D4FA")) +
      theme_classic() + xlab("") + ylab("")
    g = g + labs(title="     FINDINGS BY LEVEL")+
      scale_size_continuous(range=c(5,10)) + 
      geom_text(data = df, aes(x = Var1, y = Var2, label = value), 
                size = 5, vjust = 0.5, hjust = 0.5,colour="white", fontface="bold") + 
      theme(plot.title = element_text(size=10),legend.position="none",axis.text.x = element_text(size=10,angle = 45,hjust = 1),
            axis.text.y.right = element_text(size=10,angle = 45,vjust = 1), plot.margin = margin(0.5, 0, 0, 1, "cm"),
            plot.background = element_rect(fill = "#f4f4f4",colour = "#8c8c8c",size =1))+scale_y_discrete(position = "right") 
    coord_fixed(ratio = 0.75)
    g
    
  })
  
  
  ## Target 
  # output$others = DT::renderDataTable({
  #   
  #   validate(need(try(!is.null(Druggable()[[1]]) && nrow(Druggable()[[1]])!=0 ), "No other actionable variants found"))
  #   
  #   datatable(Druggable()[[1]],
  #     selection = "single",
  #     rownames  = F,
  #     escape    = F,
  #     options=list(scrollX = TRUE)) 
  # })
  
  ## By Levels
  output$resultsTCGA = DT::renderDataTable({
    # validate(need(try(outUPLOAD() != "" ), ""))
    print("render data0")
    validate(need(try(Druggable()[[5]] != ""), "No actionable variants found"))
    
    outTCGA()
  })
  
  ## Report PDF
  output$report = downloadHandler(
    filename =  function() {
      paste("MTBReport", ".pdf", sep="")},
    content = function(file) {
      wd = getwd()
      f = paste(wd,"/tmp/",sep="")
      s=myFun(1)
      report <- file.path(wd, "/Rnw/Report_tcga-knitr.Rnw")
      knit(report,  output =paste(f,s,".tex",sep=""))
      setwd(f) 
      texi2pdf(paste(f,s,".tex",sep=""),clean=T)
      file.copy(paste(f,s,".pdf",sep=""), file, overwrite = TRUE)
      #report <- "/home/jperera/shinyapps/app/Report_tcga.Rnw"
      # Sweave2knitr(report,output = paste(f,s,".Rnw",sep=""))
      # setwd(f)     
      # tryCatch(knit(paste(f,s,".Rnw",sep=""),  output =paste(f,s,".tex",sep="") ), error=function(e){write.csv(paste(e,paste(f,s,".tex",sep="")), file)})
      # setwd(f)
      # tryCatch(texi2pdf(paste(f,s,".tex",sep=""),clean=T), error=function(e){write.csv(paste(e,paste(f,s,".tex",sep="")), file)})
      # file.copy(paste(f,s,".pdf",sep=""), file, overwrite = TRUE)
      unlink(paste(f,s,".Rnw",sep=""))
      unlink(paste(f,s,".tex",sep=""))
      unlink(paste(f,s,".pdf",sep=""))
      setwd(wd)
    }
  )
  
  ## Report CSV
  output$csv = downloadHandler(
    filename = paste(input$ID,'_MTBReport.csv',sep=''),
    content = function(file) {
      print("write 0")
      write.csv(Druggable()[[5]], file,row.names = F)
    }
  )
  
  ##Download Metadata
  ##For the TCGA Tab, I'm not sure how to handle this yet - put into TODO
  output$metadata = downloadHandler(
    filename = function() {
      paste(input$ID,'_Metadata.yml',sep='')
    },
    content = function(file) {
      write_yaml("data:\n tcga_metadata: not yet available")
    })


  output$downloadImage <- downloadHandler(
    filename = function() {
      paste("mtbreport.tar","",sep="")
     },
     content = function(file) {
       file.copy("mtbreport.tar", file)
     }
  )

})

##increase maximum upload size to 100M
options(shiny.maxRequestSize=100*1024^2)  



#output$mafUPLOAD <- DT::renderDataTable({
#  # if(is.null(UPLOAD1())){return(datatable(data.frame(c("No SNVs have been provided yet")),rownames=F,colnames=""))}
#  validate(need(try(!is.null(UPLOAD1())),"No SNVs have been provided yet"))
#  
#  datatable(cbind(UPLOAD1()[,c(1,2)],
#                  cosmicLink(UPLOAD1()[,1],UPLOAD1()[,3]),
#                  UPLOAD1()[,setdiff(1:ncol(UPLOAD1()),1:3)]),
#            colnames   = colnames(UPLOAD1()),
#            selection  = "single",  
#            escape     = FALSE ,
#            rownames   = F,
#            options    = list(lengthMenu = list(c( 10, 20, 30, -1), c("10", "20", "30", "All")), 
#                              pageLength = 5))
#})


#})



