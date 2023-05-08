
# Paths
data.path <- "../data"

my_autocomplete_list <-sort(unique(as.character(synonyms$Shiny)))
tutorial=c("1. Clinical","2. SNVs","3. CNVs","4. Fusions","5. Databases")

metainfo_path <- paste(data.path,"/MTB_Databases_versions.csv",sep="")
metainfo <- read.csv(metainfo_path, header=TRUE, sep=",")
html_meta = ""
string_meta = ""

if(any(metainfo == "CIViC.csv")){
  civic_meta = metainfo[metainfo$Database == "CIViC.csv",]
  html_meta = paste(html_meta, 	"<li><a target='_blank' href='https://civic.genome.wustl.edu/'>", civic_meta[[1]], "</a>, ", "date: ", civic_meta[[4]], ", version: ", civic_meta[[5]] ,"</li>", sep="")
  string_meta = paste(string_meta, "CIVIC:", "https://civic.genome.wustl.edu/", "Date:", civic_meta[[4]], "Version:", civic_meta[[5]], sep = " ")
  string_meta = paste(string_meta, "\n", sep = "")
}
if(any(metainfo == "GDKD.csv")){
  gdkd_meta = metainfo[metainfo$Database == "GDKD.csv",]
  html_meta = paste(html_meta, 	"<li><a target='_blank' href='https://www.synapse.org/#!Synapse:syn2370773'>", gdkd_meta[[1]], "</a>, ", "date: ", gdkd_meta[[4]], ", version: ", gdkd_meta[[5]] ,"</li>", sep="")
  string_meta = paste(string_meta, "GDKD:", "https://www.synapse.org/#!Synapse:syn2370773", "Date:", gdkd_meta[[4]], "Version:", gdkd_meta[[5]], sep = " ")
  string_meta = paste(string_meta, "\n", sep = "")
}
if(any(metainfo == "OncoKB.csv")){
  oncokb_meta = metainfo[metainfo$Database == "OncoKB.csv",]
  html_meta = paste(html_meta, 	"<li><a target='_blank' href='https://oncokb.org/api/v1/utils/allActionableVariants.txt'>", oncokb_meta[[1]], "</a>, ", "date: ", oncokb_meta[[4]], ", version: ", oncokb_meta[[5]], "</li>", sep="")
  string_meta = paste(string_meta, "OncoKB:", "https://civic.genome.wustl.edu/", "Date:", oncokb_meta[[4]], "Version:", oncokb_meta[[5]], sep = " ")
  string_meta = paste(string_meta, "\n", sep = "")
}

sidebar <- dashboardSidebar(
	sidebarMenu(id="tab",
		tags$head(tags$style(HTML('.disclaimer_narrow {padding: 7px;background-color: #34444c;border-radius: 25px;margin: 10px;position: fixed;width: 180px;bottom: 0px;}
															 .main-sidebar .user-panel, .sidebar-menu, .sidebar-menu>li.header {white-space: normal;overflow: hidden;}'
		))), # close tags
		menuItem("Home",                  tabName = "home",     icon = icon('home')),
		#menuItem("Documentation",              tabName = "tutorial", icon = icon('info')),
		menuItem("Upload patient data",   tabName = "upload",   icon = icon('upload')),
		menuItem("Explore TCGA dataset",  tabName = "tcga",     icon = icon('navicon')),
		HTML("<p class='disclaimer_narrow'>
					<strong>Disclaimer:</strong> This resource is intended for research use only. It should not be used for medical or professional advice. We make no guarantee regarding comprehensiveness, accuracy or reliability of the content on this site. By using this site, you assume full responsibility for all risks associated with this.
					</p>"
		)
	),
	collapsed = TRUE
)

#shiny-tab-home {width:1200px;padding: 10px 80px;}

##################################
# Shiny UI Start
##################################

body <- dashboardBody(
  tags$head(tags$style(
    HTML('.h1, .h2, .h3, .h4, .h5, .h6, h1, h2, h3, h4, h5, h6 {line-height: 1.5;margin-top: 0px;}
          .h5, h5 {font-size: 14px;}
          .disclaimer {background-color: #f4f4f4;width:700px;border-radius: 25px;padding: 10px;}
          .content-wrapper {background-color: #FFF;}
          .main-header .logo {font-size:19px;}
          .main-sidebar {width: 200px;}
          .navbar { margin-left: 200px; } .content-wrapper, .main-footer, .right-side { margin-left: 200px; }
          .nav-tabs {background-color:#bfbfbf;}
          .nav-tabs-custom .nav-tabs li.active {border-top-color: #0099b4}
          .nav-tabs-custom .tab-content {padding:  10px 40px 10px 40px;background: #f4f4f4}
          .nav-tabs-custom {border: 1px solid #8c8c8c;}
          .btn-default:hover {background-color: #0099b4}
          .box {border-top: 3px solid rgba(255,255,255,.15);box-shadow: 0 1px 0px rgba(0, 0, 0, 0)}

          #x1, #x2, #mafUPLOAD, #cnvUPLOAD, #mafHDB, #cnvHDB td {line-height:70%;}
          #clinical {text-align:left;}
          #maf td {padding: 2px 6px;}
          #cnv td {padding: 2px 6px;}
          #actionableSNVs td , #actionableCNVs td {padding: 2px 6px}
          #resultsTCGA td , #others td , #UPLOAD td , #othersUPLOAD td{padding: 1px 6px;}

          .dataTables_paginate {font-size:9px;}
          .dataTables_length {display: none;}
          .dataTables_wrapper no-footer  {vertical-align: center;}
          .popover-content {
              color: gray;
              font-size: 11px;
          }'
				) #close HTML
			)
    ), #close tags head + styles 

  ##################################
  # Shiny Tab - HOME
  ##################################
  
	tabItems(
		tabItem("home",
			fluidRow(
				column(width = 9,
					box(width = 12,
						HTML("
							<center>
							<img src='Logo.jpg' style='width:120px;float: left'>
							<img src='umg-logo-colored.svg' style='width:250px;float: right'>
							<h1><strong>Molecular Tumor Board Report</strong> </h1>
							<h4><strong> - Tool to browse clinically actionable variants of cancer patients -</strong> </h1>
							<h4><i>developed by Nadine S. Kurz, J&uacutelia Perera-Bel, Charlotte Höltermann, Tim Tucholski, Jingyu Yang, Tim Bei&szligbarth and J&uuml;rgen D&ouml;nitz</h4></i>
							<hr>"),
						box(width = 12,
								actionButton("upload_button","Upload Patient Data"),
								actionButton("tcga_button","Explore TCGA dataset"),
								br(),
								br(),
								HTML("<hr>")
							),
							HTML("
								<h4>1. Upload of the patient's somatic variants (<a href='snv_legend.png' target='_blank'>SNVs</a> , <a href='cnv_legend.png' target='_blank'>CNVs</a>, <a href='tx_legend.png' target='_blank'>fusions</a>)</h4>
									<img src='man_icon.png' style='width:40px'>
									<img src='dna.png' style='width:80px'>
									<img src='upload.png' style='width:60px'>
									<br>
									<br>
								<h4>2. Public databases are browsed for clinically actionable somatic variants by the tool</h4>
									<img src='database_search.png' style='width:60px'>
									<img src='man_drug.png' style='width:40px'>
									<br>
									<br>
								<h4>3. The patient-specific results are shown in the web interface for further processing.</h4>
								<h4>Additionally, <a href='Report.pdf' target='_blank'>pdf reports</a> and other formats can be downloaded</h4>
									<img src='web_browser3.png' style='width:80px'>
									<img src='download_file2.png' style='width:60px'>
									</center>")
					) # close box 12
				), # close column 9

				column(width=3,
					verticalLayout(
						box(title="MTB Report information",width = 12,status = "primary",solidHeader = T,
							HTML("The <b>Molecular Tumor Board Report (MTB Report)</b> is a tool designed to bring together all efforts and knowledge on predictive biomarkers. The tool is oriented at identifying actionable somatic variants of a patient's genomic profile.<br> More detailed information about the method can be found in <a target='_blank' href='https://doi.org/10.1186/s13073-018-0529-2'>our publication (Perera-Bel et al., 2018).</a>")
							), # close box Info
						box(title="Supplementary Information",width = 12,status = "primary",solidHeader = T,
						 HTML(paste(
						    "Additional information and example data files can be found in the <a href='http://mtb.bioinf.med.uni-goettingen.de/downloads/'>download section</a>. 
						Download MTB Report Docker image <a target='_blank' href='http://mtb.bioinf.med.uni-goettingen.de/downloads/mtbreport.tar'>here</a>.",
						 "",sep="")))  ,box(
						 title="Current databases supported:",width=12, status="primary",solidHeader=T,
							HTML(paste("<ul>",
								"<li><a target='_blank' href='https://www.synapse.org/#!Synapse:syn2370773'>Gene Drug Knowledge Database</a>, v20.0</li>",
								"<li><a target='_blank' href='https://civic.genome.wustl.edu/'>CIViC</a>, 1. Jan. 2019</li>",
								"<li><a target='_blank' href='http://archive.broadinstitute.org/cancer/cga/target'>TARGET</a>, v3</li>",
								"<li><a target='_blank' href='https://doi.org/10.1093/jnci/djv098'>Meric-Bernstam et al., 2015</a></li>",
								"<li><a target='_blank' href='https://oncokb.org/api/v1/utils/allActionableVariants.txt'>OncoKB</a>, v1</li>",
							"</ul>",
						 "<hr>",
						 "Loaded Databases:",
							"<ul>",
						    html_meta,
						   "</ul>",
						   sep="") 
						) # close HTML
					 ) # close box News and Updates
					) # close VerticalLayout
				) # close column 3
			) # close Fluidrow
		), # close tabItem home

		
		##################################
		# Shiny Tab - Tutorial
		##################################
		
		#tabItem("tutorial",
		#	"Work in progress"
		#), # close tabItem tutorial

		
		##################################
		# Shiny Tab - Upload
		##################################
		
		tabItem("upload",
		        
			HTML("<h2><strong>Generate MTB-Report of user defined data</strong> </h2>"),
			
			#sidebarLayout(
			# sidebarPanel(width=3,
			box(width = 12,
				#HTML("<h5>1. Please upload clinical and genomic data of the patient.</h5>
				#	<br>
				#	<font color='red'>*</font>
				#	You have to specify <b>cancer type</b> and upload <b>at least one type of genomic data</b> (<i>SNVs, CNVs or fusions</i>):</h4>"),
				 # h5("Cancer type is the only information required by the method, which is used to classify actionable variants into Levels of Evidence.
				#					Any other information is merely for completeness of the final PDF report."),
			#),
			#div(margin=0,padding=0,width = 12,
				fluidRow(
					tabBox(id="uploadtab",width = 4,side="left",title=strong("Upload Data"),
						tabPanel(title="1. Clinical",
							fluidRow(
								column(width=12,
											selectInput("cancer",label=HTML("Select cancer type<font color='red'>*</font>"),choices=my_autocomplete_list,width = "250px",selected="unspecified"),
											textInput("ID",label="Type in patient or cell line ID",width = "250px"),
											radioButtons("gender",label="Select gender",choices=c("Female","Male","Other"),width = "250px",selected="none"),
											selectInput("stage",label="Select cancer stage",choices=c("not specified","I","II","III","IV"),width = "130px"),
								#),
								#column(width=6,
											 textInput("prev_ther",label="Type in previous therapies",width = "250px"),
											 radioButtons("tissue",label="Select tissue type",choices=c("Fresh Frozen","FFPE","Other"),width = "250px",selected="none"),
											 textInput("tumor_content",label="Type in tumor content (%)",width = "170px")
								)
							)
						), # close clinical
						tabPanel(title=HTML("2. SNVs"),
							HTML("<b> Upload file with <u>S</u>ingle <u>N</u>ucleotide <u>V</u>ariants (SNVs). Check CSV format <a href='snv_legend.png' target='_blank'>here</a> </b>
										<br>
										For SNVs, MAF (Mutation Annotation Format) and VCF (Variant Call Format) files are also accepted."),
							br(),
							br(),
							h5("This tool does not check the quality of the variants, so please make sure they are validated and fullfilling your own requirements"),
							fluidRow(
								column(width=12,
											fileInput("inputSNV",label=HTML("<a href='snv.csv' target='_blank'>Click to get example CSV file</a>"),width = "250px"),
											h5(em("File extensions accepted are: .xls, .xlsx, .csv, .dat, .tsv, .maf, .vcf")),
								#			),
								#column(width=3,
											radioButtons('sep', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
											checkboxInput(inputId = 'header1', label = 'Header', value = F),
								#			),
								#column(width=4,HTML('<b>Mandatory columns:</b> Gene Symbol, variant type, Aminoacid change')
											),
								HTML("")
						)), # close snv
						tabPanel(title="3. CNVs",
							HTML("<b> Upload file with <u>C</u>opy <u>N</u>umber <u>V</u>ariations (CNVs). Check CSV format <a href='cnv_legend.png' target='_blank'>here</a> </b>"),
							br(),
							br(),
							h5("This tool does not check the quality of the variants, so please make sure they are validated and fullfilling your own requirements"),
							fluidRow(
								column(width=5,
											fileInput("inputCNV",label=HTML("<a href='cnv.csv' target='_blank'>Click to get example CSV file</a>"),width = "250px"),
											h5(em("File extensions accepted are: .xls, .xlsx, .csv, .dat, .tsv"))
											),
								column(width=3,
											radioButtons('sep2', 'Separator',c(Comma=',',Semicolon=';', Tab='\t'),','),
											checkboxInput(inputId = 'header2', label = 'Header', value = F)
											),
								column(width=4,
											column(width=4,HTML('Mandatory columns: Gene Symbol, variant type'))
											)
						)), # close cnv
						tabPanel(title="4. Fusions",
							HTML("<b> Upload file with fusions. Check CSV format <a href='tx_legend.png' target='_blank'>here</a> </b>"),
							br(),
							br(),
							h5("This tool does not check the quality of the variants, so please make sure they are validated and fullfilling your own requirements"),
							fluidRow(
								column(width=8,
											fileInput("inputTX",label=HTML("<a href='fusions.csv' target='_blank'>Click to get example CSV file</a>"),width = "250px"),
											h5(em("File extensions accepted are: .xls, .xlsx, .csv, .dat, .tsv"))
											),
								column(width=4,
											radioButtons('sep3', 'Separator',c(Comma=',',Semicolon=';',Tab='\t'),','),
											checkboxInput(inputId = 'header3', label = 'Header', value = F)
											)
						))#, # close fusions
						#tabPanel(title="5. Databases",
							#h5(strong(c("Please select the databases used by the report"))),
							#br(),
							#fluidRow(
								#column(width=8, string_meta
								
								       #c("Add Used Databases. Format: Database-Name, Filename, Date, Version. Maybe with Buttons if they should be used or not.")
								# 	radioButtons('gdkd', em('Gene-Drug Knowledge Database (GDKD)'),
								# 							 c('v20.0'='default',none='no',other='other'),'default'),
								# 	conditionalPanel(
								# 		condition = "input.gdkd == 'other'",
								# 		fileInput("file_gdkd",label=h6(a(href='https://www.synapse.org/#!Synapse:syn2627707','Please upload a file from repository',target="_blank")),width = "250px")),
								# 	br(),
								# 	br(),
								# 	radioButtons('oncokb', em('PRecision Oncolocy Knowledge Base (OncoKB)'),
								# 							 c('v1'='default',none='no',other='other'),'default'),
								# 	conditionalPanel(
								# 		condition = "input.oncokb == 'other'",
								# 		fileInput("file_oncokb",label=h6(a(href='https://www.oncokb.org/','Please upload a file from repository',target="_blank")),width = "250px")
								# )),
								# column(width=4,
								# 	radioButtons('civic', em('Clinical Interpretations of Variants in Cancer (CIViC)'),
								# 							 c('1 June 2019'='default',none='no',other='other'),'default'),
								# 	conditionalPanel(
								# 		condition = "input.civic == 'other'",
								# 		fileInput("file_civic",label=h6(a(href='https://civicdb.org/releases','Please upload a file from repository',target="_blank")),width = "250px")),
								# 	br(),
								# 	br(),
								# 	radioButtons('target', em('TARGET-MERIC'),
								# 							 c('v3'='default',none='no'),'default'),
								#)
							#)
						#)
					),
					#) # close FluidRow 
				#) # close div
			# ), # close sidebar panel # close Upload Data tabBox
				# mainPanel(
				  
				  
			  tabBox(width=8,
						side="left",
						title = strong("Browse Genomic Data"),
						tabPanel(title="SNVs ",width=6,h5("Single nucleotide variants of the patient's tumor:"),DT::dataTableOutput("mafUPLOAD"),br()),
						tabPanel(title="CNVs ",h5("Copy number variants of the patient's tumor:"),DT::dataTableOutput("cnvUPLOAD"),br()),
						tabPanel(title="Fusions ",h5("Gene fusions of the patient's tumor:"),DT::dataTableOutput("txUPLOAD"),br())
					), # close Genomic Data tabBox
				#) # close mainPanel
				)
			),#, # close box 12
			conditionalPanel(condition=
				"typeof output.mafUPLOAD != 'undefined' ||
				 typeof output.cnvUPLOAD != 'undefined' ||
				 typeof output.txUPLOAD != 'undefined'",
				box(width = 12,
						HTML(" <h5>2. If the genomic data looks correct, <b>launch the MTB report</b> by clicking the button: </h5>"),
						actionButton("action","Launch/Refresh MTB report",icon = icon("refresh"))
				)
			), # close conditionalPanel to display ActionButtion
			br(),
			
			#fluidRow(
				conditionalPanel(
					condition = "(input.action > 0 )",
					box(solidHeader = T, width = 12,collapsible = F,
						h3(strong("  Browse clinically relevant genomic data of the selected patient:")),
						HTML("<h5> You can now explore the filtered variants divided into 6
									<a target='_blank' href='levels.png'>levels of evidence</a>, which determine the actionability of the variant.
									<br>
									You can also download a .pdf or .csv report with the results.
									For more information on the method we forward you to our
									<a target='_blank' href='https://doi.org/10.1186/s13073-018-0529-2'>publication.</a> </h5> "
						),
						br(),
						box(title = "Summary of Actionable Variants: Filter & Selection",width = 12,collapsible = F,status = 'success',solidHeader = T,
							column(width=6,
								p(strong( "Select Actionable Variants " )) ,
						#     p(strong( "Single Nucleotide Variants " )) ,
								DT::dataTableOutput('actionableSNVs'),
								br(),
						#    p(strong( "Copy Number Variants " )) ,
								DT::dataTableOutput('actionableCNVs')
							),
							column(width=3,
								p(strong( "Select Levels of Evidence" )),
								HTML('<font color="gray">Evidence on SAME cancer type</font>'),
								checkboxGroupInput(
												"LevelAUPLOAD","",
												c("A1) FDA & Guidelines"="A1","A2) Clinical Trials"="A2","A3) Pre-clinical"="A3"),selected="none"
												),
								HTML('<font color="gray">Evidence on OTHER cancer types</font>'),
								checkboxGroupInput(
												"LevelBUPLOAD","",
												c("B1) FDA & Guidelines"="B1","B2) Clinical Trials"="B2","B3) Pre-clinical"="B3"),selected="none"
												),
								checkboxInput("checkothersUPLOAD","Display other genes without evidence level",value = FALSE)
							),
							column(width=3,
								br(),
								plotOutput(outputId = "figureUPLOAD",height = "170px" , width = "80%"),
								conditionalPanel(
												condition=" typeof output.othersUPLOAD != 'undefined' || typeof output.figureUPLOAD != 'undefined'",
												br(),
												radioButtons('sortA',"Sort results by:",c("Genes"="genes","Drugs"="drug","Levels of Evidence"="levels"),selected="genes"),
												h5(strong("Download results")),
												downloadButton('reportUPLOAD','Download report (.pdf)'),br(),
												downloadButton('csvUPLOAD','Download report (.csv)'),br(),
												downloadButton('metaUPLOAD', 'Download metadata (.yml)')
												)
							)
						),
						br(),
						box(title = "Details of Actionable Variants",width = 12,collapsible = F,status = 'success',solidHeader = T,
							column(width=12,
								#  textOutput('cancerGDKD'),
								# Display Results
								DT::dataTableOutput('UPLOAD'),
								# Display Other Results
								conditionalPanel(
										condition = "input.checkothersUPLOAD",
										DT::dataTableOutput('othersUPLOAD')
								),
								br()
							)
						)
					) # close box with results
				) # close conditionalPanel responding to ActionButtion
			#) # close FluidRow
	
    # close mainPanel	
	#) # close tabItem
	# ) # close main Panel
  
  ), #close tabItem upload

		
		##################################
		# Shiny Tab - TGCA
		##################################
		
		tabItem("tcga",
			h2(strong("Explore the Cancer Genome Atlas Dataset")),
			br(),
			box(width = 10,
				HTML("<h4>1. You can browse up to 34 cancer cohorts. Data is dowloaded using the
							<a href='https://github.com/mariodeng/FirebrowseR' target='_blank'>FirebrowseR</a> package
							(R client for <a href='http://firebrowse.org/api-docs/' target='_blank'>Broad Firehose Web API</a>)</h4>")
			),
			box(width = 12,
				fluidRow(
					tabBox(title = strong("Select Patient"),
						tabPanel(title="TCGA",
							h5("Please select one cohort from the list below:"),
							selectInput("cohort","Select cohort",choices = paste(cohorts$V1,cohorts$V2, sep = " - " )),
							HTML("<h5>Please select one row from the table to explore its clinically relevant data:</h5>"),
							DT::dataTableOutput('x1')
						)
					), # close tabBox Select Patient
					conditionalPanel(
						condition = "input.x1_rows_selected.length >0 ",
						tabBox(
							side="left",
							title = strong("Patient Data"),
							tabPanel(title="SNVs ",h5("Single nucleotide variants of the patient's tumor:"),DT::dataTableOutput('maf')),
							tabPanel(title="CNVs ",h5("Copy number variants of the patient's tumor:"),DT::dataTableOutput('cnv')),
							tabPanel(title="Clinical",h5("Main clinical features of the patient's tumor:"),uiOutput('clinical'))
						)
					) # close ConditionalPanel display genomic
				) # close FluidRow
			), # close box 12
			conditionalPanel(
				condition="input.x1_rows_selected.length >0 && (typeof output.maf != 'undefined' || typeof output.cnv != 'undefined')",
				box(width = 10,
					HTML(" <h4>2. The tool will search in public databases for clinically relevant information. To see the results,
								launch the MTB report by pressing the button below. <br>After selecting another patient, please remember to refresh the results by
								pressing the button again. </h4>"),
					actionButton("action2","Launch/Refresh MTB report",icon = icon("refresh"))
				)
			), # close ConditionalPanel Launch MTB button
			br(),
			br(),
			br(),
			br(),
			fluidRow(
				conditionalPanel(
					condition = "input.x1_rows_selected.length >0 && input.action2 > 0 && (typeof output.maf != 'undefined' || typeof output.cnv != 'undefined')",
					box(solidHeader = T, width = 12,collapsible = F,
						h3(strong("  Browse clinically relevant genomic data of the selected patient:")),
						HTML("<h4> You can now explore the filtered variants divided into 6
									<a target='_blank' href='levels.png'>levels of evidence</a>,
									which determine the actionability of the variant.
									<br>
									You can also download a .pdf or .csv report with the results.
									For more information on the method we forward you to our
									<a target='_blank' href='https://doi.org/10.1186/s13073-018-0529-2'>publication.</a> </h4> "
								),
						br(),
						box(solidHeader = T, width = 12,collapsible = F,
							column(width=3,
								checkboxGroupInput(
									"LevelA","Evidence on SAME cancer type:",
									c("A1) FDA & Guidelines"="A1","A2) Clinical Trials"="A2","A3) Pre-clinical"="A3"),selected="none"
									),
								checkboxGroupInput(
									"LevelB","Evidence on OTHER cancer types:",
									c("B1) FDA & Guidelines"="B1","B2) Clinical Trials"="B2","B3) Pre-clinical"="B3"),selected="none"
									),
								plotOutput(outputId = "figure",height = "170px" , width = "80%"),
								br(),
								checkboxInput("checkothers","Display other genes without evidence level",value = FALSE),
								conditionalPanel(
									condition=" typeof output.other != 'undefined' || typeof output.figure != 'undefined'",
									h5(strong("Download results")),
									downloadButton('report','Download report (.pdf)'),br(),
									downloadButton('csv','Download report (.csv)'),br(),
									downloadButton('metadata', 'Download metadata (.yml)')
								)
							),
							column(width=9,
								# Display Results
								DT::dataTableOutput('resultsTCGA'),
								# Display Other Results
								conditionalPanel(
									condition = "input.checkothers",
									DT::dataTableOutput('others')
								)
							)
						)
					) # close box with results
				) # close conditionalPanel responding to ActionButton
			#) # close sidepanel
			) # close FluidRow
		) #close tabItem tcga

	) #close tabItems
) #close dashboardbody

# Put them together into a dashboardPage
dashboardPage(
	skin = "green",
	dashboardHeader(title="MTB Report",titleWidth = "200px"),
	sidebar,
	body
)
