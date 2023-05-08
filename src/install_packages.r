mtbreport_packages <- c('DT', 'ggplot2', 'xtable', 'timeDate','timeSeries', 'pander', 'knitr', 'RCurl', 'httr', 'openssl', 'curl', 'XML', 'covr', 'stringr', 'stringi', 'plyr','dplyr', 'BiocManager', 'parallel', 'RestRserve', 'stats4', 'testthat', 'jsonlite', 'tools','devtools','gdtools','xml2','shinyjs','openxlsx','maftools')
mtbreport_bm_packages <- c('shinydashboard','S4Vectors', 'AnnotationDbi', 'XVector', 'Biostrings', 'Biobase', 'GenomeInfoDb', 'VariantAnnotation', 'BiocGenerics', 'EnsDb.Hsapiens.v86', 'IRanges', 'GenomicRanges', 'GenomicFeatures', 'gwascat', 'biomaRt', 'liftOver', 'EnsDb.Hsapiens.v75','maftools','EnsDb.Hsapiens.v75','gwascat','liftOver','AnnotationHub','interactiveDisplayBase')

not_installed <- mtbreport_packages[!(mtbreport_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed)) install.packages(not_installed, repos='http://cran.us.r-project.org', dependencies=TRUE)

not_installed_bm <- mtbreport_bm_packages[!(mtbreport_bm_packages %in% installed.packages()[ , "Package"])]
if(length(not_installed_bm)) BiocManager::install(not_installed_bm)

library(devtools)
install_github('mariodeng/FirebrowseR')