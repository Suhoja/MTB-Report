FROM		rocker/shiny:4.2.1

RUN             apt-get update
RUN             apt-get upgrade -y
ARG DEBIAN_FRONTEND=noninteractive
RUN             apt-get install -fy --no-install-recommends texlive-base texlive-fonts-recommended texlive-latex-recommended texlive-latex-extra texinfo
RUN             apt-get update
RUN             apt-get install -y -f default-jdk openssl curl xml2 libcurl4-openssl-dev libssl-dev libxml2-dev file dirmngr gnupg2 apt-transport-https ca-certificates software-properties-common

RUN		apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN		add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN		apt-get update
RUN		apt-get install -y r-base

RUN		Rscript -e "install.packages(c('DT', 'ggplot2', 'xtable', 'timeDate','timeSeries', 'pander', 'knitr', 'RCurl', 'httr', 'openssl', 'curl', 'XML', 'covr', 'stringr', 'stringi', 'plyr','dplyr', 'BiocManager', 'parallel', 'RestRserve', 'stats4', 'testthat', 'jsonlite', 'tools','devtools','gdtools','xml2','magrittr','maftools','latexpdf'), dependencies=TRUE)"
RUN		Rscript -e "BiocManager::install(c('shinydashboard','S4Vectors', 'AnnotationDbi', 'XVector', 'Biostrings', 'Biobase', 'GenomeInfoDb', 'VariantAnnotation', 'BiocGenerics', 'EnsDb.Hsapiens.v86', 'IRanges', 'GenomicRanges', 'GenomicFeatures', 'gwascat', 'biomaRt', 'liftOver', 'EnsDb.Hsapiens.v75','maftools','org.Hs.eg.db'))"
RUN		Rscript -e "devtools::install_github('mariodeng/FirebrowseR')"
RUN		Rscript -e "install.packages(c('shinyjs','openxlsx'))"
RUN		Rscript -e "options(timeout=1000); BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19','BSgenome.Hsapiens.UCSC.hg38'))"

VOLUME		/work
RUN		mkdir files

WORKDIR		/work
RUN		apt-get install -y wget lsb-release
COPY		install_shiny.sh /work/install_shiny.sh
RUN		chmod 755 /work/install_shiny.sh
RUN		/work/install_shiny.sh

COPY		src/ /opt/MTB/src/
COPY		data/ /opt/MTB/data/
RUN		chmod -R 777 /opt/MTB/data
COPY		helpers/ /opt/MTB/helpers/
COPY		tests/ /opt/MTB/tests/
COPY		example_data/ /opt/MTB/example_data/
RUN		chmod -R 777 /opt/MTB/example_data
COPY		shiny-server.conf /etc/shiny-server/

RUN		chmod -R 777 /opt/MTB/data/

RUN		chmod +x /opt/MTB/src/starter_s_c.sh
WORKDIR		/opt/MTB

RUN 		Rscript -e 'install.packages(c("vctrs"))'
RUN 		Rscript -e "BiocManager::install('GenomicFeatures')"

RUN		Rscript src/txDBdata_preprocess.r

RUN		Rscript src/preprocess.R

RUN		chmod -R 777 /opt/MTB/src/tmp
RUN		chown -R shiny:shiny /var/lib/shiny-server
RUN		chmod -R 755 /var/lib/shiny-server
RUN		su shiny

RUN     Rscript "/opt/MTB/src/install_packages.r"

EXPOSE 3838

ENTRYPOINT	["/opt/MTB/src/starter_s_c.sh"]
