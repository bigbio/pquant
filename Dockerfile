# get shiny server plus tidyverse packages image
FROM rocker/shiny-verse:latest

maintainer Douer douerww@gmail.com

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-openssl-dev \
    libnode-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    curl \
    libxml2-dev \
    libudunits2-0 \
    libudunits2-dev \
    xtail \
    wget \
    gdebi-core \
	vim

# Install library
RUN install2.r --error --skipinstalled \
	--repos 'http://cran.rstudio.com' \
    shinydashboard DT pheatmap shinyjs shinyWidgets IDPmisc devtools htmlwidgets log4r gplots RcppArmadillo BiocManager reshape2 checkmate lme4 ggrepel rhandsontable viridis proteus downloader igraph

RUN R -e "devtools::install_github('Vitek-Lab/MSstats', ref = 'hotfix-fractions-check')"
RUN R -e "if (!library(MSstats, logical.return=T)) quit(status=10)"

RUN R -e "devtools::install_github('Vitek-Lab/MSstatsConvert', ref = 'hotfix-techreplicate')"
RUN R -e "if (!library(MSstats, logical.return=T)) quit(status=10)"


RUN wget https://github.com/Vitek-Lab/MSstatsTMT/archive/refs/heads/master.zip -O /tmp/MSstatsTMT.zip
RUN unzip /tmp/MSstatsTMT.zip
RUN R -e "devtools::install_local('MSstatsTMT-master')"
RUN rm /tmp/MSstatsTMT.zip
RUN R -e "if (!library(MSstatsTMT, logical.return=T)) devtools::install_github('Vitek-Lab/MSstatsTMT', repos='http://cran.rstudio.com')"
RUN R -e "if (!library(MSstatsTMT, logical.return=T)) quit(status=10)"

RUN R -e "BiocManager::install('AnnotationDbi')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('org.Sc.sgd.db')"
RUN R -e "BiocManager::install('org.Rn.eg.db')"
RUN R -e "BiocManager::install('org.EcK12.eg.db')"
RUN install2.r --error --skipinstalled -r http://bioconductor.org/packages/3.0/bioc --deps TRUE \
    AnnotationDbi \
    org.Hs.eg.db \
    org.Sc.sgd.db \
    org.Rn.eg.db \
    && rm -rf /tmp/downloaded_packages/

RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('DOSE')"
RUN R -e "BiocManager::install('GOSemSim')"

RUN R -e "devtools::install_github('YuLab-SMU/enrichplot', repos='http://cran.rstudio.com')"
RUN R -e "if (!library(enrichplot, logical.return=T)) quit(status=10)"

RUN wget https://github.com/YuLab-SMU/clusterProfiler/archive/refs/heads/master.zip -O /tmp/clusterProfiler.zip
RUN unzip /tmp/clusterProfiler.zip
RUN R -e "devtools::install_local('clusterProfiler-master')"
RUN rm /tmp/clusterProfiler.zip
RUN R -e "if (!library(clusterProfiler, logical.return=T)) devtools::install_github('YuLab-SMU/clusterProfiler', repos='http://cran.rstudio.com')"
RUN R -e "if (!library(clusterProfiler, logical.return=T)) quit(status=10)"

RUN installGithub.r bartongroup/Proteus \
    && rm -rf /tmp/downloaded_packages/
RUN R -e "if (!library(proteus, logical.return=T)) devtools::install_github('bartongroup/Proteus')"
RUN R -e "if (!library(proteus, logical.return=T)) quit(status=10)"

RUN install2.r --error --skipinstalled \
	--repos 'http://cran.rstudio.com' \
    R.utils

RUN R -e "BiocManager::install('UniProt.ws')"
RUN R -e "if (!library(UniProt.ws, logical.return=T)) quit(status=10)"

### --------------------------------------

# copy shiny-server.sh to image
COPY shiny-server.sh /usr/bin/

# copy shiny server config to image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

# copy the contents of app folder to image
COPY ./pquantr /srv/shiny-server/pquantr

# select port
EXPOSE 3838

# allow permission for user ‘shiny’ to run
RUN sudo chown -R shiny:shiny /srv/shiny-server

# install linux programs to enable conversion of ms dos file to unix file
RUN apt-get update && apt-get install -y dos2unix

# we do this so that the shiny-server.sh file is recognized by the linux machine
RUN dos2unix /usr/bin/shiny-server.sh && apt-get --purge remove -y dos2unix && rm -rf /var/lib/apt/lists/*

# Change access permissions to shiny-server.sh - did not need this for my purposes
RUN chmod -R 755 /usr/bin/shiny-server.sh

# run app
CMD ["/usr/bin/shiny-server.sh"]
