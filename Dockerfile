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
RUN R -e "install.packages(c('shinydashboard','DT','pheatmap','shinyjs','shinyWidgets','IDPmisc','devtools','htmlwidgets','log4r','gplots','RcppArmadillo','BiocManager','reshape2','checkmate','lme4','ggrepel','rhandsontable'),dependencies=TRUE,repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstats')"
RUN R -e "BiocManager::install('AnnotationDbi')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('org.Sc.sgd.db')"
RUN R -e "BiocManager::install('org.Rn.eg.db')"


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
