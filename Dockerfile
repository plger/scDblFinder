FROM bioconductor/bioconductor_docker:devel
WORKDIR /home/rstudio/package
COPY . .
RUN R -e "install.packages('remotes'); options(repos = BiocManager::repositories()); \
          remotes::install_local('.', dependencies = TRUE, upgrade = 'always')"
