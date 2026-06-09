FROM bioconductor/bioconductor_docker:devel
WORKDIR /home/rstudio/package
RUN sudo apt-get install libmagick++-dev
COPY . .
RUN R -e "install.packages(c('remotes','magick')); options(repos = BiocManager::repositories()); \
          remotes::install_local('.', dependencies = TRUE, upgrade = 'always')"
