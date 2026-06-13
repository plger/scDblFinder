FROM bioconductor/bioconductor_docker:devel
WORKDIR /home/rstudio/package
RUN sudo apt-get install libmagick++-dev
COPY . .
RUN R -e "install.packages(c('remotes','magick')); options(repos = BiocManager::repositories()); \
          remotes::install_local('.', dependencies = TRUE, upgrade = 'always')"
RUN R -e "options(ExperimentHub.ask=FALSE, AnnotationHub.ask=FALSE); sce <- scRNAseq::BachMammaryData(samples=c('G_1','G_2'))"
