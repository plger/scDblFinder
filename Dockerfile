FROM bioconductor/bioconductor_docker:devel

MAINTAINER pl.germain@gmail.com

WORKDIR /home/build/package

COPY . /home/build/package 

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN Rscript -e "devtools::install('.', dependencies=TRUE, repos = BiocManager::repositories(), build_vignettes = TRUE)"
