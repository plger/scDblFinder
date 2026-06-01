# Introduction to the scDblFinder package

Abstract

An introduction to the various methods included in the scDblFinder
package.

## Introduction

The `scDblFinder` package gathers various methods for the detection and
handling of doublets/multiplets in single-cell sequencing data
(i.e. multiple cells captured within the same droplet or reaction
volume). This vignette provides a brief overview of the different
approaches (which are each covered in their own vignettes) for
single-cell RNA sequencing. *For doublet detection in genomic data, see
the [scATACseq
vignette](https://plger.github.io/scDblFinder/articles/scATAC.md)*. For
a more general introduction to the topic of doublets, refer to the [OCSA
book](https://osca.bioconductor.org/doublet-detection.html).

All methods require as an input either a matrix of counts or a
*[SingleCellExperiment](https://bioconductor.org/packages/3.23/SingleCellExperiment)*
containing count data. With the exception of
[findDoubletClusters](https://plger.github.io/scDblFinder/articles/findDoubletClusters.md),
which operates at the level of clusters (and consequently requires
clustering information), all methods try to assign each cell a score
indicating its likelihood (broadly understood) of being a doublet.

The approaches described here are *complementary* to doublets identified
via cell hashes and SNPs in multiplexed samples: while hashing/genotypes
can identify doublets formed by cells of the same type (homotypic
doublets) from two samples, which are often nearly undistinguishable
from real cells transcriptionally (and hence generally unidentifiable
through the present package), it cannot identify doublets made by cells
of the same sample, even if they are heterotypic (formed by different
cell types). Instead, the methods presented here are primarily geared
towards the identification of heterotypic doublets, which for most
purposes are also the most critical ones.

  

### computeDoubletDensity

The `computeDoubletDensity` method (formerly
[`scran::doubletCells`](https://rdrr.io/pkg/scran/man/defunct.html))
generates random artificial doublets from the real cells, and tries to
identify cells whose neighborhood has a high local density of articial
doublets. See
[computeDoubletDensity](https://plger.github.io/scDblFinder/articles/computeDoubletDensity.md)
for more information.

### recoverDoublets

The `recoverDoublets` method is meant to be used when some doublets are
already known, for instance through genotype-based calls or cell hashing
in multiplexed experiments. The function then tries to identify
intra-sample doublets that are neighbors to the known inter-sample
doublets. See
[recoverDoublets](https://plger.github.io/scDblFinder/articles/recoverDoublets.md)
for more information.

### scDblFinder

The `scDblFinder` method combines both known doublets (if available) and
cluster-based artificial doublets to identify doublets. The approach
builds and improves on a variety of earlier efforts, and is at present
the most accurate approach included in this package. See
[scDblFinder](https://plger.github.io/scDblFinder/articles/scDblFinder.md)
for more information.

### directDblClassification

The `directDblClassification` method identifies doublets by training a
classifier directly on gene expression. This follows the same procedure
as `scDblFinder` for doublet generation and iterative training, but
skips the *k*-nearest neighbor step and directly uses the matrix of real
cells and artificial doublets. This is computationally more intensive
and generally leads to worse predictions than `scDblFinder`, and it is
included chiefly for comparative purposes. See
[`?directDblClassification`](https://plger.github.io/scDblFinder/reference/directDblClassification.md)
for more information.

### findDoubletClusters

The `findDoubletClusters` method identifies clusters that are likely to
be composed of doublets by estimating whether their expression profile
lies between two other clusters. See
[findDoubletClusters](https://plger.github.io/scDblFinder/articles/findDoubletClusters.md)
for more information.

  

## Installation

``` r

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scDblFinder")

# or, to get that latest developments:
BiocManager::install("plger/scDblFinder")
```

## Which method to choose?

A benchmark of the main methods available in the package is presented in
the [scDblFinder paper](https://f1000research.com/articles/10-979/).
While the different methods included here have their values, overall the
`scDblFinder` method had the best performance (also superior to other
methods not included in this package), and should be used by default.

## Session information

``` r

sessionInfo()
```

    ## R version 4.6.0 (2026-04-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] BiocStyle_2.40.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
    ##  [4] bookdown_0.46       fastmap_1.2.0       xfun_0.57          
    ##  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
    ## [10] rmarkdown_2.31      lifecycle_1.0.5     cli_3.6.6          
    ## [13] sass_0.4.10         pkgdown_2.2.0       textshaping_1.0.5  
    ## [16] jquerylib_0.1.4     systemfonts_1.3.2   compiler_4.6.0     
    ## [19] tools_4.6.0         ragg_1.5.2          bslib_0.11.0       
    ## [22] evaluate_1.0.5      yaml_2.3.12         BiocManager_1.30.27
    ## [25] otel_0.2.0          jsonlite_2.0.0      rlang_1.2.0        
    ## [28] fs_2.1.0            htmlwidgets_1.6.4
