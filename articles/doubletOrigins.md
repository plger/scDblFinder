# Doublet origins

Abstract

Functions to identify the origins of doublets, i.e. the cell types
composing them, and perform analysis on these origins.

## Introduction

This vignette addresses the question of doublet origins,
i.e. identifying which combination of cell-types or clusters generated a
given doublet. The first part of the vignette concerns the identifying
doublet origins, while the second part covers test for enrichment in
specific kinds of doublets.

## Identifying doublet origins

As a preliminary, note that *doublet origins called in versions prior to
1.27.8 were mostly wrong*.

### kNN-based approach

Identifying doublet origins requires running
[`scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.md)
on your data passing your cell type labels, i.e. using the `clusters`
argument of
[`scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.md)
to pass cluster labels (the cluster labels of doublets will be ignored,
so no need to worry about those). For example:

``` r

# we first generate example data:
sce <- mockDoubletSCE(ncells=c(100,150,200), ngenes=500, dbl.rate=0.15)
# we call doublets (somewhat quickly):
sce <- scDblFinder(sce, artificialDoublets = 300, clusters="cluster", threshold=0.5)
```

    ## 3 clusters

    ## Creating ~300 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 16 cells excluded from training.

    ## iter=1, 8 cells excluded from training.

    ## iter=2, 7 cells excluded from training.

    ## Threshold found:0.5

    ## 11 (2.2%) doublets called

Ran in this fashion, the output already includes a rough kNN-based guess
of the doublet origin:

``` r

table(call=sce$scDblFinder.class,
      mostLikelyOrigin=sce$scDblFinder.mostLikelyOrigin)
```

    ##          mostLikelyOrigin
    ## call      cluster1+cluster2 cluster1+cluster3 cluster2+cluster3
    ##   singlet                16                15                13
    ##   doublet                 5                 4                 2

Note that since this is independent of the classifier (and based purely
on the kNN), some non-doublets may have a combination of clusters as
their `mostLikelyOrigin`, we recommend ignoring those.
(`scDblFinder.originAmbiguous` in colData also contains information
about whether the call is ambiguous or not)

We can compare this to the ground truth (focusing on actual doublets):

``` r

isDoublet <- which(sce$type=="doublet")
table(prediction=sce$scDblFinder.mostLikelyOrigin[isDoublet],
      truth=sce$origin[isDoublet])
```

    ##                    truth
    ## prediction          cluster1+cluster2 cluster1+cluster3 cluster2+cluster3
    ##   cluster1+cluster2                17                 0                 0
    ##   cluster1+cluster3                 0                14                 0
    ##   cluster2+cluster3                 0                 0                15

In this easy case, all true doublet origins are accurately identified.
However, a slightly more powerful approach is to identify doublets using
a specific multi-label classifier.

### Multi-label classifier

The
[`identifyDoubletOrigins()`](https://plger.github.io/scDblFinder/reference/identifyDoubletOrigins.md)
function trains a multi-label classifier on doublet origins from
artificial doublets, which can then be used to predict doublet origins.

``` r

# to keep things quick here, with cap it to 20 rounds: 
clf <- identifyDoubletOrigins(sce, clusters="cluster", max_rounds = 20, verbose = FALSE)
```

Here the function uses the `sce$scDblFinder.class` doublets, but one
could specify manually for which cells predictions should be made.

We can compare the predictions to the ground truth:

``` r

table(prediction=clf$calls, truth=colData(sce)[names(clf$calls), "origin"])
```

    ##                    truth
    ## prediction          cluster1+cluster2 cluster1+cluster3 cluster2+cluster3
    ##   cluster1+cluster2                 5                 0                 0
    ##   cluster1+cluster3                 0                 4                 0
    ##   cluster2+cluster3                 0                 0                 2

Note that in real data, due to the very large variations in library
sizes, it is not uncommon for doublets to contain a large fraction of
reads from one cell type, and conversely a small fraction from the other
cell type. In those circumstances, while one of the two originating cell
type is clear, the other can however be mis-assigned.

## Doublet enrichment analysis

Doublet enrichment analysis aims to investigate whether certain doublet
types are found more often than would be expected by chance. Such
enrichment could, for instance, indicate cell-to-cell interactions.

When ran with the `clusters` argument,
[`scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.md)
saves some stats on the doublet origins in
`metadata(sce)$scDblFinder.stats`. Similarly, the output of
[`identifyDoubletOrigins()`](https://plger.github.io/scDblFinder/reference/identifyDoubletOrigins.md)
contains a `stats` slot:

``` r

clf$stats
```

    ##         combination observed expected deviation prop.deviation        FNR
    ## 1 cluster1+cluster2        5  0.34584   4.65416       3.599394 0.00000000
    ## 2 cluster1+cluster3        4  0.41920   3.58080       2.769288 0.01075269
    ## 3 cluster2+cluster3        2  0.52800   1.47200       1.138403 0.00000000
    ##   difficulty
    ## 1 0.05194089
    ## 2 0.05148621
    ## 3 0.04395504

This table contains the observed and expected types of doublets, along
with the difficulty of their identification. The
[`plotDoubletMap()`](https://plger.github.io/scDblFinder/reference/plotDoubletMap.md)
function can be used to visualize enrichment. Most importantly, these
compiled statistics can then be used to test for enrichment.

Specifically, we define two forms of doublet enrichment, each with their
own testing function:

- The
  [`clusterStickiness()`](https://plger.github.io/scDblFinder/reference/clusterStickiness.md)
  function tests whether each cluster forms more doublet than would be
  expected given its abundance, by default using a single quasi-binomial
  model fitted across all doublet types. (Only applicable with \>=4
  clusters)
- The
  [`doubletPairwiseEnrichment()`](https://plger.github.io/scDblFinder/reference/doubletPairwiseEnrichment.md)
  function separately tests whether each specific doublet type
  (i.e. combination of clusters) is more abundant than expected, by
  default using a poisson model.

For example:

``` r

doubletPairwiseEnrichment(clf$stats)
```

    ## theta=0.188806055483142

    ##         combination  log2enrich    p.value       FDR
    ## 1 cluster1+cluster2  0.55183959 0.09357226 0.2807168
    ## 2 cluster1+cluster3  0.06689129 0.32700281 0.6540056
    ## 3 cluster2+cluster3 -1.08032686 1.00000000 1.0000000

This shows no significant enrichment in any of the doublet types, which
makes sense because they were randomly generated.

For more detail, see the appropriate section of [Germain et al.,
2022](https://f1000research.com/articles/10-979).

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] scDblFinder_1.27.9          SingleCellExperiment_1.35.1
    ##  [3] SummarizedExperiment_1.43.0 Biobase_2.73.1             
    ##  [5] GenomicRanges_1.65.0        Seqinfo_1.3.0              
    ##  [7] IRanges_2.47.2              S4Vectors_0.51.3           
    ##  [9] BiocGenerics_0.59.7         generics_0.1.4             
    ## [11] MatrixGenerics_1.25.0       matrixStats_1.5.0          
    ## [13] BiocStyle_2.41.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1         viridisLite_0.4.3        vipor_0.4.7             
    ##  [4] dplyr_1.2.1              farver_2.1.2             viridis_0.6.5           
    ##  [7] S7_0.2.2                 Biostrings_2.81.3        bitops_1.0-9            
    ## [10] fastmap_1.2.0            RCurl_1.98-1.19          scrapper_1.7.3          
    ## [13] bluster_1.23.0           GenomicAlignments_1.49.0 XML_3.99-0.23           
    ## [16] digest_0.6.39            rsvd_1.0.5               lifecycle_1.0.5         
    ## [19] cluster_2.1.8.2          magrittr_2.0.5           compiler_4.6.0          
    ## [22] rlang_1.2.0              sass_0.4.10              tools_4.6.0             
    ## [25] igraph_2.3.2             yaml_2.3.12              data.table_1.18.4       
    ## [28] rtracklayer_1.73.0       knitr_1.51               S4Arrays_1.13.0         
    ## [31] htmlwidgets_1.6.4        xgboost_3.2.1.1          curl_7.1.0              
    ## [34] DelayedArray_0.39.3      RColorBrewer_1.1-3       abind_1.4-8             
    ## [37] BiocParallel_1.47.0      desc_1.4.3               grid_4.6.0              
    ## [40] beachmat_2.29.0          ggplot2_4.0.3            scales_1.4.0            
    ## [43] MASS_7.3-65              cli_3.6.6                rmarkdown_2.31          
    ## [46] crayon_1.5.3             ragg_1.5.2               otel_0.2.0              
    ## [49] httr_1.4.8               rjson_0.2.23             BiocBaseUtils_1.15.1    
    ## [52] scuttle_1.23.1           ggbeeswarm_0.7.3         cachem_1.1.0            
    ## [55] parallel_4.6.0           BiocManager_1.30.27      XVector_0.53.0          
    ## [58] restfulr_0.0.17          vctrs_0.7.3              Matrix_1.7-5            
    ## [61] jsonlite_2.0.0           bookdown_0.46            BiocSingular_1.29.0     
    ## [64] BiocNeighbors_2.7.2      ggrepel_0.9.8            beeswarm_0.4.0          
    ## [67] irlba_2.3.7              scater_1.41.1            systemfonts_1.3.2       
    ## [70] jquerylib_0.1.4          glue_1.8.1               pkgdown_2.2.0           
    ## [73] codetools_0.2-20         gtable_0.3.6             GenomeInfoDb_1.49.1     
    ## [76] BiocIO_1.23.3            UCSC.utils_1.9.0         ScaledMatrix_1.21.0     
    ## [79] tibble_3.3.1             pillar_1.11.1            htmltools_0.5.9         
    ## [82] R6_2.6.1                 textshaping_1.0.5        evaluate_1.0.5          
    ## [85] lattice_0.22-9           Rsamtools_2.29.0         cigarillo_1.3.0         
    ## [88] bslib_0.11.0             Rcpp_1.1.1-1.1           gridExtra_2.3           
    ## [91] SparseArray_1.13.2       xfun_0.58                fs_2.1.0                
    ## [94] pkgconfig_2.0.3
