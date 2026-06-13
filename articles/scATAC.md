# Doublet identifiation in single-cell ATAC-seq

Abstract

An introduction to the methods implemented for doublet detection in
single-cell ATAC-seq.

## Introduction

Analyses in single-cell RNAseq are typically limited to a relatively
small (e.g. one or two thousands) set of features that are most
informative; these are often the genes with a higher expression (and
hence more chances of being quantified). In contrast, single-cell
ATACseq (scATACseq) data is considerably more sparse, with most reads
being spread across hundreds of thousands of regions. In this context,
selecting a subset of genes is highly ineffective, and therefore many of
the methods developed for single-cell RNAseq are not easily applicable,
and need to be adapted. Methods have therefore been developed
specifically for scATACseq data (Granja et al. 2021; Thibodeau et
al. 2021).

This vignette presents different approaches to doublet detection in
single-cell ATAC-seq implemented in the package: the first is an
adaptation of `scDblFinder`, the second a reimplementation of the AMULET
method from Thibodeau et al. (2021). The latter has the advantage of
capturing homotypic doublets, but does not perform well in all datasets,
and especially requires the cells to have a high library size. We
therefore next present two ways of combining the two.

## Applying the scDblFinder method

With default parameters, the `scDblFinder` method performs very poorly
on scATACseq data due to the increase spread of the reads across many
features. Since working with all features (i.e. tiles or peaks) is
computationally very expensive, an alternative approach is to begin by
reducing the size of the dataset, not through *selection* (as in
scRNAseq), but by *aggregating* correlated features into a relatively
small set. This has the advantage of using all information, as well as
making the count data more continuous. This method yields comparable
performance to specialized single-cell ATACseq software (Germain et al.,
2021).

The feature aggregation can be triggered using the
`aggregateFeatures=TRUE` argument, which will aggregate peak or tile
counts into the number of meta-features defined by the `nfeatures`. If
the number of meta-features is low (which we recommend), the
meta-features can be directly used to calculated distances rather than
going through the SVD step (which can be triggered with the `processing`
argument). Such an example would be:

``` r

suppressPackageStartupMessages(library(scDblFinder))
# we use a dummy SingleCellExperiment as example:
sce <- mockDoubletSCE(ngenes=300)
# setting low number of artificial doublets (same as ncells) just for speedup:
sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
```

    ## Aggregating features...

    ## Creating ~526 artificial doublets...

    ## Dimensional reduction

    ## Evaluating kNN...

    ## Training model...

    ## iter=0, 17 cells excluded from training.

    ## iter=1, 36 cells excluded from training.

    ## iter=2, 42 cells excluded from training.

    ## Threshold found:0.764

    ## 23 (4.4%) doublets called

If you encounter problems running the aggregation-based approach on
large datasets, first make sure you have the `mbkmeans` package
installed.

## Using the Amulet method

The AMULET method from Thibodeau et al. (2021) is based on the
assumption that, in a diploid cell, any given genomic region should be
captured at most twice. Therefore, cells with loci covered by more than
two fragments are indicative of the droplet being a doublet. Of note,
this approach has the advantage of capturing homotypic doublets, which
instead tend to be missed by other methods. Since it was only available
in the form of a mixture of java and python scripts, we re-implemented
the method in `scDblFinder` (see
[`?amulet`](https://plger.github.io/scDblFinder/reference/amulet.md)),
leading to equal or superior results to the original implementation
(Germain et al. 2021).

As in the original implementation, we recommend excluding the
mitochondrial and sex chromosomes, as well as repetitive regions. This
can be specified with the `regionsToExclude` argument (see the
underlying
[`?getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md)).
It can be used as follows:

``` r

# here we use a dummy fragment file for example:
fragfile <- system.file("extdata", "example_fragments.tsv.gz", package="scDblFinder")

# we might also give a GRanges of repeat elements, so that these regions are excluded:
suppressPackageStartupMessages(library(GenomicRanges))
repeats <- GRanges("chr6", IRanges(1000,2000))
# it's better to combine these with mitochondrial and sex chromosomes
otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
# here since I don't know what chromosome notation you'll be using I've just put them all,
# although this will trigger a warning when combining them:
toExclude <- suppressWarnings(c(repeats, otherChroms))
# we then launch the method
res <- amulet(fragfile, regionsToExclude=toExclude)
```

    ## Fragment file is not tabix-indexed, requiring thewhole file to be imported in memory.

    ## 02:51:43 PM - Splitting and subsetting barcodes...

    ## 02:51:43 PM - Obtaining overlaps...

``` r

res
```

    ##          nFrags uniqFrags nAbove2 total.nAbove2     p.value     q.value
    ## barcode1    878       878       1             1 0.475069053 0.791781755
    ## barcode2   2401      2401       0             0 0.798103482 0.798103482
    ## barcode3   2325      2325       1             1 0.475069053 0.791781755
    ## barcode4   1882      1882       0             0 0.798103482 0.798103482
    ## barcode5   1355      1355       6             6 0.001335761 0.006678806

The results is a data.frame with statistics for each barcode, including
a p-value. In contrast to the `scDblFinder` score, a lower p-value here
is indicative of the droplet being more likely to be a doublet (as in
the original method). By default, only the barcodes with a minimum
number of reads are considered, but it is possible to specify the
droplets for which to gather statistics using the `barcodes` argument.

While the package includes an implementation that works based on
peak/tile count matrices (see
[`?amuletFromCounts`](https://plger.github.io/scDblFinder/reference/amuletFromCounts.md)),
it has a much lower performance with respect to the one based directly
on the fragment files (see
[`?amulet`](https://plger.github.io/scDblFinder/reference/amulet.md)),
and we therefore discourage its use.

The workhorse behind the `amulet` function is the `getFragmentOverlaps`,
which also includes all of the relevant arguments. If the fragment files
are not Tabix-indexed, the whole fragment file will have to be loaded in
memory for processing; while this ensures relatively rapid computation,
it has high memory requirements. Therefore, if the fragment file is
Tabix-indexed (as is for instance done as part of the ArchR pipeline),
it will be read and processed per chromosome, which is a little slower
due to overhead, but keeps memory requirements rather low. This behavior
can be disabled by specifying `fullInMemory=TRUE`.

## Combining mehtods

While the `scDblFinder`-based approach generally performs well, none of
the two approach is optimal across all datasets tested. We therefore
investigated two strategies for combining the rationales of each
approach.

The Amulet method tends to perform best with datasets that have
homotypic doublets and where cells have a high library size (i.e. median
library size per cell of 10-15k reads), while the `scDblFinder`-based
approach works better for heterotypic doublets. Until an optimal
solution is found, we recommend using multiple approaches to inform
decisions, in particular using the p-value combination method below.

### The Clamulet method

The `clamulet` method (Classification-powered Amulet-like method)
operates similarly to the `scDblFinder` method, but generates artificial
doublets by operating on the fragment coverages. This has the advantage
that the number of loci covered by more than two reads can be computed
for artificial doublets, enabling the use of this feature (along with
the kNN-based ones) in a classification scheme. It however has the
disadvantage of being rather slow and memory hungry, and appears to be
outperformed by a simple p-value combination of the two methods (see
below). We therefore *do not* recommend its usage.

The `clamulet` method uses the aforementioned aggregation approach, and
its usage includes a number of arguments from both the `scDblFinder` and
`amulet` method (see in particular
[`?getFragmentOverlaps`](https://plger.github.io/scDblFinder/reference/getFragmentOverlaps.md)):

``` r

# not run
d <- clamulet("path/to/fragments.tsv.gz")
```

Since our dummy fragment file is so small (5 barcodes), here we’ll have
to adjust the arguments for an example to run:

``` r

d <- clamulet(fragfile, k=2, nfeatures=3)
```

    ## 02:51:44 PM - Reading full fragments...

    ## 02:51:44 PM - Splitting and subsetting barcodes...

    ## 02:51:44 PM - Computing coverages

    ## 02:51:45 PM - Obtaining windows

    ## 02:51:45 PM - Obtaining window counts

    ## 02:51:45 PM - Aggregating features

    ## Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth =
    ## TRUE, : You're computing too large a percentage of total singular values, use a
    ## standard svd instead.

    ## 02:51:45 PM - Computing features for artificial doublets

    ## 02:51:45 PM - Counting overlaps for real cells

    ## 02:51:45 PM - Counting overlaps for artificial doublets

    ## 02:51:45 PM - Scoring network

    ## 02:51:45 PM - Iterative training

    ## iter=0, 5 cells excluded from training.

    ## iter=1, 5 cells excluded from training.

    ## 02:51:45 PM Done!

``` r

d
```

    ##          total nAbove2 total.nAbove2  weighted ratio.k2 include.in.training
    ## barcode1    19       1             1 1.0000000      1.0               FALSE
    ## barcode2     8       0             0       NaN      0.5               FALSE
    ## barcode3     7       1             1 0.4373488      0.5               FALSE
    ## barcode4     8       0             0       NaN      1.0               FALSE
    ## barcode5    14       6             6 0.4956057      0.5               FALSE
    ##              score
    ## barcode1 0.9999989
    ## barcode2 0.9999989
    ## barcode3 0.9999989
    ## barcode4 0.9999989
    ## barcode5 0.9999989

The score can then be interpreted as for `scDblFinder`. We however note
that this method proved *inferior to alternatives*.

### Simple p-value combination

The amulet and scDblFinder scores above can be simply combined by
treating them as p-values and aggregating them (here using Fisher’s
method from the
*[aggregation](https://CRAN.R-project.org/package=aggregation)* package,
but see also the *[metap](https://CRAN.R-project.org/package=metap)*
package):

``` r

res$scDblFinder.p <- 1-colData(sce)[row.names(res), "scDblFinder.score"]
res$combined <- apply(res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
  x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
  suppressWarnings(aggregation::fisher(x))
})
```

We found this to perform better than averaging the scores or their
ranks, and while it is not the very best method in any of the datasets
tested, it has a more robust performance overall (see Germain et al.,
2021).

## References

Jeffrey M. Granja et al., “ArchR Is a Scalable Software Package for
Integrative Single-Cell Chromatin Accessibility Analysis,” Nature
Genetics, February 25, 2021, 1–9,
<https://doi.org/10.1038/s41588-021-00790-6>

Asa Thibodeau et al., “AMULET: A Novel Read Count-Based Method for
Effective Multiplet Detection from Single Nucleus ATAC-Seq Data,” Genome
Biology 22, no. 1 (December 2021): 252,
<https://doi.org/10.1186/s13059-021-02469-x>

Pierre-Luc Germain et al., “Doublet Identification in Single-Cell
Sequencing Data Using ScDblFinder” (F1000Research, September 28, 2021),
<https://doi.org/10.12688/f1000research.73600.1>

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
    ##  [1] scDblFinder_1.27.4          SingleCellExperiment_1.35.1
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
