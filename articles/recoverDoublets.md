# Recovering intra-sample doublets

## tl;dr

See the relevant section of the [OSCA
book](https://osca.bioconductor.org/doublet-detection.html#doublet-detection-in-multiplexed-experiments)
for an example of the
[`recoverDoublets()`](https://plger.github.io/scDblFinder/reference/recoverDoublets.md)
function in action on real data. A toy example is also provided in
[`?recoverDoublets`](https://plger.github.io/scDblFinder/reference/recoverDoublets.md).

## Mathematical background

Consider any two cell states $`C_1`$ and $`C_2`$ forming a doublet
population $`D_{12}`$. We will focus on the relative frequency of
inter-sample to intra-sample doublets in $`D_{12}`$. Given a vector
$`\vec p_X`$ containing the proportion of cells from each sample in
state $`X`$, and assuming that doublets form randomly between pairs of
samples, the expected proportion of intra-sample doublets in $`D_{12}`$
is $`\vec p_{C_1} \cdot \vec p_{C_2}`$. Subtracting this from 1 gives us
the expected proportion of inter-sample doublets $`q_{D_{12}}`$.
Similarly, the expected proportion of inter-sample doublets in $`C_1`$
is just $`q_{C_1} =1 - \| \vec p_{C_1} \|_2^2`$.

Now, let’s consider the observed proportion of events $`r_X`$ in each
state $`X`$ that are known doublets. We have $`r_{D_{12}} = q_{D_{12}}`$
as there are no other events in $`D_{12}`$ beyond actual doublets. On
the other hand, we expect that $`r_{C_1} \ll q_{C_1}`$ due to presence
of a large majority of non-doublet cells in $`C_1`$ (same for $`C_2`$).
If we assume that $`q_{D_{12}} \ge q_{C_1}`$ and $`q_{C_2}`$, the
observed proportion $`r_{D_{12}}`$ should be larger than $`r_{C_1}`$ and
$`r_{C_2}`$. (The last assumption is not always true but the $`\ll`$
should give us enough wiggle room to be robust to violations.)

The above reasoning motivates the use of the proportion of known doublet
neighbors as a “doublet score” to identify events that are most likely
to be themselves doublets.
[`recoverDoublets()`](https://plger.github.io/scDblFinder/reference/recoverDoublets.md)
computes the proportion of known doublet neighbors for each cell by
performing a $`k`$-nearest neighbor search against all other cells in
the dataset. It is then straightforward to calculate the proportion of
neighboring cells that are marked as known doublets, representing our
estimate of $`r_X`$ for each cell.

## Obtaining explicit calls

While the proportions are informative, there comes a time when we need
to convert these into explicit doublet calls. This is achieved with
$`\vec S`$, the vector of the proportion of cells from each sample
across the entire dataset (i.e., `samples`). We assume that all cell
states contributing to doublet states have proportion vectors equal to
$`\vec S`$, such that the expected proportion of doublets that occur
between cells from the same sample is $`\| \vec S\|_2^2`$. We then solve

``` math
\frac{N_{intra}}{(N_{intra} + N_{inter}} = \| \vec S\|_2^2
```

for $`N_{intra}`$, where $`N_{inter}`$ is the number of observed
inter-sample doublets. The top $`N_{intra}`$ events with the highest
scores (and, obviously, are not already inter-sample doublets) are
marked as putative intra-sample doublets.

## Discussion

The rate and manner of doublet formation is (mostly) irrelevant as we
condition on the number of events in $`D_{12}`$. This means that we do
not have to make any assumptions about the relative likelihood of
doublets forming between pairs of cell types, especially when cell types
have different levels of “stickiness” (or worse, stick specifically to
certain other cell types). Such convenience is only possible because of
the known doublet calls that allow us to focus on the inter- to
intra-sample ratio.

The most problematic assumption is that required to obtain $`N_{intra}`$
from $`\vec S`$. Obtaining a better estimate would require, at least,
the knowledge of the two parent states for each doublet population. This
can be determined with some simulation-based heuristics but it is likely
to be more trouble than it is worth.

In this theoretical framework, we can easily spot a case where our
method fails. If both $`C_1`$ and $`C_2`$ are unique to a given sample,
all events in $`D_{12}`$ will be intra-sample doublets. This means that
no events in $`D_{12}`$ will ever be detected as inter-sample doublets,
which precludes their detection as intra-sample doublets by
`recoverDoublets`. The computational remedy is to augment the
predictions with simulation-based methods (e.g.,
[`scDblFinder()`](https://plger.github.io/scDblFinder/reference/scDblFinder.md))
while the experimental remedy is to ensure that multiplexed samples
include technical or biological replicates.

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
    ## [1] BiocStyle_2.41.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
    ##  [4] bookdown_0.46       fastmap_1.2.0       xfun_0.58          
    ##  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
    ## [10] rmarkdown_2.31      lifecycle_1.0.5     cli_3.6.6          
    ## [13] sass_0.4.10         pkgdown_2.2.0       textshaping_1.0.5  
    ## [16] jquerylib_0.1.4     systemfonts_1.3.2   compiler_4.6.0     
    ## [19] tools_4.6.0         ragg_1.5.2          bslib_0.11.0       
    ## [22] evaluate_1.0.5      yaml_2.3.12         BiocManager_1.30.27
    ## [25] otel_0.2.0          jsonlite_2.0.0      rlang_1.2.0        
    ## [28] fs_2.1.0            htmlwidgets_1.6.4
